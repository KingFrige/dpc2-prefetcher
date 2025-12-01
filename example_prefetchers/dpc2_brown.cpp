//
// Sandbox Based Optimal Offset Estimation (SBOOE)
// Data Prefetching Championship 2
// Paper #11
//
// Nathan Brown and Resit Sendag
// Department of Electrical, Computer and Biomedical Engineering
// University of Rhode Island, Kingston, RI, USA
// (nattayb, sendag}@ele.uri.edu
//	
// To compile:
//	g++ -Wall -o dpc2sim prefetcher.cpp lib/dpc2sim.a
// 

#include <stdio.h>
#include <math.h>
#include <list>
extern "C" {
#include "../inc/prefetcher.h"
}


// Prefetch Buffer (PB) size in entries.
#define PREFETCH_BUFFER_SIZE				64

// Sandbox size in entries.
#define L2_ACCESS_COUNT						1024

// Latency averages in entires.
#define LATENCY_AVERAGES					32
#define LATENCY_AVERAGES_SHIFT				(unsigned short)log2(LATENCY_AVERAGES)

// Number of prefetchers (total).
#define PREFETCHER_COUNT					9


//
// Abstract generic prefetcher class definition. Any prefetchers must implement these
// functions because they're pure virtual. It makes a clear interface if we have more
// than one type of prefetcher.
//
class BasePrefetcher
{
  public:
    BasePrefetcher() {}
    virtual ~BasePrefetcher() {}
    virtual unsigned long long int predict(unsigned long long int pc,
        unsigned long long int address) = 0;
    virtual int getOffset() = 0;
};

//
// Sequential offset prefetcher class definition. This prefetcher always predicts
// Address access Ai + offset.
//
class SequentialPrefetcher : public BasePrefetcher
{
  public:
    SequentialPrefetcher(int offset) : BasePrefetcher(),
    offset(offset) {}
    virtual ~SequentialPrefetcher() {}
    int getOffset() {return this->offset;};
    void setOffset(int offset) {this->offset = offset;};
    unsigned long long int predict(unsigned long long int pc,
        unsigned long long int address)
    {
      if(offset > 0)
        return (address + (offset * CACHE_LINE_SIZE));
      else
        return (address - ((offset * -1) * CACHE_LINE_SIZE));
    }
  protected:
    int offset;
};


// Prefetchers and resources (580 bits).
SequentialPrefetcher * prefetchers[PREFETCHER_COUNT];
unsigned long long int predictions[PREFETCHER_COUNT];			// 9 x 64 bits = 576 bits.
int primary;													// 4 bits.

// Resources for calculuting average L2 to lower level memory 
// access latency (898 bits).
std::list<unsigned short> mshrAddresses, 						// 16 bits of address + 1 valid bit each x 16 + 4 bits head index + 
                                                    // 4 bits tail index = 280 bits
  mshrCycles, 							// 16 bits of cycle each x 16 (head/tail index shared with cycles) = 256 bits
  mshrLatencies;						// 10 bits of cycle count each x 32 + 5 bits tail + 5 bits head = 330 bits 
                            // (this could be implemented with less resources if needed).
unsigned short lastL2Access;									// 16 bits.
unsigned short averageLatencyCycles;							// 16 bits.

// Sandbox and resources (240,138 bits bits).
std::list<unsigned short> sandbox[PREFETCHER_COUNT];			// 16 bits of address each x 1024 + 9 + 10 bits for index (shared
                                                          // with cyclesToArrival) x 9 = 147,546 bits.
std::list<unsigned short> cyclesToArrival[PREFETCHER_COUNT];	// 10 bits of cycles each x 1024 x 9 = 92,160 bits.
unsigned short sandboxScoreTotal[PREFETCHER_COUNT],				// 16 bits of score each x 9 = 144 bits.
               lateScoreTotal[PREFETCHER_COUNT],				// 16 bits of score each x 9 = 144 bits.
               usefulScoreTotal[PREFETCHER_COUNT];				// 16 bits of score each x 9 = 144 bits.

// Prefetch buffers (1036 bits).
std::list<unsigned short> prefetchBuffer;						// 16 bits of address each x 64 + 10 bits index.


//
// Prefetch filter function.
//
//	Returns true if a potential prefetch address passed the filtering criteria and
//	false if it did not. It's written as individual tests so debug information can
//	be printed.
//
bool filter(unsigned long long int baseAddress,
    unsigned long long int prefetchAddress)
{
  std::list<unsigned short>::iterator it;

  // Check if prefetch is in MSHR.
  for(it = mshrAddresses.begin(); it != mshrAddresses.end(); ++it)
  {
    if(*it == ((prefetchAddress >> 6) & 0xFFFF))
      return false;
  }
  // Check if prefetch is in prefetch buffer.
  for(it = prefetchBuffer.begin(); it != prefetchBuffer.end(); ++it)
  {
    if(*it == ((prefetchAddress >> 6) & 0xFFFF))
      return false;
  }
  // Check that the prefetch is in the same page as the base address.
  if((prefetchAddress >> 12) != (baseAddress >> 12))
  {
    return false;
  }
  // Check that the read queue isn't too full.
  else if(get_l2_read_queue_occupancy(0) > (L2_READ_QUEUE_SIZE - 6))
  {
    return false;
  }

  // If we got this far return true.
  return true;
}

void l2_prefetcher_initialize(int cpu_num)
{  
  // Instantiate our prefetchers (-1 through +8, 9 in total).
  prefetchers[0] = new SequentialPrefetcher(-1);
  prefetchers[1] = new SequentialPrefetcher(1);
  prefetchers[2] = new SequentialPrefetcher(2);
  prefetchers[3] = new SequentialPrefetcher(3);
  prefetchers[4] = new SequentialPrefetcher(4);
  prefetchers[5] = new SequentialPrefetcher(5);
  prefetchers[6] = new SequentialPrefetcher(6);
  prefetchers[7] = new SequentialPrefetcher(7);
  prefetchers[8] = new SequentialPrefetcher(8);

  // The starting "old" last l2 access is the current cycle.
  lastL2Access = (unsigned short)(get_current_cycle(0) & 0xFFFF);

  // Set the starting average L2 to LLC/memory latency to 100 cycles.
  averageLatencyCycles = 100;

  // Set the best prefetcher indices to invalid.
  primary = -1;

  // Loop through the score registers clearing them.
  for(int pf = 0;  pf < PREFETCHER_COUNT; pf++)
    sandboxScoreTotal[pf] = lateScoreTotal[pf] = usefulScoreTotal[pf] = 0;
}

//
// This function is called once for each Mid Level Cache read, and is the entry point for participants' prefetching algorithms.
// 		cpu_num - always pass 0 for competition.
// 		addr - the byte address of the current cache read.
// 		ip - the instruction pointer (program counter) of the instruction that caused the current cache read.
// 		cache_hit - 1 for an L2 cache hit, 0 for an L2 cache miss.
//
// This is only called for L1 cache misses, L2 cache reads (including those caused by an L1 miss). In other words
// the prefetcher lives at the L2 level of cache.
//
void l2_prefetcher_operate(int cpu_num, 
    unsigned long long int addr, 
    unsigned long long int ip, 
    int cache_hit)
{
  std::list<unsigned short>::iterator sandboxIT,
    cyclesToArrivalIT;
  unsigned short elapsedCycles;
  bool scoresReady = false;

  //
  // Calculate how many cycles since the previous L2 access and then shift down the new
  // last L2 access.
  //
  elapsedCycles = (unsigned short)(get_current_cycle(0) & 0xFFFF) - lastL2Access;
  lastL2Access = (unsigned short)(get_current_cycle(0) & 0xFFFF);

  //
  // Check if the L2 access was a miss and update the user MSHR.
  //
  if(cache_hit == 0)
  {
    // Ensure it's not alrady in the MSHR.
    std::list<unsigned short>::iterator addressIT;
    for(addressIT = mshrAddresses.begin(); addressIT != mshrAddresses.end(); ++addressIT)
      // Break if we find the address.
      if(*addressIT == ((addr >> 6) & 0xFFFF))
        break;

    // If the iterator isn't pointing to the end it means it contains the address. If it is
    // then add the address and cycle to the lists.
    if(addressIT == mshrAddresses.end())
    {
      // Insert MISS address into our buffer. If there are too many entries then
      // remove the oldest.
      mshrAddresses.push_back((addr >> 6) & 0xFFFF);
      if(mshrAddresses.size() > L2_MSHR_COUNT)
        mshrAddresses.pop_front();

      // Insert the MISS address cycle (now) into our buffer. If there are too many
      // entries then remove the oldest.
      mshrCycles.push_back(get_current_cycle(0) & 0xFFFF);
      if(mshrCycles.size() > L2_MSHR_COUNT)
        mshrCycles.pop_front();
    }
  }

  //
  // Loop through the prefetchers updating which sandbox predictions have now
  // been filled and making new ones.
  //
  for(int pf = 0; pf < PREFETCHER_COUNT; pf++)
  {
    //
    // Loop through the predictions in the sandbox updating which predictions have "arrived" in
    // the L2 cache (i.e. cyclesToArrive reached zero).
    //
    cyclesToArrivalIT = cyclesToArrival[pf].begin();
    for(sandboxIT = sandbox[pf].begin(); sandboxIT != sandbox[pf].end(); ++sandboxIT)
    {
      // The arrival time value is the number of cycles left until the sandbox prediction would have been filled. This is based on
      // the average latency for L2 demand MISSes.
      if(*cyclesToArrivalIT > elapsedCycles)
        // If it still hasn't arrived then decrement the remaining cycles.
        *cyclesToArrivalIT -= elapsedCycles;
      else
        // The elapsed cycles is equal to or more than the reaminging cycles to arrival so
        // the prefetch has "arrived" in the L2 cache.
        *cyclesToArrivalIT = 0;

      // Increment our cycles to arrival list iterator.
      ++cyclesToArrivalIT;
    }

    //
    // Try and find the current cache access in the sandbox.
    //
    cyclesToArrivalIT = cyclesToArrival[pf].begin();
    for(sandboxIT = sandbox[pf].begin(); sandboxIT != sandbox[pf].end(); ++sandboxIT)
    {
      // Check if the L2 access address is in the sandbox for this prefetcher. If it is break
      // leaving the iterator pointed to it.
      if(((addr >> 6) & 0xFFFF) == *sandboxIT)
        break;

      // Increment our cycle list iterator.
      ++cyclesToArrivalIT;
    }

    //
    // Check if the iterator is pointing to past the end meaning that we looped through the entire
    // list and didn't find the entry.
    //
    if(sandboxIT != sandbox[pf].end())
    {
      // The score for the pure sandbox is the current latency in cycles. As in how many cycles are being
      // saved since this value was prefetched.
      sandboxScoreTotal[pf] += 1;

      // The tooLateScore is the remainging number of cycles until the prefetch would have been filled.
      if(*cyclesToArrivalIT > 0)
        lateScoreTotal[pf] += 1;
    }

    //
    // Check if the sandbox is full.
    //
    if(sandbox[pf].size() == L2_ACCESS_COUNT)
    {
      // Calculate final score.
      usefulScoreTotal[pf] = sandboxScoreTotal[pf] - lateScoreTotal[pf];

      // Clear sandbox and cycles to arrival lists.
      sandbox[pf].clear();
      cyclesToArrival[pf].clear();

      // Mark the evaluation period over.
      scoresReady = true;
    }

    //
    // Have the prefetcher make a sandbox prediction and enter the average latency in cycles
    // into the arrival time buffer. If we have over L2_ACCESS_COUNT entries in either then remove
    // the oldest entries.
    //		
    predictions[pf] = prefetchers[pf]->predict(ip, addr);

    // Insert prediction into the sandbox.
    sandbox[pf].push_back((predictions[pf] >> 6) & 0xFFFF);
    if(sandbox[pf].size() > L2_ACCESS_COUNT)
      sandbox[pf].pop_front();

    // If the value was over our limit (1024) then max it out.
    if(averageLatencyCycles > 0x3FF)
      cyclesToArrival[pf].push_back(0x3FF);
    else
      cyclesToArrival[pf].push_back(averageLatencyCycles & 0x3FF);
    if(cyclesToArrival[pf].size() > L2_ACCESS_COUNT)
      cyclesToArrival[pf].pop_front();
  }

  //
  // If we have new scores calculated then pick the primary prefetcher.
  //
  if(scoresReady == true)
  {
    // Find the primary (optimal distance) prefetcher.
    primary = -1;
    for(int pf = 0; pf < PREFETCHER_COUNT; pf++)
    {
      // Save potential new best prefetcher.
      if(primary == -1 || usefulScoreTotal[pf] > usefulScoreTotal[primary])
        primary = pf;
    }

    // Thresholds to decide if the scores are good enough to prefetch.
    if(sandboxScoreTotal[primary] < ((L2_ACCESS_COUNT >> 2) * 1))
      primary = -1;

    // Clear scores.
    for(int pf = 0; pf < PREFETCHER_COUNT; pf++)
      sandboxScoreTotal[pf] = lateScoreTotal[pf] = usefulScoreTotal[pf] = 0;
  }

  //
  // Ensure the primary prefetcher is valid.
  //
  if(primary != -1)
  {
    // Filter the prediction. If it returns true we can proceed.
    if(filter(addr, predictions[primary]) == true)
    {
      // Check MSHR occupancy level (greater than half occupancy).
      if(get_l2_mshr_occupancy(0) > (L2_MSHR_COUNT >> 1)) 
      {
        // Prefetch into the LLC (L3) because MSHRs are scarce.
        if(l2_prefetch_line(0, addr, predictions[primary], FILL_LLC) == 1)
        {
          // Add to prefetch buffer.
          prefetchBuffer.push_back((predictions[primary] >> 6) & 0xFFFF);
          if(prefetchBuffer.size() > PREFETCH_BUFFER_SIZE)
            prefetchBuffer.pop_front();
        }
      }
      else 
      {
        // MSHRs not too busy, so prefetch into L2.
        if(l2_prefetch_line(0, addr, predictions[primary], FILL_L2) == 1)
        {
          // If we were successful at issuing the prefetch then add the address to our
          // global prefetch buffer.
          prefetchBuffer.push_back((predictions[primary] >> 6) & 0xFFFF);
          if(prefetchBuffer.size() > PREFETCH_BUFFER_SIZE)
            prefetchBuffer.pop_front();
        }
      }
    }
  }
}


void l2_cache_fill(int cpu_num, unsigned long long int addr, int set, int way, int prefetch, unsigned long long int evicted_addr)
{	
  // Iterators.
  std::list<unsigned short>::iterator addressIT,
    cycleIT,
    averageIT;

  //
  // Loop through MSHR address buffer looking for the filled address.
  //
  cycleIT = mshrCycles.begin();
  for(addressIT = mshrAddresses.begin(); addressIT != mshrAddresses.end(); ++addressIT)
  {
    // Check if the current iterator address is equal to the filled address.
    if(*addressIT == ((addr >> 6) & 0xFFFF))
    {	
      // Calculate the number of cycles between when the demand MISS occured and when it was
      // filled (now). This is the instantaneous latency so add it to the buffer.
      mshrLatencies.push_back(((get_current_cycle(0) & 0xFFFF) - *cycleIT) & 0x3FF);
      if(mshrLatencies.size() > LATENCY_AVERAGES)
        mshrLatencies.pop_front();

      // Since the address and it's corresponding cycle have been filled we remove
      // them from our lists.
      mshrAddresses.erase(addressIT);
      mshrCycles.erase(cycleIT);

      // Check if the latency buffer is full and an average can be calculated.
      if(mshrLatencies.size() == LATENCY_AVERAGES)
      {
        // Clear the old average latency.
        averageLatencyCycles = 0;

        // Loop through the latency buffer summing the latencies.
        for(averageIT = mshrLatencies.begin(); averageIT != mshrLatencies.end(); ++averageIT)
        {
          // Add the cycle count.
          averageLatencyCycles += *averageIT;
        }

        // Calculate new average latency.
        averageLatencyCycles = averageLatencyCycles >> LATENCY_AVERAGES_SHIFT;
      }

      // We found the address so break out of the loop.
      break;
    }
    ++cycleIT;
  }
}


void l2_prefetcher_heartbeat_stats(int cpu_num)
{
}


void l2_prefetcher_warmup_stats(int cpu_num)
{
}


void l2_prefetcher_final_stats(int cpu_num)
{
}



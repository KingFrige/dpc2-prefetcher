// Implementation of Temp0, for DPC-2-
//  - Mark Sutherland, Ajaykumar Kannan
//    University of Toronto, 2015
//
//    COMPILATION COMMAND:
//      g++ -Wall -o dpc2sim Tempo_Paper9_Sutherland.cc lib/dpc2sim.a
#include <stdio.h>
#include <assert.h>

// cpp libraries
#include <bitset>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <utility>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <inttypes.h>

extern "C" {
#include "../inc/prefetcher.h"
}

// types and helper functions. change nop_fprintf to fprintf to have 
// all debug prints with the aforementioned header turned ON
void nop_fprintf(FILE* fp,...) { }
#define debug nop_fprintf
#define trace nop_fprintf
#define MSHRfprintf nop_fprintf
#define EVICTdebug nop_fprintf
#define SIGNED_COUNTER long long

typedef unsigned long long int PC_t;
typedef unsigned long long int CacheAddr_t;
#define MTHB_DEPTH 32
#define MTHB_WIDTH 3
#define GHB_WIDTH 7

//#define CHARACTERIZE // this will print MASSIVE data at the end of a trace.
//#define GATHER_ACCESS_PATTERNS
//#define GATHER_TEMPOS
//#define NLP

//#define THROTTLE_ACCURACY
#define THROTTLE_AMAT   // throttle the pfetcher basd on 500 request AMAT
#define RECALIBRATE_DEMANDS 500


// use for tracking the time between misses
SIGNED_COUNTER monoMissCtr;

// accumulator for past 500 demand accesses (limit of exponential moving average of AMAT)
double last_amat = -10.0;
double accumulator = 0;
double alpha = 0.5;
bool accum_init = false;
SIGNED_COUNTER num_demands = 0;

//accurate prefetches in this phase
SIGNED_COUNTER accuratePhase;
//polluted prefetches in this phase
SIGNED_COUNTER pollutionPhase;
// total prefetches in this phase
SIGNED_COUNTER totalPhase;

using std::map;
using std::vector;
using std::pair;
using std::size_t;
using std::cout;
using std::cerr;
using std::hex;
using std::dec;
using std::endl; 

const int theBlockSize(CACHE_LINE_SIZE /* defined in prefetcher.h */); 
const int blocksPerPage  = (PAGE_SIZE / CACHE_LINE_SIZE);
const bool NoRotation = false;
typedef unsigned long long pattern_t;
typedef struct block {
  CacheAddr_t addr;
  int offset;
  int hit;
  unsigned long long cycle;
} block_t;

class CacheAddrInfo_t // classed because I wanted to include a prefetch bit
{
  public:
    PC_t ip;
    CacheAddr_t LineAddr;
    bool PrefetchBit;
    bool referenced;
    int miss_counter;

    CacheAddrInfo_t(CacheAddr_t line,PC_t point,bool pfb,int mc) : 
      ip(point),
      LineAddr(line),
      PrefetchBit(pfb),
      referenced(false),
      miss_counter(mc)
  { }

    // copy constructor
    CacheAddrInfo_t(const CacheAddrInfo_t &rhs) :
      ip(rhs.ip),
      LineAddr(rhs.LineAddr),
      PrefetchBit(rhs.PrefetchBit),
      referenced(rhs.referenced),
      miss_counter(rhs.miss_counter)
  { }

    // equality 
    bool operator==(const CacheAddrInfo_t &other) const {
      return (LineAddr == other.LineAddr);
    }
};

// comparison
bool operator<(const CacheAddrInfo_t &i1, const CacheAddrInfo_t &i2) {
  return (i1.LineAddr < i2.LineAddr);
}

class PrefetchData_t {
  public:
    PC_t LastRequestAddr;
    CacheAddr_t DataAddr;
    bool hit;
    int miss_counter;

    PrefetchData_t (PC_t ip, CacheAddr_t mia, bool didwe,int mc) :
      LastRequestAddr(ip),
      DataAddr(mia),
      hit(didwe),
      miss_counter(mc)
  { }
};

struct Stats
{
  typedef std::map<std::string,long long> tCounters;
  tCounters theCounters;
  ~Stats() {
    for(tCounters::iterator i=theCounters.begin();i!=theCounters.end();++i) {
      std::cout << i->first << "," << i->second << std::endl;
    }
  }
} theStats;
#define C(n) do { ++(theStats.theCounters[(#n)]); } while(0);
#define S(n,v) do { (theStats.theCounters[(#n)])+=(v); } while(0);

// templated hash table
template<class KeyType,class ItemType>
struct Container
{
  int theHeight,theWidth,theKeyShift;
  unsigned long long theTagMask; // dictates number of bits used to tag elements (includes lower index bits for simulation purposes only)
  unsigned long long theIndexMask; // based only on Height of structure

  typedef std::pair<KeyType,ItemType> Item;
  typedef std::list<Item> ListType;
  std::vector<ListType> theItems;
  typedef typename ListType::iterator Iter;
  Container(
      int aHeight
      ,	int aWidth
      ,	int aKeyShift
      ,	int aTagBits
      )
    :	theHeight(aHeight)
      ,	theWidth(aWidth)
      ,	theKeyShift(aKeyShift)
      ,	theTagMask((1ULL<<aTagBits)-1)
      ,	theIndexMask(theHeight-1)
  {
    assert(!(theHeight & theIndexMask));
    theItems.resize(theHeight);
    for(int i=0;i<theHeight;++i)
      theItems[i].resize(theWidth);
  }

  ~Container() {
    /*for(int i=0;i<theHeight;++i)
      for(int j=0;j<theWidth;++i)
      theItems[i].push_back(new ItemType());*/
  }

  Item insert(KeyType aKey,ItemType anItem) {
    int aKeyIndex(index(aKey));
    Item anOldItem(*theItems[aKeyIndex].rbegin());
    theItems[aKeyIndex].pop_back();
    KeyType aTag = tag(aKey);
    theItems[aKeyIndex].push_front(std::make_pair(aTag,anItem));
    theItems[aKeyIndex].front().second.EVICT_DETECTOR = aKey;
    //debug(stderr,"\tinsert %llx [%d:%d/%d]\n",aKey,index(aKey),theHeight,theWidth);
    return anOldItem;
  }

  void checkSize(KeyType aKey) {
    int aKeyIndex(index(aKey));
    ListType thisList = theItems[aKeyIndex];
    std::cerr << "LIST SIZE: " << thisList.size() << std::endl;
  }


  Iter end() {
    return theItems[0].end();
  }

  static inline unsigned long long inthash(unsigned long long key)
  {
    key += (key << 12);
    key ^= (key >> 22);
    key += (key << 4);
    key ^= (key >> 9);
    key += (key << 10);
    key ^= (key >> 2);
    key += (key << 7);
    key ^= (key >> 12);
    return key;
  }

  int index(KeyType aKey) {
    return (inthash(aKey >> theKeyShift) & theIndexMask);
  }

  int tag(KeyType aKey) {
    return (inthash(aKey >> theKeyShift) & theTagMask);
  }

  Iter find(KeyType aKey) {
    ListType& aList(theItems[index(aKey)]);
    KeyType aTag = tag(aKey);
    for(Iter i=aList.begin();i!=aList.end();i++) {
      //debug(stderr,"\tsearch for %llx in %llx [%d:%d/%d]\n",aKey,i->first,index(aKey),theHeight,theWidth);
      if (i->first == aTag) {
        aList.push_front(*i);
        aList.erase(i);
        return aList.begin();
      }
    }
    return end();
  }

  bool erase(KeyType aKey) {
    ListType& aList(theItems[index(aKey)]);
    KeyType aTag = tag(aKey);
    for(Iter i=aList.begin();i!=aList.end();i++) {
      if (i->first == aTag) {
        aList.erase(i);
        aList.resize(theWidth);
        return true;
      }
    }
    return false;
  }
};

struct MSHRs {
  unsigned int theSize;
  unsigned long long int PG_MASK;
  unsigned long long int BLK_MASK;
  typedef std::map<CacheAddrInfo_t,SIGNED_COUNTER> tEntries;
  tEntries theEntries;

  MSHRs(unsigned int aSize)
    :	theSize(aSize)
  { 
    PG_MASK = ~(CACHE_LINE_SIZE - 1);
    BLK_MASK = (CACHE_LINE_SIZE -1);
  }

  bool resetPFBits(PrefetchData_t& req,SIGNED_COUNTER newCycle) {
    // if we already prefetched this address and then saw a demand hit, reset
    // the PF bit and return true.
    CacheAddrInfo_t tmp(req.DataAddr,0,0,0);
    tEntries::iterator i = theEntries.find(tmp);
    if(i == theEntries.end()) return false; // did not reset anything
    else { 
      if( (i->first).PrefetchBit == false) return false; // if mshr but no pf, return.
      else {
        CacheAddrInfo_t reInsertMe(i->first);
        theEntries.erase(i); // remove old
        reInsertMe.PrefetchBit = false;
        theEntries.insert(std::make_pair(reInsertMe,newCycle)); // add new with cleared pf bit
        return true; // hit an MSHR that was already allocated as a prefetch
      }
    }
  }

  void clean(SIGNED_COUNTER cycle) {
    if (cycle<0) cycle = -cycle;
    tEntries::iterator next;
    for(tEntries::iterator i=theEntries.begin();i!=theEntries.end();i = next) {
      next = i;
      ++next;
      if (GetPrefetchBit(i->first)) {
        MSHRfprintf(stderr,"removing MSHR entry %llx, %ld left (i->second=%lld, cycle=%lld)\n",(i->first).LineAddr,theEntries.size()-1,i->second,cycle);
        if (i->second > 0) {
          S(mem_waited,cycle - i->second)
            assert(cycle >= i->second);
        }
        theEntries.erase(i);
      } else { // cannot remove demand requests
        MSHRfprintf(stderr,"keeping MSHR entry %llx, %ld entries\n",(i->first).LineAddr,theEntries.size());
      }
    }
  }

  bool inflight(CacheAddr_t anAddr) {
    anAddr &= (PG_MASK);
    CacheAddrInfo_t tmp(anAddr,0,0,0);
    bool res = (theEntries.find(tmp) != theEntries.end());
    MSHRfprintf(stderr,"MSHR in-flight check for %llx: %d\n",anAddr,res);
    return res;
  }

  bool allocate(PC_t ip,CacheAddr_t address, SIGNED_COUNTER cycle, bool prefetch,int miss_counter) {
    address &= PG_MASK;
    CacheAddrInfo_t anAddr(address,ip,prefetch,miss_counter);
    if(prefetch) {
      clean(cycle);
    }
    //if (cycle<0 && GetPrefetchBit(anAddr)) C(L2_Prefetch_Not_Needed)

    tEntries::iterator i = theEntries.find(anAddr);
    if (i != theEntries.end()) {
      if (cycle<0) {
        CacheAddrInfo_t allocd = (*i).first;
        if(GetPrefetchBit(allocd)) {
          C(L2_Prefetch_Already_InFlight);
        } else {
          C(L2_Request_Already_InFlight);
        }
      } else if (i->second<0) {
        C(L2_Prefetch_PartialHit);
        S(mem_saved,cycle + i->second)
          i->second = cycle;
      } else {
        C(L2_OoO_Miss_Overlap);
      }
      MSHRfprintf(stderr,"MSHR for %llx exists\n",address);
      return false;
    }
    if (theEntries.size() == theSize) {
      MSHRfprintf(stderr,"MSHRs are full, can't insert %llx (size=%d)\n",address,theSize);
      C(MSHR_full)
        return false;
    }
    theEntries.insert(std::make_pair(anAddr,cycle));
    MSHRfprintf(stderr,"added MSHR entry %llx, %ld present\n",address,theEntries.size());
    return true;
  }

  bool GetPrefetchBit(CacheAddrInfo_t anAddr) { 
    return anAddr.PrefetchBit;
  }

  bool GetExternalPrefetchBit(CacheAddr_t in) {
    CacheAddrInfo_t tmp(in,0,0,0);
    tEntries::iterator i = theEntries.find(tmp);
    if(i == theEntries.end()) {
      MSHRfprintf(stderr,"Getting prefetch bit for: %llx. NOT FOUND.\n",in);
      return false;
    }
    MSHRfprintf(stderr,"Getting prefetch bit for: %llx. bit=%d.\n",in,GetPrefetchBit(i->first));
    return GetPrefetchBit(i->first);
  }

  std::pair<CacheAddrInfo_t,SIGNED_COUNTER> deallocMSHR(CacheAddr_t address, int prefetch,SIGNED_COUNTER cycle)
  {
    address &= PG_MASK;
    CacheAddrInfo_t anAddr(address,0,prefetch,0);
    tEntries::iterator i = theEntries.find(anAddr);
    if(i == theEntries.end()) {
      MSHRfprintf(stderr,"cycle=%llu ERROR: Tried to dealloc MSHR for address: %llx, matching entry not found?!\n",cycle,address);
      return (std::make_pair(anAddr,-1));
    }
    else{
      assert(cycle >= i->second);
      CacheAddrInfo_t returnMe(i->first); // copy construct this
      theEntries.erase(i);
      MSHRfprintf(stderr,"Line fill removed mshr entry: %llx, size=%ld\n",address,theEntries.size());
      return (std::make_pair(returnMe,i->second));
    }
  }

  int getNumUniquePCs(void) {
    tEntries::iterator i;
    int ret = 0;
    std::vector<PC_t> seen;

    for(i = theEntries.begin(); i != theEntries.end();i++) {
      CacheAddrInfo_t me = (*i).first;
      std::vector<PC_t>::iterator j = std::find(seen.begin(),seen.end(),me.ip);
      if(j == seen.end()) { // not found
        seen.push_back(me.ip);
        ret++;
      }
    }
    return ret;
  }

  int numOutstanding() {
    return theEntries.size();
  }
};

#ifdef GATHER_TEMPOS
class TempoBuckets {
  public:
    std::vector<int> local_vec;
    std::vector<int> global_vec;

    TempoBuckets() { 
      local_vec.assign(4,0);
      global_vec.assign(4,0);
    }
    TempoBuckets(vector<int>l,vector<int>g):local_vec(l),global_vec(g) { }
    TempoBuckets(const TempoBuckets& rhs):local_vec(rhs.local_vec),global_vec(rhs.global_vec) { }
};

class TempoPair {
  public:
    int lhb_tempo;
    int ghb_tempo;
    TempoPair() {} 
    TempoPair(int loc,int glob) : lhb_tempo(loc),ghb_tempo(glob) { }
    TempoPair(const TempoPair& rhs) : lhb_tempo(rhs.lhb_tempo),ghb_tempo(rhs.ghb_tempo) { }
};

typedef map<SIGNED_COUNTER,TempoPair> CycleToTempoMap;
typedef map<PC_t,TempoBuckets> PCToTempoMap;
PCToTempoMap tempo_map;
#endif

#ifdef GATHER_ACCESS_PATTERNS
std::map<unsigned long long, std::vector<block_t> > map_pc;
std::map<unsigned long long, std::vector<block_t> > map_page;
#endif

#ifdef CHARACTERIZE
PC_t pcToCharacterize = 0x000000000040A6F1;
class TrainingInfoClass
{
  public:
    int tempo_way;
    pattern_t pattern_trained;
    bool to_pht;

    TrainingInfoClass() :
      tempo_way(0),
      pattern_trained(0),
      to_pht(false)
  { }

    TrainingInfoClass(const TrainingInfoClass& rhs) :
      tempo_way(rhs.tempo_way),
      pattern_trained(rhs.pattern_trained),
      to_pht(rhs.to_pht)
  { }
};

class PrefetchedInfoClass 
{
  public:
    int tempo_way;
    std::vector<CacheAddr_t> pfetches;

    PrefetchedInfoClass() :
      tempo_way(0)
  { }

    PrefetchedInfoClass(const PrefetchedInfoClass& rhs) :
      tempo_way(rhs.tempo_way),
      pfetches(rhs.pfetches)
  { }
};

typedef pair<TrainingInfoClass,PrefetchedInfoClass> CycleInfoPair;
typedef map<SIGNED_COUNTER,CycleInfoPair> CycleToInfoMap;
typedef map<PC_t,CycleToInfoMap> PCToCycleInfoMap;

PCToCycleInfoMap big_pc_map;
#endif

MSHRs *mshrs;
//MSHRs *prefetchMSHRs; // these are not real mshrs, only used to track pf requests

class NextLinePrefetcher
{
  public:
    CacheAddr_t get_pf_addr(CacheAddr_t addr,int n) {
      return (((addr>>6)+n)<<6);
    }
};

class PrefetchBitsArray
{
  public:
    int pfbits[L2_SET_COUNT][L2_ASSOCIATIVITY];

    PrefetchBitsArray() {
      for(int i = 0; i < L2_SET_COUNT; i++)
        for(int j = 0; j < L2_ASSOCIATIVITY;j++)
          pfbits[i][j] = 0;
    }
};

template<int NENTRIES,int BUFFER_WIDTH>
class MissTimeHistoryBuffers
{
  private:
    uint32_t ghbs[NENTRIES][BUFFER_WIDTH];
    uint32_t last_times[NENTRIES];
    uint32_t hashMask;

  public:
    void init() {
      for(int i = 0;i<NENTRIES;i++) {
        for(int j = 0;j<BUFFER_WIDTH;j++) {
          ghbs[i][j] = 0;
        }
      }

      // setup hash mask
      hashMask = (NENTRIES-1);
    }

    int hashCode(PC_t key) {
      key += (key << 12);
      key ^= (key >> 22);
      key += (key << 4);
      key ^= (key >> 9);
      key += (key << 10);
      key ^= (key >> 2);
      key += (key << 7);
      key ^= (key >> 12);
      return (key & hashMask);
    }

    int updateLastTime(int row, uint32_t last) {
      int thelasttime = last_times[row];
      last_times[row] = last;
      return (last - thelasttime);
    }

    void update(PC_t key, uint32_t accessTimeDelta) {
      int row = hashCode(key);

      // ror [ghb],[ghb],1
      for(int i = BUFFER_WIDTH-1;i>0;i--) {
        //debug(stderr,"accessing [row][%d] and [row][%d]\n",i,i-1);
        ghbs[row][i] = ghbs[row][i-1];
      }
      ghbs[row][0] = accessTimeDelta;
    }

    uint32_t translateNumToBin(uint32_t in) {
      // translates a latency (in terms of misses expired) to a "bin" representing long/short
      return 0;
      /*
         if(in <=2 ) return 0;
         if(in <= 9) return 1;
         if(in <= 19) return 2;
         return 3;
         */
    }

    uint32_t getFunctDelta(PC_t key) {
      int row = hashCode(key);

      // return a function of the GHB
      uint32_t ret = 0;
      for(int i =0;i<BUFFER_WIDTH-1;i++) {
        ret += (ghbs[row][i] - ghbs[row][i+1]);
      }
      ret /= BUFFER_WIDTH; // INTERGER DIVISION
      return ret;
    }

    uint32_t getFunct(PC_t key) {
      int row = hashCode(key);

      // return a function of the GHB
      uint32_t ret = 0;
      for(int i =0;i<BUFFER_WIDTH;i++) {
        ret += ghbs[row][i];
      }
      ret /= BUFFER_WIDTH; // INTERGER DIVISION
      return ret;
    }

    void printHistBuffers() {
      for(int i = 0;i<NENTRIES;i++) {
        printf("BUF=%d: ",i);
        for(int j = 0;j<BUFFER_WIDTH;j++) {
          printf("%d ",ghbs[i][j]);
        }
        printf("\n");
      }
    }
};

class SMS
{
  public:
    unsigned long long int theRegionShift,theRegionSize,theRegionMask,theBlocksPerRegion;

    std::set<CacheAddr_t> pollutionAddrs;

    // set once in l2_prefetcher_operate
    int max_l2_prefetches;
    int max_llc_prefetches;
    int l2_issued;
    int llc_issued;
    int degree_llc,degree_l2;

    struct AGTent {
      CacheAddr_t pc;
      int offset;
      pattern_t pattern;
      AGTent(CacheAddr_t aPC = 0ULL,int aOffset = 0,pattern_t aPattern = 0ULL)
        :	pc(aPC)
          ,	offset(aOffset)
          ,	pattern(aPattern)
      { }
      CacheAddr_t EVICT_DETECTOR;

      AGTent(const AGTent& rhs) :
        pc(rhs.pc),
        offset(rhs.offset),
        pattern(rhs.pattern),
        EVICT_DETECTOR(rhs.EVICT_DETECTOR)
      { }

      void PrintMe(void) {
        debug(stderr,"Printing AGT entry address: %p\n",this);
        debug(stderr,"\tPC=%llx, offset=%d, pattern=%llu,EVICT_DETECTOR=%llx\n",pc,offset,pattern,EVICT_DETECTOR);
      }
    };
    struct PHTent {
      pattern_t pattern;
      PHTent(pattern_t aPattern = 0ULL):pattern(aPattern) { }
      CacheAddr_t EVICT_DETECTOR;

      PHTent(const PHTent& rhs) :
        pattern(rhs.pattern),
        EVICT_DETECTOR(rhs.EVICT_DETECTOR)
      { }

    };
    struct EntryData {
      PHTent pht_entry;
      SIGNED_COUNTER add_cycle;
      SIGNED_COUNTER evict_cycle;

      EntryData(const EntryData& rhs) :
        pht_entry(rhs.pht_entry),
        add_cycle(rhs.add_cycle),
        evict_cycle(rhs.evict_cycle)
      { }

      EntryData(const PHTent& copy_entry, 
          const SIGNED_COUNTER copy_cycle) :
        pht_entry(copy_entry),
        add_cycle(copy_cycle),
        evict_cycle(0)
      { }
    };
    typedef Container<CacheAddr_t,AGTent> AGT; 
    AGT* theAGT[4];
    typedef Container<CacheAddr_t,PHTent> PHT; 
    PHT* thePHT[4];

    typedef std::list<EntryData> thePHTEntries;
    typedef std::map<CacheAddr_t,thePHTEntries> characterization;
    characterization theCharacterization;

    typedef std::map<SIGNED_COUNTER,SIGNED_COUNTER> pc_counter_type;
    pc_counter_type outstandingIPs;

    typedef std::pair<int,int> tempoWayBitPairs;
    typedef std::list<tempoWayBitPairs> tempoLists;
    typedef std::map<PC_t,tempoLists> tempoMapType;
    tempoMapType tempoMap;


    SMS()
      : theRegionShift(log2(PAGE_SIZE)) // shift: 12 = 4KB region, 10 = 1KB region
        , theRegionSize(PAGE_SIZE)
        , theRegionMask(theRegionSize-1)
        , theBlocksPerRegion(theRegionSize/theBlockSize)
  { 
    for(int i = 0;i<4;i++) { // init AGT and PHT
      theAGT[i] = new AGT(64,2,theRegionShift,64-theRegionShift);
      thePHT[i] = new PHT(256,2,NoRotation?0:2,14);
    }
    //degree_llc = 4;
    degree_l2 = 2;
  }

    ~SMS() {
      // free all AGT and PHT
      for(int i = 0;i<4;i++) {
        delete theAGT[i];
        delete thePHT[i];
      }
    }

    void printCharacterization() {
      // foreach key, go over the list and print all PHT entries
      characterization::iterator i;
      for(i = theCharacterization.begin(); i != theCharacterization.end(); i++) {
        CacheAddr_t thisPC = i->first;
        thePHTEntries listOfPHTents = i->second;
        thePHTEntries::iterator j;
        // print header for this pc
        cout << "\tPC: " << hex << thisPC << "\tNum PHT's: " << dec << listOfPHTents.size();
        for(j = listOfPHTents.begin(); j != listOfPHTents.end(); j++) {
          EntryData cur = *j;  
          std::bitset<blocksPerPage> patternbits(cur.pht_entry.pattern);
          cout << endl << patternbits << "\tAdd_cyc=" << cur.add_cycle << " Evict_cyc=" << cur.evict_cycle;
        }
        cout << endl;
      } 

      cout << "--------------------------------------------------------------------" << endl;

      pc_counter_type::iterator j;
      for(j = outstandingIPs.begin();j != outstandingIPs.end();j++) {
        SIGNED_COUNTER cyc = (*j).first;
        SIGNED_COUNTER pcs = (*j).second;

        cout << "Cycle " << cyc << ": " << pcs << endl;
      }

      cout << "*********************************************************************" << endl;
      // print time evolution of all tempo entries
      tempoMapType::iterator mapIt;
      for(mapIt = tempoMap.begin(); mapIt != tempoMap.end(); mapIt++){
        cout << "PC: 0x" << (*mapIt).first << " ############## " << endl;
        tempoLists curList = (*mapIt).second;
        tempoLists::iterator listIt;
        for(listIt = curList.begin();listIt != curList.end(); listIt++) {
          tempoWayBitPairs curPair = *listIt;
          cout << " <Way: " << curPair.second << ",Bits: " << curPair.first << ">"<<endl;
        }
        cout << endl;
      }
    }

    void insertToCharacterization(const CacheAddr_t& key,const EntryData& insert_me) {
      //cout << "i'm in characterization insert" << endl;
      size_t index = PHT::inthash(key);
      thePHTEntries retList = theCharacterization[index];
      retList.push_back(insert_me);
      theCharacterization[index] = retList;
    }

    void updatePHTEviction(const CacheAddr_t& key, const SIGNED_COUNTER eCyc) {
      size_t index = PHT::inthash(key);
      thePHTEntries retList = theCharacterization[index]; // this must be set already
      if(retList.size() == 0)
        return;
      // get iterator to END (last entry) and set eviction time
      thePHTEntries::reverse_iterator i = retList.rbegin();
      EntryData lastEntry = *i;
      lastEntry.evict_cycle = eCyc;
      retList.pop_back();
      retList.push_back(lastEntry);
      theCharacterization[index] = retList;
    }

    // once you are IN the characterization, YOU CAN'T LEAVE!
    // void deleteFromCharacterization(const CacheAddr_t) {
    // }

    pattern_t rotate(int aBitIndex,int anOffset) {
      debug(stderr,"rotate(%d,%d [%llu]) : ",aBitIndex,anOffset,theBlocksPerRegion);
      pattern_t res = 1ULL<<((aBitIndex + anOffset) % theBlocksPerRegion);
      debug(stderr,"%llx > %llx\n",1ULL<<aBitIndex,res);
      assert(!(res & (res-1)));
      return res;
    }

    bool replace( SIGNED_COUNTER cycle, CacheAddr_t addr, int agtWay) {
      CacheAddr_t region_tag = addr & ~theRegionMask;
      int region_offset = (addr & theRegionMask)>>6;
      AGT::Item agt_evicted;
      //DBG_(Dev, ( << std::dec << theId << "-evict: group=" << std::hex << region_tag << "  offset=" << std::dec << region_offset ) );
      bool erased_something(false);

      AGT::Iter agt_ent = theAGT[agtWay]->find(region_tag);
      if (agt_ent != theAGT[agtWay]->end()) {
        pattern_t new_bit = rotate(region_offset,agt_ent->second.offset);
        if(NoRotation) new_bit = 1ULL << region_offset;
        if (agt_ent->second.pattern & new_bit) {
          agt_evicted = *agt_ent;
          //DBG_(Dev, ( << std::dec << theId << "-end: group=" << std::hex << region_tag << "  key=" << agt_evicted.second.pc << "  " << agt_evicted.second.pattern ) );
          theAGT[agtWay]->erase(region_tag);
          erased_something = true;
        }
      }

      if (agt_evicted.second.pattern) {
        if (thePHT[agtWay]->erase(agt_evicted.second.pc)) {
          C(PHT_erased_previous); // if replacing or if it's a singleton
                                  // look through PHT table for this matching key and update its eviction #
#ifdef CHARACTERIZE
          debug(stderr,"looking for pc=%llx in pht to update eviction number. current region tag (pc) is=%llx\n",agt_evicted.second.pc,region_tag);
          updatePHTEviction(agt_evicted.second.pc,cycle);
#endif
        }
        if ((agt_evicted.second.pattern-1)&agt_evicted.second.pattern) {// not singleton
          debug(stderr,"learned pattern (block eviction) %llx into PHT, pc=%llx\n",agt_evicted.second.pattern,agt_evicted.second.pc);
          thePHT[agtWay]->insert(agt_evicted.second.pc,PHTent(agt_evicted.second.pattern));
          // every time you insert to the PHT, add an entry into this other table
          // that stores the pattern, and the cycle it was added
#ifdef CHARACTERIZE
          EntryData newPHTEntry(PHTent(agt_evicted.second.pattern),cycle);
          insertToCharacterization(agt_evicted.second.pc,newPHTEntry);
#endif
        }
      }

      return erased_something;
    }

    void checkEvictions(SIGNED_COUNTER cycle,MSHRs* mshrs,CacheAddr_t evicted_addr) {
      // we are GIVEN the evicted address.... so simply check if this maps to any entry
      // that is already been entered in the AGT, and if so, then move it to the PHT
      //debug(stderr,"Searching the AGT for any active entries corresponding to evicted_addr = %llx\n",evicted_addr); 
      //
      //TODO: Have to check ALL of the AGT entries (4 way, remember?)
      for(int agtWay=0;agtWay<4;agtWay++) {
        for(int n=0;n<theAGT[agtWay]->theHeight;n++) {
invalidated_list:
          //debug(stderr,"Searching new list @ index=%d\n",n);
          AGT::ListType& aList(theAGT[agtWay]->theItems[n]);
          int x=0;
          for(AGT::Iter i = aList.begin();i != aList.end();i++) {
            int offset = theBlocksPerRegion-1;
            for(pattern_t pattern = i->second.pattern;pattern;--offset) {
              pattern_t mask = (1ULL<<offset);
              if (!(pattern & mask)) continue;
              EVICTdebug(stderr,"  attempt offset=%d mask=%llx on region_offset=%d with pattern %llx\n",offset,mask,i->second.offset,pattern);
              pattern &= ~mask;
              CacheAddr_t prediction = ((-i->second.offset+offset)*theBlockSize);
              if(NoRotation) prediction = (offset*theBlockSize);
              EVICTdebug(stderr,"  prediction = %llx\n",prediction);
              prediction &= theRegionMask;
              CacheAddr_t anAddress = i->second.EVICT_DETECTOR+prediction;
              EVICTdebug(stderr,"  anAddress: %llx\n",anAddress);
              EVICTdebug(stderr,"  prediction&= %llx\n",prediction);
              EVICTdebug(stderr,"  prediction!= %llx\n",i->first + prediction);

              if(anAddress == evicted_addr) {
                MSHRfprintf(stderr,"detected evicted block (n=%d, x=%d), region=%llx offset=%d, evicted addr=%llx\n",n,x,i->first,i->second.offset,anAddress);
                if (replace(cycle, anAddress,agtWay)) goto invalidated_list;
              }
            } // end pattern offset loop
            x++;
          }// end AGT list loop
        } // end loop over AGT rows
      }
    }

    bool updateAGTs(SIGNED_COUNTER cycle, PrefetchData_t* Data, MSHRs* mshrs, int locHash) {
      CacheAddr_t pc(Data->LastRequestAddr);
      CacheAddr_t region_tag = Data->DataAddr & ~theRegionMask;
      int region_offset = (Data->DataAddr & theRegionMask)>>6;
      CacheAddr_t key = pc;
      if(NoRotation) key = (pc << (theRegionShift-6)) | region_offset;
      bool miss = (! Data->hit);

      //debug(stderr,"Called IssuePrefetches with theID: %d-access: group= %llx, key= %llx, offset= %llu\n",theId,region_tag,key,region_offset);
      //DBG_(Dev, ( << std::dec << theId << "-access: group=" << std::hex << region_tag << "  key=" << key << "  offset=" << std::dec << region_offset ) );
      debug(stderr,"In updateAGTs with: region=%llx offset=%d bit=%llx pc=%llx lochash=%d\n",region_tag,region_offset,1ULL<<region_offset,pc,locHash);
#ifdef CHARACTERIZE
      bool sent_to_pht = false;
      CycleToInfoMap curMap;
      CycleInfoPair curInfoPair;
      if (pc == pcToCharacterize) {
        curMap = big_pc_map[pc];
        curInfoPair = curMap[monoMissCtr];
      }
#endif
      bool new_gen = false;
      //for(int j=0;j<4;j++) { // do this for all agt regions (to keep consistent state)
      debug(stderr,"Checking for evict/replace for way# %d\n",locHash);
      AGT::Iter agt_ent = theAGT[locHash]->find(region_tag);
      AGT::Item agt_evicted;
      if (agt_ent == theAGT[locHash]->end()) {
        C(AGT_evict_replacement);
        pattern_t new_bit = 1ULL<<(theBlocksPerRegion-1);
        if(!NoRotation) {
          agt_evicted = theAGT[locHash]->insert(region_tag,AGTent(key,theBlocksPerRegion-region_offset-1,new_bit));
        } else {
          new_bit = 1ULL << region_offset;
          agt_evicted = theAGT[locHash]->insert(region_tag,AGTent(key,0,new_bit));
        }
        new_gen = true;
        debug(stderr,"new pattern (from scratch) new_bit=%llx offset=%d->%llu\n",new_bit,region_offset,theBlocksPerRegion-region_offset);
      } else {
        pattern_t new_bit = rotate(region_offset,agt_ent->second.offset);
        if(NoRotation) new_bit = 1ULL << region_offset;
        if ((agt_ent->second.pattern & new_bit) && miss) {
          C(AGT_samebit_replacement)
            debug(stderr,"collided on pattern %llx (new_bit=%llx)\n",agt_ent->second.pattern,new_bit);
        } else {
          C(AGT_addbit);
        }
        agt_ent->second.pattern |= new_bit;
        debug(stderr,"update pattern:%llx new_bit=%llx offset=%d\n",agt_ent->second.pattern,new_bit,agt_ent->second.offset);
      }

      if (agt_evicted.second.pattern) {
        if (thePHT[locHash]->erase(agt_evicted.second.pc)){
          C(PHT_erased_previous); // if replacing or if it's a singleton
#ifdef CHARACTERIZE_OLD
          debug(stderr,"looking for pc=%llx in pht to update eviction number. current region tag (pc) is=%llx\n",agt_evicted.second.pc,region_tag);
          updatePHTEviction(agt_evicted.second.pc,cycle);
#endif
        }
        if ((agt_evicted.second.pattern-1)&agt_evicted.second.pattern) {
          // not singleton
          debug(stderr,"learned pattern (AGT eviction) %llx into PHT, pc=%llx\n",agt_evicted.second.pattern,agt_evicted.second.pc);
          thePHT[locHash]->insert(agt_evicted.second.pc,PHTent(agt_evicted.second.pattern));
#ifdef CHARACTERIZE
          sent_to_pht = true;
#endif
#ifdef CHARACTERIZE_OLD
          EntryData newPHTEntry(PHTent(agt_evicted.second.pattern),cycle);
          insertToCharacterization(agt_evicted.second.pc,newPHTEntry);
#endif
        }
      }
      // store data at this cycle
#ifdef CHARACTERIZE
      TrainingInfoClass newTrainingInfo;
      newTrainingInfo.tempo_way = locHash;
      newTrainingInfo.pattern_trained = agt_ent->second.pattern;
      newTrainingInfo.to_pht = sent_to_pht;
      // map it up
      curInfoPair.first = newTrainingInfo;
      curMap[monoMissCtr] = curInfoPair;
      big_pc_map[pc] = curMap;
#endif
      return new_gen;
    }

    void IssuePrefetches( SIGNED_COUNTER cycle,PrefetchData_t *Data,int ghbHash,bool new_gen ) {
      CacheAddr_t pc(Data->LastRequestAddr);
      CacheAddr_t region_tag = Data->DataAddr & ~theRegionMask;
      int region_offset = (Data->DataAddr & theRegionMask)>>6;
      CacheAddr_t key = pc;
      if(NoRotation) key = (pc << (theRegionShift-6)) | region_offset;
      debug(stderr,"In IssuePrefetches with: region=%llx region_offset=%d bit=%llx pc=%llx,ghbHash=%d,missCtr=%d\n",region_tag,region_offset,1ULL<<region_offset,pc,ghbHash,monoMissCtr);
      // now we are done working with the AGT and PHT (for all ways)... only search the PHT
      // that matches this tempo
#ifdef CHARACTERIZE
      CycleToInfoMap curMap;
      CycleInfoPair curInfoPair;
      if (pc == pcToCharacterize) {
        curMap = big_pc_map[pc];
        curInfoPair = curMap[monoMissCtr];
      }
#endif
      if (new_gen) {
        PHT::Iter pht_ent = thePHT[ghbHash]->find(key);
        if (pht_ent != thePHT[ghbHash]->end()) {
          //DBG_(Dev, ( << std::dec << theId << "-predict: group=" << std::hex << region_tag << "  key=" << key << "  " << pht_ent->second.pattern ) );
          C(L1_Found_Pattern);
          debug(stderr,"prediction pattern for pc=%llx is %llx and region_offset %d\n",key,pht_ent->second.pattern,region_offset);
          //TODO int offset = theBlocksPerRegion-2; // extra -1 to avoid prefetch of trigger
          int offset = 1;
          if(NoRotation) offset += 1;

          //assert((pht_ent->second.pattern-1)&pht_ent->second.pattern);
#ifdef CHARACTERIZE_OLD
          // count the average number of bits in a "tempo way"
          std::bitset<blocksPerPage> curTemp0Way(pht_ent->second.pattern);
          int numBits = curTemp0Way.count();
          tempoWayBitPairs newBitPair = std::make_pair(numBits,ghbHash);
          //size_t index = PHT::inthash(pc);
          tempoLists retList = tempoMap[pc]; // this is already set
          retList.push_back(newBitPair);
          tempoMap[pc] = retList;
#endif
#ifdef CHARACTERIZE
          std::vector<CacheAddr_t> blocksPrefetched;
#endif
          //for(pattern_t pattern = pht_ent->second.pattern;pattern;--offset) {
          for(pattern_t pattern = pht_ent->second.pattern;pattern;++offset) { // INCREMENT OVER PG rather than DEC
                                                                              //debug(stderr,"  prediction at offset %d, pattern is %llx and region_offset %d\n",offset,pht_ent->second.pattern,region_offset);
            pattern_t mask = (1ULL<<offset);
            debug(stderr,"  attempt offset=%d mask=%llx on region_offset=%d with pattern %llx\n",offset,mask,region_offset,pattern);
            if (!(pattern & mask)) continue;
            pattern &= ~mask;
            if (NoRotation && (offset == region_offset)) continue;
            debug(stderr,"  prediction at offset %d, pattern is %llx and region_offset %d\n",offset,pht_ent->second.pattern,region_offset); CacheAddr_t prediction = ((region_offset+offset+1)*theBlockSize);
            if(NoRotation) prediction = (offset*theBlockSize);
            debug(stderr,"  prediction = %llx\n",prediction);
            prediction &= theRegionMask;
            debug(stderr,"  prediction&= %llx\n",prediction);
            debug(stderr,"  prediction!= %llx\n",region_tag + prediction);

            bool pfret = false;
            bool inflight;
            if (l2_issued < max_l2_prefetches ) { // TODO
              debug(stderr,"\tTRYING NEW PREFETCH: Base Address=0x%llx, Prefetch Address=0x%llx, into L2\n",Data->DataAddr,(region_tag+prediction));
              // don't prefetch a line where there's already a demand pending
              inflight = mshrs->inflight(region_tag+prediction);
              // don't prefetch a line where there's already a prefetch pending
              //inflight = prefetchMSHRs->inflight(region_tag+prediction);
              if(!inflight) {
                pfret = l2_prefetch_line(0,Data->DataAddr,(region_tag+prediction),FILL_L2);
              } else {
                C(Already_Outstanding);
              }
              if(pfret) {
                //prefetchMSHRs->allocate(pc,region_tag+prediction,cycle,true,Data->miss_counter); these were just to track stats without cleaning MSHRs every cycle
                mshrs->allocate(pc,region_tag + prediction, cycle,true,Data->miss_counter);
                l2_issued++;
                totalPhase++;
                C(L2_Prefetches_Issued); 
#ifdef CHARACTERIZE
                blocksPrefetched.push_back(region_tag+prediction);
#endif
              } else  {
                C(L2_Prefetches_Dropped);
              }
              debug(stderr,"  PREF returned to L2: %d\n",pfret);
            }
            else {
              if(l2_issued == max_l2_prefetches && ((l2_issued+llc_issued) < degree_l2)) { // TODO
                debug(stderr,"\tTRYING NEW PREFETCH: Base Address=0x%llx, Prefetch Address=0x%llx, into LLC\n",Data->DataAddr,(region_tag+prediction));
                pfret = l2_prefetch_line(0, Data->DataAddr , (region_tag + prediction),FILL_LLC);
                if(pfret) {
                  C(LLC_Prefetches_Issued);
                  llc_issued++;
#ifdef CHARACTERIZE
                  blocksPrefetched.push_back(region_tag+prediction);
#endif
                } else  {
                  C(LLC_Prefetches_Dropped);
                }
                debug(stderr,"  PREF returned to LLC: %d\n",pfret);
              }
            }
          } // end loop over patterns
            // store data at this cycle
#ifdef CHARACTERIZE
          PrefetchedInfoClass newPrefetchedInfo;
          newPrefetchedInfo.tempo_way = ghbHash;
          newPrefetchedInfo.pfetches = blocksPrefetched; // includes llc lines
                                                         // map it up
          curInfoPair.second = newPrefetchedInfo;
          curMap[monoMissCtr] = curInfoPair;
          big_pc_map[pc] = curMap;
#endif
        } // end if-PHT entry found
        } // end if new_generation
      }
    };

    SMS* sms;
    PrefetchBitsArray* pfArray;
    MissTimeHistoryBuffers<MTHB_DEPTH,MTHB_WIDTH>* lhbs;
    MissTimeHistoryBuffers<MTHB_DEPTH,GHB_WIDTH>* ghbs;
    std::list<CacheAddrInfo_t> prefetchesInCache;
    NextLinePrefetcher* nlp;
    SIGNED_COUNTER last_interval;

    int accessToWayDecoder(int accessFromGHB)
    {
      if(accessFromGHB <= 4) return 0;
      if(accessFromGHB <= 9) return 1;
      if(accessFromGHB <= 15) return 2;
      return 3;
    }

    void updatePrefetcher() 
    {
#ifdef THROTTLE_ACCURACY
      if (accuratePhase >= pollutionPhase) { // turn up prefetcher
        if(sms->degree_l2 < 8) {
          sms->degree_l2 += 1;
        }
        if (sms->degree_llc < 16) {
          sms->degree_llc+=1;
        } 
      } else {
        if(sms->degree_l2 > 1) {
          sms->degree_l2 -= 1;
        } 
        if (sms->degree_llc > 1) {
          sms->degree_llc -= 1;
        }  
      }

      float pacc = (float)accuratePhase / (float)totalPhase;
      float ppol = (float) pollutionPhase / (float) totalPhase;
      printf("Total in phase=%lld, Phase accuracy=%0.5f, phase pollution=%0.5f\n",totalPhase,pacc,ppol);
#endif
#ifdef THROTTLE_AMAT
      //printf("Cur AMAT: %0.5f, Last AMAT: %0.5f\n",accumulator,last_amat);
      if ( accumulator <= last_amat ) { // turn up prefetcher, decreased the AMAT
        if(sms->degree_l2 < 8) {
          sms->degree_l2++;
        }
        if (sms->degree_llc < 12) {
          sms->degree_llc++;
        } 
      } else { // turn down, increased the AMAT
        if(sms->degree_l2 > 1) {
          sms->degree_l2--;
        } 
        if (sms->degree_llc > 0) {
          sms->degree_llc--;
        }  
      }
#endif
      //printf("New SMS L2=%d, SMS LLC=%d\n",sms->degree_l2,sms->degree_llc);
      accuratePhase = 0;
      pollutionPhase = 0;
      totalPhase = 0;
      // reset accumulator, and accum_init, and num_demands
      last_amat = accumulator;
      //accumulator = 0;
      num_demands = 0;
      //accum_init = false;
    }


    void l2_prefetcher_initialize(int cpu_num)
    {
      printf("Temp0 Prefetcher\n");
      // you can inspect these knob values from your code to see which configuration you're runnig in
      printf("\nKnobs visible from prefetcher are: %d %d %d\n\n", knob_scramble_loads, knob_small_llc, knob_low_bandwidth);

      sms = new SMS();
      nlp = new NextLinePrefetcher();
      pfArray = new PrefetchBitsArray();
      mshrs = new MSHRs(L2_MSHR_COUNT); // these mshr's are for demands
                                        //prefetchMSHRs = new MSHRs(32); // these are not "real" mshrs, they are just used to track prefetching stats
      lhbs = new MissTimeHistoryBuffers<MTHB_DEPTH,MTHB_WIDTH>();
      ghbs = new MissTimeHistoryBuffers<MTHB_DEPTH,GHB_WIDTH>();
      lhbs->init();
      ghbs->init();
      last_interval = 0;
      monoMissCtr = 0;
      accuratePhase = 0;
      pollutionPhase = 0;
      totalPhase = 0;
    }

    void l2_prefetcher_operate(int cpu_num, unsigned long long int addr, unsigned long long int ip, int cache_hit)
    {
#ifdef GATHER_ACCESS_PATTERNS
      block_t tmp_pc;
      tmp_pc.addr = addr & (~(unsigned long long)(0xFFF)); // page addr
      tmp_pc.offset = addr & ((unsigned long long)(0xFFF));
      tmp_pc.offset /= 64;
      tmp_pc.hit = cache_hit;
      //tmp_pc.cycle = get_current_cycle(cpu_num);
      tmp_pc.cycle = monoMissCtr;
      if (map_pc.find(ip) == map_pc.end()) {
        vector<block_t> blockentry;
        blockentry.clear();
        map_pc[ip] = blockentry;
      }
      map_pc[ip].push_back(tmp_pc);

      block_t tmp_page;
      unsigned long long pageaddr = addr & (~(unsigned long long)(0xFFF)); // page addr
      tmp_page.addr = ip; // pc addr
      tmp_page.offset = addr & ((unsigned long long)(0xFFF));
      tmp_page.offset /= 64;
      tmp_page.hit = cache_hit;
      //tmp_page.cycle = get_current_cycle(cpu_num);
      tmp_page.cycle = monoMissCtr;
      if (map_page.find(pageaddr) == map_page.end()) {
        vector<block_t> blockentry;
        blockentry.clear();
        map_page[pageaddr] = blockentry;
      }
      map_page[pageaddr].push_back(tmp_page);
#endif
      // uncomment this line to see all the information available to make prefetch decisions
      trace(stderr,"ACCESS: 0x%llx 0x%llx %d %d %d\n", addr, ip, cache_hit, get_l2_mshr_occupancy(0),get_l2_read_queue_occupancy(0));

      // get cpu cycle
      SIGNED_COUNTER cycle = get_current_cycle(cpu_num);
      if(!cache_hit) { // miss
        monoMissCtr++;
      } 
#ifdef THROTTLE_ACCURACY
      if ( cycle > (last_interval + 1000000) ) {
        last_interval = cycle;
        // reset counters and update prefetch aggressiveness
        updatePrefetcher();
      }
#endif
#ifdef THROTTLE_AMAT
      if( num_demands > RECALIBRATE_DEMANDS ) {
        updatePrefetcher();
      }
#endif
      std::vector<int> waysToFetch;

      // hash into the local history buffers and update with the current miss latency delta
      lhbs->update((PC_t)ip,monoMissCtr);
      int LocalAv = lhbs->getFunctDelta((PC_t) ip); // the average access delta for THIS PC
      int locHash = accessToWayDecoder(LocalAv);

      int GlobalAv = ghbs->getFunct((PC_t)ip); // the average miss delta seen by the memory system
      int globHash = accessToWayDecoder(GlobalAv);
      switch(globHash) {
        case 0: 
          C(gHash_Short);
          break;
        case 1: 
          C(gHash_Mid);
          break;
        case 2: 
          C(gHash_Long);
          break;
        case 3: 
          C(gHash_VLong);
          break;
      }
      switch(locHash) {
        case 0: 
          C(lHash_Short);
          break;
        case 1: 
          C(lHash_Mid);
          break;
        case 2: 
          C(lHash_Long);
          break;
        case 3: 
          C(lHash_VLong);
          break;
      }

#ifdef GATHER_TEMPOS
      if(tempo_map.find(ip) == tempo_map.end()) { // not found
        TempoBuckets newBuckets;
        newBuckets.local_vec[locHash]++;
        newBuckets.global_vec[globHash]++;
        tempo_map[ip] = newBuckets;
      } else {
        TempoBuckets prevBuckets = tempo_map[ip];
        prevBuckets.local_vec[locHash]++;
        prevBuckets.global_vec[globHash]++;
        tempo_map[ip] = prevBuckets;
      }
#endif

      while( globHash < 4) { // Now... Were you rushing, or were you dragging?
                             // fetch this "way"
        debug(stderr,"adding %d to waysToFetch\n",globHash);
        waysToFetch.push_back(globHash++);
      }


      // characterization for how many PC's currently have outstanding requests
#ifdef CHARACTERIZE_OLD
      {
        // counter mask
        if ( cycle < (last_interval + 1000000) ) {
          SIGNED_COUNTER mapThis = mshrs->getNumUniquePCs();

          sms->outstandingIPs[last_interval] = mapThis;
          last_interval = cycle;
          //pc_counter_type outstandingIPs;
        }
      } 
#endif

      // create data struct
      PrefetchData_t fifteenLove(ip,addr,cache_hit,monoMissCtr);

      // per cycle we only have 1 incoming request (and this gets called AFTER THE L2 LATENCY)
      bool resetPFBitInMSHR = mshrs->resetPFBits(fifteenLove,cycle); // TODO: do this no matter what??
      if(!cache_hit) {
        bool mshrAllocated = mshrs->allocate(ip,(CacheAddr_t)addr,cycle,0,monoMissCtr);
        if (mshrAllocated) {
          C(L2_Total_Misses);
        }
        else {
          C(L2_Total_MSHR_Hits);
          assert(!cache_hit);
          //L1Data[i].hit = 1;
        }
        //debug(stderr,"cycle %lld: %llx missed on %llx.\n",cycle,fifteenLove.LastRequestAddr,fifteenLove.DataAddr);
        if (resetPFBitInMSHR) {
          C(L2_Total_Prefetch_Hits_MSHR);
          debug(stderr,"cycle %lld: %llx prefetch hit in MSHR on %llx.\n",cycle,fifteenLove.LastRequestAddr,fifteenLove.DataAddr);
          assert(!cache_hit); // should only happen when we miss
        }
      } else {
        C(L2_Total_Hits);
        debug(stderr,"cycle %lld: %llx hit line %llx.\n",cycle,fifteenLove.LastRequestAddr,fifteenLove.DataAddr);
        // check if this hit a prefetch
      }
      bool new_gen = sms->updateAGTs(cycle,&fifteenLove,mshrs,locHash);

      int mshr_occupancy = get_l2_mshr_occupancy(0);
      int MAX_L2_PREFETCHES = 12 - mshr_occupancy;
      MAX_L2_PREFETCHES = (MAX_L2_PREFETCHES >= sms->degree_l2) ? sms->degree_l2 : MAX_L2_PREFETCHES;
      sms->max_l2_prefetches = MAX_L2_PREFETCHES;
      sms->l2_issued = 0;
      sms->llc_issued = 0;

      std::vector<int>::iterator it;
      for(it = waysToFetch.begin(); it != waysToFetch.end(); it++) {
        //for(size_t i = 0;i<waysToFetch.size();i++) {
        int fetchMe = *it;
        sms->IssuePrefetches(cycle,&fifteenLove,fetchMe,new_gen);
      }
#ifdef NLP
      if(sms->l2_issued == 0 && sms->llc_issued == 0) {
        if(mshr_occupancy < 12 ) {
          // next line
          for(int i = 1;i<=1;i++) {
            int ret = l2_prefetch_line(0,addr,nlp->get_pf_addr(addr,i),FILL_L2);
            if(ret) {
              C(NLP_Prefetches_Issued);
              totalPhase++;
            } else {
              C(NLP_Prefetches_Dropped);
            }
          }
        }
      }
#endif
      }

      void l2_cache_fill(int cpu_num, unsigned long long int addr, int set, int way, int prefetch, unsigned long long int evicted_addr)
      {
        // uncomment this line to see the information available to you when there is a cache fill event
        trace(stderr,"LINE FILL: 0x%llx %d %d %d 0x%llx\n", addr, set, way, prefetch, evicted_addr);
        if(!prefetch) num_demands++;

        SIGNED_COUNTER cycle = get_current_cycle(cpu_num);
        std::pair<CacheAddrInfo_t,SIGNED_COUNTER> returnedPair = mshrs->deallocMSHR(addr,prefetch,cycle);
        CacheAddrInfo_t returnedEntry = returnedPair.first;
        SIGNED_COUNTER mshr_allocate_cycle = returnedPair.second;

        trace(stderr,"\tRETURN FROM MSHR: pc=0x%llx, addr=0x%llx, pf=%d, ref=%d, mc=%d\n",returnedEntry.ip,returnedEntry.LineAddr,returnedEntry.PrefetchBit,returnedEntry.referenced,returnedEntry.miss_counter);

        /* // old stats code
           if (pfArray->pfbits[set][way] == 1) { // this set and way were filled with a pf, but evicted before use
           pollutionPhase++;
           } 
           */

        PC_t returnedPC = returnedEntry.ip;
        int miss_delta = monoMissCtr - returnedEntry.miss_counter;
        SIGNED_COUNTER amat = cycle - (mshr_allocate_cycle);

        // hash into the global memory timer and update with the current return delta
        if(returnedPC != 0 && amat > 0) {
          ghbs->update(returnedPC,miss_delta);

          // UPDATE EMA ONLY WITH AMAT FOR DEMAND REQUESTS
          if(!accum_init) { // first time dont' multiply by 0
            accumulator = (double) amat;
            accum_init = true;
          } else {
            accumulator = (alpha * (double) amat)  + ((double)(1.0-alpha)*accumulator);
          }
          //cout << "AMAT EMA: " << accumulator << endl;
        } 
        /* // old code used to generate stats
         * else {
         returnedPair = prefetchMSHRs->deallocMSHR(addr,prefetch,cycle);
         returnedEntry = returnedPair.first;
         returnedPC = returnedEntry.ip;
         miss_delta = monoMissCtr - returnedEntry.miss_counter;
         if(!returnedPC) { C(MSHR_Zero_PC_Returns); }
         else { ghbs->update(returnedPC,miss_delta); }
         }
         std::list<CacheAddrInfo_t>::iterator it;
         if(prefetch) {
        // update the prefetch bits array
        pfArray->pfbits[set][way] = 1;
        debug(stderr,"Adding addr=%llx, prefetch=%d to list\n",addr,prefetch);
        CacheAddrInfo_t me(addr,0,prefetch,0);
        it = find(prefetchesInCache.begin(),prefetchesInCache.end(),me);
        if( it == prefetchesInCache.end() ) 
        prefetchesInCache.push_back(me);
        } else {
        pfArray->pfbits[set][way] = 0;
        }
        */

        CacheAddrInfo_t evict(evicted_addr,0,0,0);
        if(evicted_addr) {
          sms->checkEvictions(cycle,mshrs,evicted_addr);
        }
        /*
        //for(int i=0;i<prefetchesInCache.size();i++) {
        for(it = prefetchesInCache.begin(); it != prefetchesInCache.end();it++) {
        //CacheAddrInfo_t inCache = prefetchesInCache[i];
        CacheAddrInfo_t inCache = *it;
        debug(stderr,"checking prefetched line %llx in cache.\n",inCache.LineAddr);
        if(inCache == evict) {
        if(inCache.referenced == false) {
        C(L2_Prefetches_Pollution)
        // add this to the magic pollution set
        sms->pollutionAddrs.insert(evicted_addr);
        }
        // do not count referenced.true because I count it in l2_pref_operate
        prefetchesInCache.erase(it);
        return;
        }
        }
        */
      }

      void l2_prefetcher_heartbeat_stats(int cpu_num)
      {
        //TODO: print out some stats every heartbeat
      }

      void l2_prefetcher_warmup_stats(int cpu_num)
      {
        printf("Prefetcher warmup complete stats\n\n");
      }

      void l2_prefetcher_final_stats(int cpu_num)
      {
        printf("Prefetcher final stats\n");
#ifdef CHARACTERIZE_OLD
        sms->printCharacterization();
#endif
#ifdef CHARACTERIZE
        // print per pc data.
        printf("*************************************\n************************************* PC DATA NEXT LEVEL\n******************************\n*******************************\n");

        PCToCycleInfoMap::iterator PCMapIt;
        for(PCMapIt = big_pc_map.begin();PCMapIt != big_pc_map.end(); PCMapIt++) {
          CycleToInfoMap::iterator CycleMapIt;
          PC_t curPC = PCMapIt->first;
          CycleToInfoMap curInfoMap = PCMapIt->second;
          printf("xxxxx PC=0x%llx\n",curPC);
          for(CycleMapIt = curInfoMap.begin(); CycleMapIt != curInfoMap.end(); CycleMapIt++) {
            SIGNED_COUNTER cycle = CycleMapIt->first;
            CycleInfoPair curInfoPair = CycleMapIt->second;
            printf("\tCYCLE=%lld\n",cycle);
            std::bitset<64> patternToPrint = curInfoPair.first.pattern_trained;
            cout << "\t\tTRAINING: agt_way=" << curInfoPair.first.tempo_way << ", pattern=" << patternToPrint << ", to_pht=" << curInfoPair.first.to_pht << endl;
            //printf("\t\tTRAINING: agt_way=%d, pattern=%llu,to_pht=%d\n",curInfoPair.first.tempo_way,curInfoPair.first.pattern_trained,curInfoPair.first.to_pht);
            printf("\t\tPREFETCHES: pht_way=%d\n",curInfoPair.second.tempo_way);
            for(size_t i = 0;i< curInfoPair.second.pfetches.size();i++) {
              printf("\t\t\tADDR PFETCHED=0x%llx\n",curInfoPair.second.pfetches[i]);
            } // end loop over pfetch vector
          } // end loop over cycle info map
        } // end loop over big pc map
#endif

#ifdef GATHER_ACCESS_PATTERNS
        while (map_pc.size() > 0) {
          map<unsigned long long, vector<block_t> >::iterator top_it = map_pc.end();
          unsigned int top_size = 0;

          for (map<unsigned long long, vector<block_t> >::iterator it = map_pc.begin(); it != map_pc.end(); it++) {
            if (it->second.size() > top_size) {
              top_it = it;
              top_size = it->second.size();
            }
          }
          if (top_size < 12) break;

          printf("---PC 0x%016llX  NUMACCESSES %lu---\n", top_it->first, top_it->second.size());
          for (unsigned int i = 0; i < top_it->second.size(); i++) {
            printf("  PAGE 0x%016llX  CYCLE %016llu  OFFSETS", top_it->second[i].addr, top_it->second[i].cycle);
            for (int j = 0; j < 64; j++) {
              if (j == top_it->second[i].offset) {
                if (top_it->second[i].hit) {
                  printf(" 1");
                } else {
                  printf(" 0");
                }
              } else {
                printf(" _");
              }
            }
            printf("\n");
          }
          printf("\n");

          map_pc.erase(top_it);
        }

        while (map_page.size() > 0) {
          map<unsigned long long, vector<block_t> >::iterator top_it = map_page.end();
          unsigned int top_size = 0;

          for (map<unsigned long long, vector<block_t> >::iterator it = map_page.begin(); it != map_page.end(); it++) {
            if (it->second.size() > top_size) {
              top_it = it;
              top_size = it->second.size();
            }
          }
          if (top_size < 12) break;

          printf("---PAGE 0x%016llX  NUMACCESSES %lu---\n", top_it->first, top_it->second.size());
          for (unsigned int i = 0; i < top_it->second.size(); i++) {
            printf("  PC 0x%016llX  CYCLE %016llu  OFFSETS", top_it->second[i].addr, top_it->second[i].cycle);
            for (int j = 0; j < 64; j++) {
              if (j == top_it->second[i].offset) {
                if (top_it->second[i].hit) {
                  printf(" 1");
                } else {
                  printf(" 0");
                }
              } else {
                printf(" _");
              }
            }
            printf("\n");
          }
          printf("\n");

          map_page.erase(top_it);
        }
#endif
#ifdef GATHER_TEMPOS
        cout << "############################# TEMPO INFORMATION ########################"<< endl;
        PCToTempoMap::iterator pcTempoMapIt;
        for(pcTempoMapIt = tempo_map.begin() ; pcTempoMapIt != tempo_map.end();pcTempoMapIt++) {
          PC_t curPC = pcTempoMapIt->first;
          cout << "@@PC = 0x" << curPC << endl;
          TempoBuckets curBuckets = pcTempoMapIt->second;
          vector<int> local_tempo_buckets = curBuckets.local_vec;
          vector<int> global_tempo_buckets = curBuckets.global_vec;
#if 0
          uint64_t local_tempo_buckets[4];
          uint64_t global_tempo_buckets[4];
          for(int i = 0;i < 4;i++) {
            local_tempo_buckets[i] = 0;
            global_tempo_buckets[i] = 0;
          }
          for(level2It = curMap.begin(); level2It != curMap.end(); level2It++) {
            SIGNED_COUNTER cycle = level2It->first;
            TempoPair thisCyclesInfo = level2It->second;
            cout << "\tCYCLE: " << cycle << ", LTempo: " << thisCyclesInfo.lhb_tempo << ", GTempo: " << thisCyclesInfo.ghb_tempo << endl;
            local_tempo_buckets[thisCyclesInfo.lhb_tempo]++;
            global_tempo_buckets[thisCyclesInfo.ghb_tempo]++;
          }
#endif
          cout << "\tShort: Local=" << local_tempo_buckets[0] << ", Global=" << global_tempo_buckets[0] << endl;
          cout << "\tMed: Local=" << local_tempo_buckets[1] << ", Global=" << global_tempo_buckets[1] << endl;
          cout << "\tLong: Local=" << local_tempo_buckets[2] << ", Global=" << global_tempo_buckets[2] << endl;
          cout << "\tVLong: Local=" << local_tempo_buckets[3] << ", Global=" << global_tempo_buckets[3] << endl;
        }
#endif
        delete sms;
        delete mshrs;
        //delete prefetchMSHRs;
        delete lhbs;
        delete ghbs;
        delete nlp;
        delete pfArray;
      }

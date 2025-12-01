//ID 10
//Paper Title: An Optimized AMPM-based Prefetcher Combined with Configurable Block Size
//Author: Qi Jia, Maulik Bakulbhai Padia, Kashyap Amboju and Huiyang Zhou
//compile command: g++ -o dpc2sim prefetcher_file lib/dpc2sim.a
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define AMPM_PAGE_COUNT 64
#define UINT32 unsigned int
#define INT32 int
#define ADDRINT unsigned long long
#define FILL_L2 (1<<2)
#define FILL_LLC (1<<3)
#define CPT_SET 8
#define CPT_ASSOC 64

using namespace std;

extern "C" {
  int get_l2_mshr_occupancy(int);
  int get_l2_read_queue_occupancy(int);
  int l2_prefetch_line(int, long long int, long long int, int);
  void l2_prefetcher_initialize(int);
  void l2_prefetcher_operate(int, unsigned long long int, unsigned long long int, int);
  unsigned long long int get_current_cycle(int cpu_num);
  int l2_get_set(unsigned long long int addr);
  int l2_get_way(int cpu_num, unsigned long long int addr, int set);
  void l2_cache_fill(int cpu_num, unsigned long long int addr, int set, int way, int prefetch, unsigned long long int evicted_addr);
  void l2_prefetcher_heartbeat_stats(int cpu_num);
  void l2_prefetcher_warmup_stats(int cpu_num);
  void l2_prefetcher_final_stats(int cpu_num);
}
extern int knob_low_bandwidth;
extern int knob_small_llc;
extern int knob_scramble_loads;

int floorLog2(unsigned x) {
  assert(x > 0);

  int y = 0;

  if(x & 0xffff0000) { y+=16; x>>=16;}
  if(x & 0x0000ff00) { y+= 8; x>>= 8;}
  if(x & 0x000000f0) { y+= 4; x>>= 4;}
  if(x & 0x0000000c) { y+= 2; x>>= 2;}
  if(x & 0x00000002) { y+= 1;}

  return y;
}

class L2SamplerEntry {
  public:
    ADDRINT tag;
    bool prefetch;
    bool valid;
    L2SamplerEntry() {
      tag = 0;
      prefetch = false;
      valid = false;
    }
    void invalidate() {
      tag = 0;
      prefetch = false;
      valid = false;
    }
};

class L2Sampler
{
  public:
    L2SamplerEntry ** L2Model;

    L2Sampler(int set, int assoc) {
      L2Model = new L2SamplerEntry*[set];
      for(int i = 0; i < set; i++) {
        L2Model[i] = new L2SamplerEntry[assoc];
        for(int j = 0; j < assoc; j++) {
          L2Model[i][j].invalidate();
        }
      } 
    }
    void unSetPrefetchBit(ADDRINT _addr, int set, int way) {
      bool hit = (_addr == L2Model[set][way].tag && L2Model[set][way].valid);
      //assert(hit);
      if(hit)
        L2Model[set][way].prefetch = false;

    }

    bool getPrefetchBit(ADDRINT _addr, int set, int way) {
      bool hit = (_addr == L2Model[set][way].tag && L2Model[set][way].valid) ? true : false;
      //assert(hit);
      if(hit)
        return (L2Model[set][way].prefetch);
      else
        return false;
    }
    void updateSampler(ADDRINT _addr, ADDRINT _evicted_addr, bool _prefetch, int set, int way) {
      //printf("addr %llx access sampler, is prefetch %d\n", _addr, _prefetch);
      bool hit = (_addr == L2Model[set][way].tag && L2Model[set][way].valid) ? true : false;
      //assert(!hit);
      if(!hit || hit) {
        if(L2Model[set][way].valid) {
          // assert(L2Model[set][way].tag == _evicted_addr);
        }
        L2Model[set][way].tag = _addr;
        L2Model[set][way].prefetch = _prefetch;
        L2Model[set][way].valid = true;
      }
    }
    bool accessSampler(ADDRINT _addr, int set, int way) {
      bool hit = (_addr == L2Model[set][way].tag && L2Model[set][way].valid) ? true : false;
      return hit;
    }
    bool Hit(ADDRINT _addr) {
      bool hit = false;
      for(int i = 0; i < 256; i++) {
        for(int j = 0; j < 8; j++) {
          if(L2Model[i][j].valid && L2Model[i][j].tag == _addr) {
            hit = true;
            return hit;
          }
        }
      }
      return hit;
    }
    void invalidate(ADDRINT _addr, int set) {
      bool hit = false;
      for(int i = 0; i < 8; i++) {
        if(L2Model[set][i].valid && L2Model[set][i].tag == _addr) {
          hit = true;
          L2Model[set][i].valid = false;
          L2Model[set][i].tag = 0;
          L2Model[set][i].prefetch = false;       
          break;
        }
      }
      assert(hit);
    }

};


typedef struct ampm_page
{
  // page address
  unsigned long long int page;

  // The access map itself.
  // Each element is set when the corresponding cache line is accessed.
  // The whole structure is analyzed to make prefetching decisions.
  // While this is coded as an integer array, it is used conceptually as a single 64-bit vector.
  int access_map[64];

  // This map represents cache lines in this page that have already been prefetched.
  // We will only prefetch lines that haven't already been either demand accessed or prefetched.
  int pf_map[64];

  //This map record addr prefetched to LLC
  int pf_LLC_map[64];

  //This map record addr prefetched by large_blk mech, pay attention large_blk
  //can only prefetch into L3
  int pf_LB_map[64];

  //This map record addr prefetched by page_offset(PO)
  int pf_PO_map[64];

  // used for page replacement
  unsigned long long int lru;
} ampm_page_t;
//Qi: code page table entry
typedef struct cpt_entry_t {
  int lru;
  int counter;
  int offset;
  bool prefetch;
} cpt_entry;

//mshr entry
typedef struct mshr_t {
  ADDRINT addr;
  int lru;
} mshr;

//ATD entry
class monitor_entry {
  public:
    bool valid;
    ADDRINT tag;

    monitor_entry() : valid(false), tag(-1)
  {}

    void invalidate() {
      valid = false;
      tag = -1;
    }
};

//ATD set
class monitor_set {
  public:
    int assoc;
    int set;
    monitor_entry **entries;

    monitor_entry* findEntry(ADDRINT tag, int& _hit_assoc) {
      for(int i = 0; i < assoc; ++i) {
        if(entries[i]->tag == tag && entries[i]->valid) {
          _hit_assoc = i;
          return entries[i];
        }
      }
      return NULL;
    }

    void moveToHead(monitor_entry* entry) {
      if(entries[0] == entry) 
        return;

      int i = 0;
      monitor_entry *next = entry;

      do {
        assert(i < assoc);
        monitor_entry *tmp = entries[i];
        entries[i] = next;
        next = tmp;
        ++i;
      } while(next != entry);
    }

    void moveToTail(monitor_entry* entry) { 
      if(entries[assoc - 1] == entry)
        return;

      int i = assoc - 1;
      monitor_entry *next = entry;

      do {
        assert(i >= 0);
        monitor_entry *tmp = entries[i];
        entries[i] = next;
        next = tmp;
        --i;
      } while(next != entry);
    }

};

//cache ATD
class monitor {
  public:
    monitor_set *sets;
    monitor_entry *entries;
    int Lscore;
    int Mscore;
    int assoc;
    unsigned numSets;
    int blkSize;
    int tagShift;
    int setShift;
    unsigned int setMask;

    monitor(int _assoc, unsigned _total_numSets, unsigned _blk_size, int sample_ratio) {
      assoc = _assoc;
      numSets = (_total_numSets*64) / (sample_ratio * _blk_size);
      blkSize = _blk_size;
      setShift = floorLog2(_blk_size);
      tagShift = floorLog2((_total_numSets*64) / _blk_size) + setShift;
      setMask = ((_total_numSets*64) / _blk_size) - 1;

      sets = new monitor_set[numSets];
      entries = new monitor_entry[numSets * assoc];
      Lscore = 0;
      Mscore = 0;

      unsigned entryIndex = 0;
      for(unsigned i = 0; i < numSets; ++i) {
        sets[i].assoc = assoc;
        sets[i].entries = new monitor_entry*[assoc];
        sets[i].set = numSets * 8 * ((i * blkSize) / 512) + (i % (512 / blkSize));
        //printf("assign %d set\n", sets[i].set);
        for(int j = 0; j < assoc; ++j) {
          monitor_entry *entry = &entries[entryIndex];
          ++entryIndex;

          entry->invalidate();

          entry->tag = j;
          sets[i].entries[j] = entry;
        }
      }
    }

    ADDRINT extractTag(ADDRINT addr) const {
      return (addr >> tagShift);
    }

    int extractSet(ADDRINT addr) const {
      return ((addr >> setShift) & setMask);
    }

    void updateMonitor(ADDRINT addr) {
      int temp_set = extractSet(addr);
      ADDRINT temp_tag = extractTag(addr);
      int index = -1;
      for(unsigned i = 0; i < numSets; i++) {
        if(sets[i].set == temp_set) {
          index = i;
          break;
        }
      }

      if(index != -1) {
        int hit_assoc = -1;
        monitor_entry *entry = sets[index].findEntry(temp_tag, hit_assoc);

        if(entry != NULL) {
          assert(hit_assoc != -1);
          assert(hit_assoc < assoc);
          Lscore++;
          sets[index].moveToHead(entry);
        }
        else {
          monitor_entry *victim = sets[index].entries[assoc - 1];
          victim->invalidate();
          victim->tag = temp_tag;
          victim->valid = true;
          Mscore++;
          sets[index].moveToHead(victim);
        }
      }
    }

    void resetScore() {
      Lscore = Lscore / 2;
      Mscore = Mscore / 2;
    }

    int getScore() {
      return (Lscore - Mscore*blkSize/64);
    }


};


//block size predictor, new added
monitor **online_monitor;
int blkSize;
int PREFETCH_DEGREE;
ampm_page_t ampm_pages[AMPM_PAGE_COUNT];
L2Sampler *sampler;
int conflict_table[AMPM_PAGE_COUNT];
int useful_pref;
int total_pref;
int overlap_pref;
int ampm_useful_pref;
int ampm_total_pref;
int page_pref_degree;
int L2_access;
int L2_miss;
cpt_entry cpts[CPT_SET][CPT_ASSOC];
ADDRINT page_service[CPT_SET];
mshr pmshr[32];
int control_state;
int total_conflict;
bool conflict_mode;
bool disable_lb;
int param;

//used to align the sub block addr to the large block addr
ADDRINT alignment(ADDRINT addr) {
  return (addr & ~((ADDRINT)(blkSize - 1)));
}


void l2_prefetcher_initialize(int cpu_num)
{
  printf("AMPM Opt Prefetcher\n");
  // you can inspect these knob values from your code to see which configuration you're runnig in
  printf("Knobs visible from prefetcher: %d %d %d\n", knob_scramble_loads, knob_small_llc, knob_low_bandwidth);
  sampler = new L2Sampler(256,8);
  int i;
  for(i=0; i<AMPM_PAGE_COUNT; i++)
  {
    ampm_pages[i].page = 0;
    ampm_pages[i].lru = 0;
    conflict_table[i] = 0;

    int j;
    for(j=0; j<64; j++)
    {
      ampm_pages[i].access_map[j] = 0;
      ampm_pages[i].pf_map[j] = 0;
      ampm_pages[i].pf_LLC_map[j] = 0;
      ampm_pages[i].pf_LB_map[j] = 0;
      ampm_pages[i].pf_PO_map[j] = 0;
    }
  }
  //Qi: init cpt
  for(int j=0;j<CPT_SET;j++) {
    page_service[j] = 0;
    for(int m=0;m<CPT_ASSOC;m++) {
      cpts[j][m].offset = 0;
      cpts[j][m].counter = 0;
      cpts[j][m].prefetch = false;
      cpts[j][m].lru = j;
    }
  }
  //Qi: init mshr
  for(int i=0; i<32;i++) {
    pmshr[i].lru = i;
    pmshr[i].addr = 0;
  }
  online_monitor = new monitor*[4];
  for(int k = 0; k < 4; k++) {
    online_monitor[k] = new monitor(8, 256, 64*((int)(pow(2,k))), 16);
    //printf("finish init %d monitor\n", k);
  }
  blkSize = 64;
  useful_pref = 0;
  total_pref = 0;
  overlap_pref = 0;
  page_pref_degree = 2;
  L2_access = 0;
  L2_miss = 0;
  control_state = 0;
  total_conflict = 0;
  conflict_mode =false;
  disable_lb = false;
  ampm_total_pref = 0;
  ampm_useful_pref = 0;
  PREFETCH_DEGREE = 2;
  if(knob_small_llc) {
    param = 16;
  }
  else {
    param = 64;
  }
}

void l2_prefetcher_operate(int cpu_num, unsigned long long int addr, unsigned long long int ip, int cache_hit)
{
  // uncomment this line to see all the information available to make prefetch decisions
  //printf("(0x%llx 0x%llx %d %d %d) ", addr, ip, cache_hit, get_l2_read_queue_occupancy(0), get_l2_mshr_occupancy(0));
  //printf("cylce %llu, addr %llx, page %llx, offset %llx, mshr %d, read queue %d, cache %s\n", get_current_cycle(0), addr, (addr>>12)<<12, addr-((addr>>12)<<12),get_l2_mshr_occupancy(0), get_l2_read_queue_occupancy(0),cache_hit?"hit":"miss");
  L2_access++;

  //Qi: generate blkSize
  for(int i = 0; i < 4; i++) {
    online_monitor[i]->updateMonitor(addr);
  }
  if(L2_access % 8192 == 0) {
    // generate new blk size
    int max_score = online_monitor[0]->getScore();
    blkSize = 64;
    for(int i = 0; i < 4; i++) {
      if((online_monitor[i]->getScore()) > max_score) {
        max_score = online_monitor[i]->getScore();
        blkSize = 64 * ((int)(pow(2,i)));
      }
      online_monitor[i]->resetScore();
    }  
    //printf("new block size is %d\n", blkSize);
  }
  int current_pref_degree = 0;

  //Qi: check if we need to disable the large block size and detect if many bits in map has been set
  if(L2_access % 1024 == 0) {
    //check if miss num is too high
    if(L2_miss >= (param*1024/blkSize)) {
      disable_lb = true;
    }
    else {
      disable_lb = false;
    }
    L2_miss = 0;
    //Qi: check if some of the pages has been accessed/pref L2 for most of the part
    for(int i=0;i<AMPM_PAGE_COUNT;i++) {
      int access_num = 0;
      int pref_num = 0;
      for(int j=0;j<64;j++) {
        if(ampm_pages[i].access_map[j] == 1) {
          access_num++;
        }
        if(ampm_pages[i].pf_map[j] == 1) {
          pref_num++;
        }
      }
      if((access_num>=48 || pref_num>=48) && conflict_table[i]>=12) {
        ampm_pages[i].page = 0;
        ampm_pages[i].lru = 0;
        conflict_table[i] = 0;
        for(int j=0;j<64;j++) {
          ampm_pages[i].access_map[j] = 0;
          ampm_pages[i].pf_map[j] = 0;
          ampm_pages[i].pf_LLC_map[j] = 0;
          ampm_pages[i].pf_LB_map[j] = 0;
          ampm_pages[i].pf_PO_map[j] = 0;
        }
      }
    }
  }
  //Qi: update CPT with the offset
  int cpt_id = (addr>>12) % CPT_SET;
  int blk_offset = addr - ((addr>>12)<<12);
  bool cpt_hit = false;
  int cpt_hit_assoc = -1;
  //check whether the offset hit in cpt
  for(int i=0; i<CPT_ASSOC;i++) {
    if(cpts[cpt_id][i].offset == blk_offset) {
      cpt_hit = true;
      cpt_hit_assoc = i;
      cpts[cpt_id][i].counter = (cpts[cpt_id][i].counter == 31) ? 31 : (cpts[cpt_id][i].counter + 1);
      break;
    }
  }
  //if hit, update the LRU
  if(cpt_hit) {
    for(int i=0;i<CPT_ASSOC;i++) {
      if(cpts[cpt_id][i].lru < cpts[cpt_id][cpt_hit_assoc].lru) {
        cpts[cpt_id][i].lru++;
      }
    }
    cpts[cpt_id][cpt_hit_assoc].lru = 0;
  }
  //if miss, evict LRU
  else {
    //printf("offset %x miss in CPT set %d\n", blk_offset, cpt_id);
    for(int i=0; i<CPT_ASSOC;i++) {
      if(cpts[cpt_id][i].lru == CPT_ASSOC-1) {
        cpts[cpt_id][i].counter = 1;
        cpts[cpt_id][i].offset = blk_offset;
        cpts[cpt_id][i].prefetch = false;
        cpts[cpt_id][i].lru = 0;
      }
      else {
        cpts[cpt_id][i].lru++;
      }
    }
  }

  if(L2_access % 4096 == 0) {
    //check prefetch offset and decrease by 2
    for(int i=0; i<CPT_SET;i++) {
      //printf("set id %d\n",i);
      for(int j=0; j<CPT_ASSOC;j++) {
        cpts[i][j].counter /= 2;
        if(cpts[i][j].prefetch) {
          cpts[i][j].counter = (cpts[i][j].counter <= 2) ? 0 : (cpts[i][j].counter-2);
          cpts[i][j].prefetch = false;
        }
        //printf("id %d offset %x counter %d ... ", j, cpts[i][j].offset, cpts[i][j].counter);
      }
      //printf("\n");
    }
  }

  if(L2_access % 8192 == 0) {
    //printf("total pref %d, useful pref %d, overlap pref %d\n",total_pref,useful_pref,overlap_pref);
    if(page_pref_degree == 0) {
      control_state++;
      if(control_state == 6) {
        page_pref_degree++;
        control_state = 0;
      }
    }
    else if(total_pref > 0) {
      float accuracy = float(useful_pref)/float(total_pref);
      float overlap = float(overlap_pref)/float(useful_pref);
      if(accuracy < 0.6 || overlap > 0.8) {
        page_pref_degree /= 2;
      }
      else if(accuracy < 0.8 || overlap > 0.6) {
        page_pref_degree = page_pref_degree;
      }
      else {
        page_pref_degree = (page_pref_degree == 4) ? 4 : (page_pref_degree * 2);
      }
    }
    else {
      page_pref_degree = page_pref_degree;
    }
    total_pref /= 2;
    useful_pref /= 2;
    overlap_pref /= 2;
  }

  //check if it is late prefetch miss
  if(!cache_hit) {
    bool late_prefetch = false;
    for(int i=0; i<32; i++) {
      if(pmshr[i].addr == addr) {
        //it is late prefetch
        late_prefetch = true;
        break; 
      }
    }
    if(late_prefetch) {
      //printf("offset %x is late prefetch\n", blk_offset);
      for(int i=0;i<CPT_ASSOC;i++) {
        if(cpts[cpt_id][i].offset == blk_offset && cpts[cpt_id][i].prefetch) {
          cpts[cpt_id][i].counter = (cpts[cpt_id][i].counter <= 2) ? 0 : (cpts[cpt_id][i].counter-2);
          cpts[cpt_id][i].prefetch = false;
          break;
        }
      }
    }
  }

  //check whether our L2 model works correctly
  bool _cache_hit;
  if(l2_get_way(0, addr, l2_get_set(addr)) == -1) {
    _cache_hit = false;
  }
  else {
    _cache_hit = sampler->accessSampler(addr, l2_get_set(addr), l2_get_way(0, addr, l2_get_set(addr)));
  }
  // inclusion policy applied, invalidate corresponding entry in L2 model
  if(cache_hit == false && _cache_hit == true) {
    sampler->invalidate(addr, l2_get_set(addr));
  }
  bool isPrefetch = false;
  if(cache_hit) {
    isPrefetch = sampler->getPrefetchBit(addr, l2_get_set(addr), l2_get_way(0, addr, l2_get_set(addr)));
  }

  if(cache_hit && _cache_hit && isPrefetch) {
    //check if the hit comes from page prefetch, if yes increase counter;
    for(int i=0; i<CPT_ASSOC; i++) {
      if(cpts[cpt_id][i].offset == blk_offset && cpts[cpt_id][i].prefetch) {
        cpts[cpt_id][i].counter = (cpts[cpt_id][i].counter >= 29) ? 31 : (cpts[cpt_id][i].counter+2);
        cpts[cpt_id][i].prefetch = false;
        break;
      }
    }
    sampler->unSetPrefetchBit(addr, l2_get_set(addr), l2_get_way(0, addr, l2_get_set(addr)));
  }



  int prefetch_enable = false;
  if(!cache_hit || (cache_hit && isPrefetch)) {
    prefetch_enable = true;
  }
  //printf("cylce %llu, pc %llx, addr %llx, page %llx, offset %llx, mshr %d, read queue %d, cache %s, prefetch %d\n", get_current_cycle(0), ip, addr, (addr>>12)<<12, addr-((addr>>12)<<12),get_l2_mshr_occupancy(0), get_l2_read_queue_occupancy(0),cache_hit?"hit":"miss", isPrefetch);

  unsigned long long int cl_address = addr>>6;
  unsigned long long int page = cl_address>>6;
  unsigned long long int page_offset = cl_address&63;



  // check to see if we have a page hit
  bool page_hit = false;
  int page_index = -1;
  int i;
  for(i=0; i<AMPM_PAGE_COUNT; i++)
  {
    if(ampm_pages[i].page == page)
    {
      page_index = i;
      page_hit = true;
      break;
    }
  }

  if(page_index == -1)
  {
    // the page was not found, so we must replace an old page with this new page

    // find the oldest page
    int lru_index = 0;
    unsigned long long int lru_cycle = ampm_pages[lru_index].lru;
    int i;
    for(i=0; i<AMPM_PAGE_COUNT; i++)
    {
      if(ampm_pages[i].lru < lru_cycle)
      {
        lru_index = i;
        lru_cycle = ampm_pages[lru_index].lru;
      }
    }
    page_index = lru_index;
    //printf("page miss, evict %x\n", ampm_pages[page_index].page);
    // reset the oldest page
    ampm_pages[page_index].page = page;
    conflict_table[page_index] = 0;
    for(i=0; i<64; i++)
    {
      ampm_pages[page_index].access_map[i] = 0;
      ampm_pages[page_index].pf_map[i] = 0;
      ampm_pages[page_index].pf_LLC_map[i] = 0;
      ampm_pages[page_index].pf_LB_map[i] = 0;
      ampm_pages[page_index].pf_PO_map[i] = 0;
    }
  }

  else {
    //check how accurate the access_map bit is
    //if((cache_hit && (ampm_pages[page_index].access_map[page_offset]!=1 && ampm_pages[page_index].pf_map[page_offset]!=1)) || 
    if((!cache_hit && (ampm_pages[page_index].access_map[page_offset]==1 || ampm_pages[page_index].pf_map[page_offset]==1))){
      conflict_table[page_index]++;
    } 

    //Qi: check if it is a good prefetch
    if(cache_hit) {
      if((ampm_pages[page_index].pf_map[page_offset]==1 || ampm_pages[page_index].pf_LLC_map[page_offset]==1) && ampm_pages[page_index].access_map[page_offset] != 1) {
        ampm_useful_pref++;
      }
    }
  }

  // Qi: check here if the prefetch hit is init by page offset
  if(cache_hit) {
    if(ampm_pages[page_index].access_map[page_offset] != 1 && ampm_pages[page_index].pf_PO_map[page_offset] == 1) {
      //printf("useful pref, addr %llx, page %llx, offset %llx\n", addr, (addr>>12)<<12, page_offset<<6);
      useful_pref++;
      if(ampm_pages[page_index].pf_map[page_offset] == 1 || ampm_pages[page_index].pf_LLC_map[page_offset] == 1) {
        overlap_pref++;
        //printf("overlap pref, addr %llx, page %llx, offset %llx\n", addr, (addr>>12)<<12, page_offset<<6); 
      }
    }
  }

  // update LRU
  ampm_pages[page_index].lru = get_current_cycle(0);

  // mark the access map
  ampm_pages[page_index].access_map[page_offset] = 1;
  //pay attention! only prefetch when miss,also if miss too frequently, means we may prefetch too many blocks, disable
  if(!cache_hit && !disable_lb) {
    int prefetch_degree = blkSize / 64;
    ADDRINT align_addr = alignment(addr);
    for(int i = 0; i < prefetch_degree; i++) {     
      ADDRINT pref_addr = align_addr + 64 * i;
      int pf_map_index = (pref_addr-((addr>>12)<<12))>>6;
      if(ampm_pages[page_index].pf_map[pf_map_index]!=1 && ampm_pages[page_index].access_map[pf_map_index]!=1 && ampm_pages[page_index].pf_LLC_map[pf_map_index]!=1 && ampm_pages[page_index].pf_LB_map[pf_map_index]!=1){
        int suc;
        suc = l2_prefetch_line(0, addr, pref_addr, FILL_LLC);
        if(suc){
          ampm_pages[page_index].pf_LB_map[pf_map_index] = 1;
          //do not do this thing, this will make 2-level prefetch inefficient  
        }
      }
    }
  }



  //Qi: selet page_pref_degree most accessed offset from CPT, if we first see the page, which means the access map in the
  // ampm page should not have more than one 1

  if(page_service[cpt_id] != page && !page_hit && prefetch_enable){
    //printf("new page\n");
    int pref_offset[4] = {-1,-1,-1,-1};
    int pref_way[4] = {-1,-1,-1,-1};
    int pref_counters[4] = {5,5,5,5};
    for(int i=0; i<page_pref_degree; i++) {
      for(int j=0;j<CPT_ASSOC;j++) {
        if((i == 0) || ((i!=0) && (j!=pref_way[0]) && (j!=pref_way[1]) && (j!=pref_way[2]))) {
          if(cpts[cpt_id][j].counter > pref_counters[i]) {
            pref_offset[i] = cpts[cpt_id][j].offset;
            pref_way[i] = j;
            pref_counters[i] = cpts[cpt_id][j].counter; 
          }
        }
      }

      //prefetch if we have valid offset
      if(pref_offset[i] != -1) {
        int temp_ampm_index = pref_offset[i] >> 6;
        //check if pref_addr is already in pmshr
        bool pmshr_hit = false;
        int pmshr_hit_assoc = -1;
        for(int j=0; j<32; j++) {
          if(pmshr[j].addr == ((page<<12)+pref_offset[i])) {
            pmshr_hit = true;
            pmshr_hit_assoc = j;
            break;
          }
        }
        //if hit, update LRU 
        if(pmshr_hit) {
          for(int j=0;j<32;j++) {
            if(pmshr[j].lru < pmshr[pmshr_hit_assoc].lru) {
              pmshr[j].lru++;
            }
          }
          pmshr[pmshr_hit_assoc].lru = 0;
        }
        //if it has sit in pf_map, set pf_PO_map bit also
        if(ampm_pages[page_index].pf_map[temp_ampm_index] == 1) {
          ampm_pages[page_index].pf_PO_map[temp_ampm_index] = 1;
        }
        //only prefetch if it is not in pmshr
        if((ampm_pages[page_index].pf_map[temp_ampm_index] != 1) && (ampm_pages[page_index].access_map[temp_ampm_index] != 1) && (ampm_pages[page_index].pf_PO_map[temp_ampm_index] != 1) && !pmshr_hit) {
          int suc;
          suc = l2_prefetch_line(0,addr,(page<<12)+pref_offset[i], FILL_L2);
          if(suc) {
            //printf("new pref add %llx, offset %x\n", (page<<12) + pref_offset[i], pref_offset[i]);
            cpts[cpt_id][pref_way[i]].prefetch = true;
            ampm_pages[page_index].pf_PO_map[temp_ampm_index] = 1;
            total_pref++;
            current_pref_degree++;
            //if succeed update the pmshr
            for(int j=0;j<32;j++) {
              if(pmshr[j].lru != 31) {
                pmshr[j].lru++;
              }
              else{
                pmshr[j].lru = 0;
                pmshr[j].addr = (page<<12)+pref_offset[i];
              }
            }

            ampm_pages[page_index].pf_PO_map[temp_ampm_index] = 1;
          }
        }
      }
    }
    page_service[cpt_id] = page;
  }



  // positive prefetching, the first cycle only search the offset which has not been pref into
  // both LLC and L2.
  int count_prefetches = 0;
  for(i=1; i<=16; i++)
  {
    int check_index1 = page_offset - i;
    int check_index2 = page_offset - 2*i;
    int pf_index = page_offset + i;

    if(!prefetch_enable) {
      break;
    }

    if(check_index2 < 0)
    {
      break;
    }

    if(pf_index > 63)
    {
      break;
    }

    if(count_prefetches >= (PREFETCH_DEGREE - (current_pref_degree/2)))
    {
      break;
    }

    if(ampm_pages[page_index].access_map[pf_index] == 1)
    {
      // don't prefetch something that's already been demand accessed
      continue;
    }

    if(ampm_pages[page_index].pf_map[pf_index] == 1)
    {
      // don't prefetch something that's alrady been prefetched
      continue;
    }

    if(ampm_pages[page_index].pf_LLC_map[pf_index] == 1) {
      continue;
    }

    if((ampm_pages[page_index].access_map[check_index1]==1) && (ampm_pages[page_index].access_map[check_index2]==1))
    {
      // check if it has been prefetched by PO approach
      if(ampm_pages[page_index].pf_PO_map[pf_index] == 1) {
        ampm_pages[page_index].pf_map[pf_index] = 1;
        //printf("detect overlap, page %llx, offset %x\n",(addr>>12)<<12, pf_index<<6);
        continue;
      }
      unsigned long long int pf_address = (page<<12)+(pf_index<<6);
      //int pref_offset = pf_address - ((addr>>12)<<12);
      // check the MSHR occupancy to decide if we're going to prefetch to the L2 or LLC
      if(get_l2_mshr_occupancy(0) < 12 && ampm_pages[page_index].pf_LB_map[pf_index] != 1)
      {
        int suc;
        suc = l2_prefetch_line(0, addr, pf_address, FILL_L2);
        if(suc) {
          count_prefetches++;
          ampm_total_pref++;
          ampm_pages[page_index].pf_map[pf_index] = 1;
          //printf("pos pref addr %llx to L2, offset %x, mshr %d, read queue %d\n", pf_address, pf_index<<6, get_l2_mshr_occupancy(0), get_l2_read_queue_occupancy(0));
        }
      }
      else 
      {
        if(ampm_pages[page_index].pf_LB_map[pf_index] != 1) {
          int suc;
          suc = l2_prefetch_line(0, addr, pf_address, FILL_LLC);	      
          if(suc) {
            count_prefetches++;
            ampm_total_pref++;
            ampm_pages[page_index].pf_LLC_map[pf_index] = 1;     
            //printf("pos pref addr %llx to LLC, offset %x, mshr %d, read queue %d\n", pf_address, pf_index<<6, get_l2_mshr_occupancy(0), get_l2_read_queue_occupancy(0));
          }
        }
        else{
          ampm_total_pref++;
          ampm_pages[page_index].pf_LLC_map[pf_index] = 1;
        }
      }
    }
  }

  // then search positive to see that prefetched by LLC but have not been prefetched into L2
  for(i=1; i<=16; i++)
  {
    int check_index1 = page_offset - i;
    int check_index2 = page_offset - 2*i;
    int pf_index = page_offset + i;

    if(!prefetch_enable) {
      break;
    }

    if(check_index2 < 0)
    {
      break;
    }

    if(pf_index > 63)
    {
      break;
    }

    if(count_prefetches >= (PREFETCH_DEGREE - (current_pref_degree/2)))
    {
      break;
    }

    if(ampm_pages[page_index].access_map[pf_index] == 1)
    {
      // don't prefetch something that's already been demand accessed
      continue;
    }

    if(ampm_pages[page_index].pf_map[pf_index] == 1)
    {
      // don't prefetch something that's alrady been prefetched
      continue;
    }

    if((ampm_pages[page_index].access_map[check_index1]==1) && (ampm_pages[page_index].access_map[check_index2]==1))
    {
      // check if it has been prefetched by PO approach
      if(ampm_pages[page_index].pf_PO_map[pf_index] == 1) {
        ampm_pages[page_index].pf_map[pf_index] = 1;
        //printf("detect overlap, page %llx, offset %x\n",(addr>>12)<<12, pf_index<<6);
        continue;
      }
      unsigned long long int pf_address = (page<<12)+(pf_index<<6);
      //int pref_offset = pf_address - ((addr>>12)<<12);
      // check the MSHR occupancy to decide if we're going to prefetch to the L2 or LLC
      if(get_l2_mshr_occupancy(0) < 12)
      {
        int suc;
        suc = l2_prefetch_line(0, addr, pf_address, FILL_L2);
        if(suc) {
          ampm_total_pref++;
          ampm_pages[page_index].pf_map[pf_index] = 1;
          count_prefetches++;
          //printf("pos pref addr %llx to L2, offset %x, mshr %d, read queue %d\n", pf_address, pf_index<<6, get_l2_mshr_occupancy(0), get_l2_read_queue_occupancy(0));
        }
      }
    }
  }


  // negative prefetching
  count_prefetches = 0;
  for(i=1; i<=16; i++)
  {
    int check_index1 = page_offset + i;
    int check_index2 = page_offset + 2*i;
    int pf_index = page_offset - i;

    if(!prefetch_enable) {
      break;
    }

    if(check_index2 > 63)
    {
      break;
    }

    if(pf_index < 0)
    {
      break;
    }

    if(count_prefetches >= (PREFETCH_DEGREE-(current_pref_degree/2)))
    {
      break;
    }

    if(ampm_pages[page_index].access_map[pf_index] == 1)
    {
      // don't prefetch something that's already been demand accessed
      continue;
    }

    if(ampm_pages[page_index].pf_map[pf_index] == 1)
    {
      // don't prefetch something that's alrady been prefetched
      continue;
    }

    if(ampm_pages[page_index].pf_LLC_map[pf_index] == 1) {
      continue;
    }

    if((ampm_pages[page_index].access_map[check_index1]==1) && (ampm_pages[page_index].access_map[check_index2]==1))
    {
      if(ampm_pages[page_index].pf_PO_map[pf_index] == 1) {
        ampm_pages[page_index].pf_map[pf_index] = 1;
        //printf("detect overlap, page %llx, offset %x\n",(addr>>12)<<12, pf_index<<6);
        continue;
      }
      // we found the stride repeated twice, so issue a prefetch

      unsigned long long int pf_address = (page<<12)+(pf_index<<6);
      //int pref_offset = pf_address - ((addr>>12)<<12);
      // check the MSHR occupancy to decide if we're going to prefetch to the L2 or LLC
      if(get_l2_mshr_occupancy(0) < 12 && ampm_pages[page_index].pf_LB_map[pf_index] != 1)
      {
        int suc;
        suc = l2_prefetch_line(0, addr, pf_address, FILL_L2);
        if(suc) {
          ampm_total_pref++;
          ampm_pages[page_index].pf_map[pf_index] = 1;
          count_prefetches++;
        }
      }
      else 
      {
        if(ampm_pages[page_index].pf_LB_map[pf_index] != 1) {
          int suc;
          suc = l2_prefetch_line(0, addr, pf_address, FILL_LLC);	      
          if(suc) {
            ampm_total_pref++;
            count_prefetches++;
            ampm_pages[page_index].pf_LLC_map[pf_index] = 1;
          }
        }
        else {
          ampm_total_pref++;
          ampm_pages[page_index].pf_LLC_map[pf_index] = 1;
        }
      }
    }
  }

  for(i=1; i<=16; i++)
  {
    int check_index1 = page_offset + i;
    int check_index2 = page_offset + 2*i;
    int pf_index = page_offset - i;

    if(!prefetch_enable) {
      break;
    }

    if(check_index2 > 63)
    {
      break;
    }

    if(pf_index < 0)
    {
      break;
    }

    if(count_prefetches >= (PREFETCH_DEGREE-(current_pref_degree/2)))
    {
      break;
    }

    if(ampm_pages[page_index].access_map[pf_index] == 1)
    {
      // don't prefetch something that's already been demand accessed
      continue;
    }

    if(ampm_pages[page_index].pf_map[pf_index] == 1)
    {
      // don't prefetch something that's alrady been prefetched
      continue;
    }

    if((ampm_pages[page_index].access_map[check_index1]==1) && (ampm_pages[page_index].access_map[check_index2]==1))
    {
      if(ampm_pages[page_index].pf_PO_map[pf_index] == 1) {
        ampm_pages[page_index].pf_map[pf_index] = 1;
        //printf("detect overlap, page %llx, offset %x\n",(addr>>12)<<12, pf_index<<6);
        continue;
      }
      // we found the stride repeated twice, so issue a prefetch

      unsigned long long int pf_address = (page<<12)+(pf_index<<6);
      //int pref_offset = pf_address - ((addr>>12)<<12);
      // check the MSHR occupancy to decide if we're going to prefetch to the L2 or LLC
      if(get_l2_mshr_occupancy(0) < 12)
      {
        int suc;
        suc = l2_prefetch_line(0, addr, pf_address, FILL_L2);
        if(suc) {
          ampm_total_pref++;
          ampm_pages[page_index].pf_map[pf_index] = 1;
          count_prefetches++;
          //printf("neg pref addr %llx to L2, offset %x, mshr %d, read queue %d\n", pf_address, pf_index<<6, get_l2_mshr_occupancy(0), get_l2_read_queue_occupancy(0));
        }
      }
    }
  }

}

void l2_cache_fill(int cpu_num, unsigned long long int addr, int set, int way, int prefetch, unsigned long long int evicted_addr)
{
  // uncomment this line to see the information available to you when there is a cache fill event
  //printf("0x%llx %d %d %d 0x%llx\n", addr, set, way, prefetch, evicted_addr);
  //if(addr == 0x7b42d7600 || evicted_addr == 0x7b42d7600)
  //printf("insert addr %llx into set %d way %d, evict addr %llx, %s\n", addr, set, way, evicted_addr,prefetch?"prefetch":"demand");
  sampler->updateSampler(addr, evicted_addr, prefetch, set, way);
  //release the entry in pmshr if hit
  bool pmshr_hit = false;
  int pmshr_hit_assoc = -1;
  for(int j=0; j<32; j++) {
    if(pmshr[j].addr == addr) {
      pmshr_hit = true;
      pmshr_hit_assoc = j;
      break;
    }
  }
  if(pmshr_hit) {
    for(int j=0;j<32;j++) {
      if(pmshr[j].lru > pmshr[pmshr_hit_assoc].lru) {
        pmshr[j].lru--;
      }
    }
    pmshr[pmshr_hit_assoc].lru = 31;
    pmshr[pmshr_hit_assoc].addr = 0;
  }

}

void l2_prefetcher_heartbeat_stats(int cpu_num)
{
  printf("Prefetcher heartbeat stats\n");
}

void l2_prefetcher_warmup_stats(int cpu_num)
{
  printf("Prefetcher warmup complete stats\n\n");
}

void l2_prefetcher_final_stats(int cpu_num)
{
  printf("Prefetcher final stats\n");
}

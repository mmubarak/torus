#ifndef INC_torus_h
#define INC_torus_h

#include <ross.h>

// unit time is nano second
// assume 475MB/s bandwidth, it takes 64 ns to transfer 32Byte
// assume average packet size 32B, exponential dirstribution
// assume average processing time 1 us = 1000 ns, exponential dirstribution 
//#define PACKET_SIZE 256
// processing time of message at the MPI level-- from DCMF paper

//ARCH is set to 1 for BG/P and 2 for BG/Q
#define ARCH 2
#define HEADER_SIZE 8

#if ARCH == 1
  #define MEAN_PROCESS 750.0
  #define HOP_DELAY 46
  #define OVERHEADS 2000.0 /*MPI software overheads*/
// Total available tokens on a VC = VC buffer size / token size
  #define VC_SIZE 1024 /*Each VC has a specific number of tokens and each token is of 32 bytes */
  #define NUM_BUF_SLOTS VC_SIZE/TOKEN_SIZE
  #define BANDWIDTH 0.425 /*Link bandwidth*/
  #define N_dims 3
  #define PACKET_SIZE_LIMIT 256 /* maximum size of packet in bytes */
  static int       dim_length[] = {8,8,8};
#else
  #define MEAN_PROCESS 500.0
  #define HOP_DELAY 40 
/*somehow small message size on BG/Q takes more time than BG/P
The only thing I can suspect its because of the PAMI or MPI overheads on 
BG/Q are more than the BG/P, so I have adjusted the overheads for BG/Q accordingly*/
  #define OVERHEADS 4100.0 /*MPI software overheads*/
  #define N_dims 5
// Total available tokens on a VC = VC buffer size / token size
  #define BANDWIDTH 2.147 /*Link bandwidth*/
  #define PACKET_SIZE_LIMIT 512
  #define VC_SIZE 8192 /*Each VC has a specific number of tokens and each token is of 32 bytes */
  #define NUM_BUF_SLOTS VC_SIZE/TOKEN_SIZE
//  static int       dim_length[] = {8,4,4,4,4,4,2};
  static int dim_length[] = {8,4,4,4,2};
#endif

#define MPI_MESSAGE_LIMIT 50 /*Number of messages to be injected by each node */
#define NUM_VC 2
#define TOKEN_SIZE 32
#define PING_PONG 0 /*Set 1 for a ping pong test, 0 for a bisection test */

// finite buffer
//#define N_dims 3
#define TRACK 49
#define N_COLLECT_POINTS 20

#define TRACK_LP 0
#define DEBUG 1

#define NUM_ZONE_NODES 32

#define NUM_ROWS 32
#define NUM_COLS 32

//static dim_length[] = {8,8,8};
//static int       dim_length[] = {8, 8, 8};
//static int       dim_length[] = {64,64,64,64};
//static int       dim_length[] = {2,2,2,2,2,2,2,2,2,2};
//static int       dim_length[] = {8,8,8,8,8,8,8,8};
//static int       dim_length[] = {4,4};

typedef enum nodes_event_t nodes_event_t;
typedef struct nodes_state nodes_state;
typedef struct mpi_process mpi_process;
typedef struct nodes_message nodes_message;
typedef struct waiting_list waiting_list;

// Test RC code in serial mode
int g_test_rc = 0;

// Total number of nodes in torus, calculate in main
static int N_nodes = 1;
static int N_mpi_procs = 1;

enum nodes_event_t
{
  GENERATE = 1,
  WAIT,
  ARRIVAL, 
  SEND,
  PROCESS,
  CREDIT,
  MPI_SEND,
  MPI_RECV
};

enum traffic
{
  UNIFORM_RANDOM=1,
  DRAGONFLY_ZONES,
  TRANSPOSE,
  NEAREST_NEIGHBOR
};
struct mpi_process
{
 unsigned long long message_counter;
 tw_stime next_available_time;

// For dragonfly zones, there needs to be a zone id at the torus level so that we can identify which MPI process belongs to which zone
 int zone_id;

 /*For matrix transpose traffic, we have a row and col value for each MPI process so that the message can be sent to the corresponding transpose of 
 the MPI process */
 int row, col;
};

struct nodes_state
{
  unsigned long long packet_counter;            
  tw_stime next_available_time;                 
  tw_stime next_link_available_time[2*N_dims][NUM_VC]; 
  tw_stime next_credit_available_time[2*N_dims][NUM_VC];
  unsigned int buffer[2*N_dims][NUM_VC]; 
  int dim_position[N_dims];
  int neighbour_minus_lpID[N_dims];
  int neighbour_plus_lpID[N_dims];
  int N_wait_to_be_processed;
  int source_dim;
  int direction;

  //first element of linked list
  struct waiting_list * root;

  // pointer to the linked list
  struct waiting_list * ptr;
};

struct nodes_message
{
  tw_stime transmission_time;
  tw_stime travel_start_time;
  tw_stime saved_available_time;
  tw_stime saved_link_available_time[2*N_dims][NUM_VC];

  unsigned long long packet_ID;
  nodes_event_t	 type;

  int saved_source_dim;
  int saved_direction;
  int saved_vc;
  int dest[N_dims];

  tw_lpid dest_lp;
  tw_lpid origin_lp;
  tw_lpid sender_lp;

  int my_N_queue;
  int my_N_hop;
  int queueing_times;
  int source_dim;
  int source_direction;
  int next_stop;
  int packet_size;
  int count;
  int wait_type;
};

struct waiting_list
{
   int dim;
   int dir;
   int vc;
   nodes_message * packet;
   struct waiting_list * next;
   struct waiting_list * prev;
};
tw_stime         average_travel_time = 0;
tw_stime         total_time = 0;
tw_stime 	 total_lp_time = 0;
tw_stime         max_latency = 0;

static unsigned long long       N_finished_packets = 0;
static unsigned long long N_finished_msgs = 0;

static unsigned long long       N_finished_storage[N_COLLECT_POINTS];
static unsigned long long       N_dropped_storage[N_COLLECT_POINTS];
static unsigned long long       N_generated_storage[N_COLLECT_POINTS];

static unsigned long long       total_queue_length = 0;
static unsigned long long       queueing_times_sum = 0;
static unsigned long long       total_hops = 0;

static int       half_length[N_dims];

static int	 nlp_nodes_per_pe;
static int 	 nlp_mpi_procs_per_pe;
static int total_lps;

// run time arguments
static int	 opt_mem = 3000;
static long mpi_message_size = 256;
static int num_mpi_msgs = 50;
static int distance = 1;
static int mem_factor = 6;
static double MEAN_INTERVAL=10.0;
static int TRAFFIC = TRANSPOSE;

//static float MEAN_INTERVAL=10.0;
FILE *g_event_trace_file=NULL;
int g_enable_event_trace=1;
int num_buf_slots;
int num_packets;
int lp_hops = 0;
int num_buf_slots = 0;
int num_zones = 0;

float link_delay=0.0;
float head_delay=0.0;
float credit_delay = 0.0;

#endif

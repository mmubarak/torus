#ifndef INC_torus_h
#define INC_torus_h

#include <ross.h>

// unit time is nano second
// assume 475MB/s bandwidth, it takes 64 ns to transfer 32Byte
// assume average packet size 32B, exponential dirstribution
// assume average processing time 1 us = 1000 ns, exponential dirstribution 
//#define PACKET_SIZE 256
// processing time of message at the MPI level-- from DCMF paper
#define MEAN_PROCESS 750.0
#define MEAN_INTERVAL 1
#define MPI_MESSAGE_SIZE 256
#define MPI_MESSAGE_LIMIT 50 /*Number of messages to be injected by each node */
#define HOP_DELAY 900 /*Processing delay on each node */
#define PACKET_SIZE_LIMIT 256 /* maximum size of packet in bytes */
#define BANDWIDTH 0.374 /*Link bandwidth*/
#define OVERHEADS 3000.0 /*MPI software overheads*/
#define NUM_VC 2
// Total available tokens on a VC = VC buffer size / token size
#define NUM_BUF_SLOTS 1024/TOKEN_SIZE /*Each VC has a specific number of tokens and each token is of 32 bytes */
#define TOKEN_SIZE 32
#define PING_PONG 1 /*Set 1 for a ping pong test, 0 for a bisection test */

// finite buffer
#define N_dims 2
#define TRACK 3136
#define N_COLLECT_POINTS 20

#define TRACK_LP 0
#define DEBUG 1

//static int       dim_length[] = {4,4,4,4,2};
static int       dim_length[] = {8,8};
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

struct mpi_process
{
 unsigned long long message_counter;
 tw_stime next_available_time;
};

struct nodes_state
{
  unsigned long long packet_counter;            
  tw_stime next_available_time;                 
  tw_stime next_link_available_time[2*N_dims][NUM_VC]; 
  unsigned int buffer[2*N_dims][NUM_VC]; 
  int dim_position[N_dims];
  int neighbour_minus_lpID[N_dims];
  int neighbour_plus_lpID[N_dims];
  int N_wait_to_be_processed;
  int source_dim;
  int direction;
  int generate_counter;
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

  int sender_lp;
  int my_N_queue;
  int my_N_hop;
  int queueing_times;
  int source_dim;
  int source_direction;
  int next_stop;
  int packet_size;
  int count;
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

static int	 opt_mem = 3000;

tw_stime g_tw_last_event_ts = -1.0;
tw_lpid  g_tw_last_event_lpid = 0;
enum nodes_event_t g_tw_last_event_type=0;
FILE *g_event_trace_file=NULL;
int g_enable_event_trace=1;
int num_buf_slots;
int num_packets;

float link_delay=0.0;
float credit_delay = 0.0;

#endif

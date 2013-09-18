#include "torus.h"

tw_peid
mapping(tw_lpid gid)
{
	return (tw_peid) gid / g_tw_nlp;
}

void
torus_setup(nodes_state * s, tw_lp * lp)
{
  int i, j;
  int dim_N[N_dims+1];
  dim_N[0]=(int)lp->gid;

  s->waiting_messages=(nodes_message**)malloc(100*sizeof(nodes_message*));
  for(i=0; i<100; i++)
    s->waiting_messages[i] = (nodes_message*)sizeof(nodes_message);

  s->wait_indx=0;
  // calculate my torus co-ordinates
  for (i=0; i<N_dims; i++)
    {
      s->dim_position[i] = dim_N[i]%dim_length[i];
      dim_N[i+1] = ( dim_N[i] - s->dim_position[i] )/dim_length[i];

      half_length[i] = dim_length[i]/2;
    }

  for (i=0; i<N_dims; i++)
    {
     for(j=0; j<NUM_VC; j++)
      {	     
      	    s->node_queue_length[0+(i*2)][j] = 0;
            s->node_queue_length[1+(i*2)][j] = 0;
      }    
    }

  int factor[N_dims];
  factor[0]=1;
  for (i=1; i<N_dims; i++)
    {
      factor[i]=1;
      for (j=0; j<i; j++)
        factor[i]*=dim_length[j];
    }

  int temp_dim_pos[N_dims];
  for (i=0; i<N_dims; i++)
    temp_dim_pos[i] = s->dim_position[i];

  //printf("\n LP %d dimension: %d %d %d %d %d", (int)lp->gid, s->dim_position[0],s->dim_position[1], s->dim_position[2], s->dim_position[3], s->dim_position[4]);
  // calculate minus neighbour's lpID
  for (j=0; j<N_dims; j++)
    {
      temp_dim_pos[j] = (s->dim_position[j] -1 + dim_length[j])%dim_length[j];
      s->neighbour_minus_lpID[j]=0;
      for (i=0; i<N_dims; i++)
        s->neighbour_minus_lpID[j]+=factor[i]*temp_dim_pos[i];
      //printf(" minus neighbor LP ID is %d ", s->neighbour_minus_lpID[j]);
      temp_dim_pos[j] = s->dim_position[j];
    }

  // calculate plus neighbour's lpID
  for (j=0; j<N_dims; j++)
    {
      temp_dim_pos[j] = (s->dim_position[j] + 1 + dim_length[j])%dim_length[j];
      s->neighbour_plus_lpID[j]=0;
      for (i=0; i<N_dims; i++)
        s->neighbour_plus_lpID[j]+=factor[i]*temp_dim_pos[i];
      //printf(" plus neighbor LP ID is %d ", s->neighbour_plus_lpID[j]);
      temp_dim_pos[j] = s->dim_position[j];
    }
  // record LP time
  s->packet_counter = 0;
  s->next_available_time = 0;
  s->N_wait_to_be_processed = 0;
}


void
torus_init(nodes_state * s, tw_lp * lp)
{
  tw_event *e;
  tw_stime ts;
  nodes_message *m;

  /*
   * Set up the initial state for each LP (node)
   * according to the number of dimensions
   */
  torus_setup(s, lp);

  /*
   * Start a GENERATE event on each LP
   */
  ts = tw_rand_exponential(lp->rng, MEAN_INTERVAL);
  e = tw_event_new(lp->gid, ts, lp);
  m = tw_event_data(e);
  m->type = GENERATE;

  m->dest_lp = lp->gid;
  tw_event_send(e);

}


void
packet_arrive(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  int i;
  tw_stime ts;
  tw_event *e;
  nodes_message *m;

  bf->c2 = 1;
  bf->c6 = 1;
  bf->c10 = 1;

  if( msg->packet_ID == TRACK )
    {
      printf( "[LP %lld] packet %lld has arrived\n", 
	      lp->gid, msg->packet_ID );
      for( i = 0; i < N_dims; i++ )
	printf( "the %d dim position is %d\n",
		i,
		s->dim_position[i]);
      printf("packet %lld destination is\n",
	     msg->packet_ID);
      for( i = 0; i < N_dims; i++ )
	printf(" the %d dim position is %d\n", 
	       i,
	       msg->dest[i]);
      printf("lp time is %f travel start time is %f\n",
	     tw_now(lp),
	     msg->travel_start_time);
      printf("My hop now is %d\n",msg->my_N_hop);
      printf("\n");
    }

  // Packet arrives and accumulate # queued
  msg->queueing_times += s->N_wait_to_be_processed;
  // One more arrives and wait to be processed
  s->N_wait_to_be_processed++;

  msg->my_N_hop++;
  msg->my_N_queue+=s->node_queue_length[msg->source_direction][msg->source_dim];
    /*
     * If no packet is queueing, router is available now
     * if there are some packets queueing, then available time > tw_now(lp)
     * add one more MEAN_PROCESS to available time
     */

  msg->saved_available_time = s->next_available_time;
  s->next_available_time = max(s->next_available_time, tw_now(lp));
    // consider 1% noise on packet header parsing
  ts = tw_rand_exponential(lp->rng, (double)MEAN_PROCESS/1000)+MEAN_PROCESS;

  s->node_queue_length[msg->source_direction][msg->source_dim]++;

  e = tw_event_new(lp->gid, s->next_available_time - tw_now(lp), lp);

  s->next_available_time += ts;    

  m = tw_event_data(e);
  m->type = PROCESS;
  m->sender_lp = msg->sender_lp;
    
   //carry on the message info	
   for( i = 0; i < N_dims; i++ )
     m->dest[i] = msg->dest[i];
   m->dest_lp = msg->dest_lp;
   m->transmission_time = msg->transmission_time;
    
   m->source_dim = msg->source_dim;
   m->source_direction = msg->source_direction;
   m->saved_vc = msg->saved_vc;
    
   m->packet_ID = msg->packet_ID;	  
   m->travel_start_time = msg->travel_start_time;

   m->my_N_hop = msg->my_N_hop;
   m->my_N_queue = msg->my_N_queue;
   m->queueing_times = msg->queueing_times;
    
    tw_event_send(e);
}

void 
dimension_order_routing(nodes_state * s,nodes_message * msg, tw_lpid * dst_lp)
{
  int i;
  for( i = 0; i < N_dims; i++ )
    {
      if ( s->dim_position[i] - msg->dest[i] > half_length[i] )
	{
	  *dst_lp = s->neighbour_plus_lpID[i];
	  s->source_dim = i;
	  s->direction = 1;
	  break;
	}
      if ( s->dim_position[i] - msg->dest[i] < -half_length[i] )
	{
	  *dst_lp = s->neighbour_minus_lpID[i];
	  s->source_dim = i;
	  s->direction = 0;
	  break;
	}
      if (( s->dim_position[i] - msg->dest[i] <= half_length[i] )&&(s->dim_position[i] - msg->dest[i] > 0))
	{
	  *dst_lp = s->neighbour_minus_lpID[i];
	  s->source_dim = i;
	  s->direction = 0;
	  break;
	}
      if (( s->dim_position[i] - msg->dest[i] >= -half_length[i] )&&(s->dim_position[i] - msg->dest[i] < 0))
	{
	  *dst_lp = s->neighbour_plus_lpID[i];
	  s->source_dim = i;
	  s->direction = 1;
	  break;
	}
    }
}

void credit_send(nodes_state * s, tw_bf * bf, tw_lp * lp, nodes_message * msg, int saved_dim, int saved_dir)
{
  tw_event * buf_e;
  nodes_message *m;
  tw_stime ts;

  ts = LINK_DELAY;
  
  s->next_link_available_time[saved_dir+(saved_dim*2)][msg->saved_vc] =
	       max(s->next_link_available_time[saved_dir+(saved_dim*2)][msg->saved_vc], tw_now(lp));
  
  s->next_link_available_time[saved_dir+(saved_dim*2)][msg->saved_vc] += ts;

  buf_e = tw_event_new(msg->sender_lp, ts, lp);
  
  m = tw_event_data(buf_e);
  m->saved_vc = msg->saved_vc;
  m->source_direction=saved_dir;
  m->source_dim = saved_dim;

  m->type=CREDIT;
  tw_event_send(buf_e);
}

void
packet_send(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  int i;
  tw_lpid dst_lp;
  tw_stime ts;
  tw_event *e;
  nodes_message *m;
  /*
   * Routing in the torus, start from the first dimension
   */
  int saved_dir, saved_dim;
  saved_dim = msg->source_dim;
  saved_dir = msg->source_direction;

  if (lp->gid !=  msg->dest_lp )  
    {
      bf->c3 = 1;
      dimension_order_routing(s,msg,&dst_lp);
      msg->source_dim = s->source_dim;
      msg->source_direction = s->direction;
    }
  else
    {
      bf->c3 = 0;
      dst_lp = lp->gid;
    }
 
  bf->c2=1;

  int vc=0;
  if(lp->gid > dst_lp)
    vc=1;

  if(s->buffer[s->direction+(s->source_dim*2)][vc]<NUM_BUF_SLOTS)
  {
  if(msg->saved_vc != -1) 	  
    credit_send(s, bf, lp, msg, saved_dim, saved_dir);

  if(msg->packet_ID == TRACK)
	  printf("\n (%lf) [LP %d] Packet %lld being sent to destination %lld source dim %d source dir %d ", 
			  tw_now(lp), (int)lp->gid, msg->packet_ID,dst_lp, msg->source_dim, msg->source_direction);
   ts = LINK_DELAY+PACKET_SIZE;//msg->transmission_time;

   s->next_link_available_time[s->direction+(s->source_dim*2)][vc] = 
     max(s->next_link_available_time[s->direction+(s->source_dim*2)][vc], tw_now(lp));

   s->buffer[s->direction+(s->source_dim*2)][vc]++;
   // consider 1% noise on packet header parsing
   ts = tw_rand_exponential(lp->rng, (double)LINK_DELAY/1000)+
     LINK_DELAY+PACKET_SIZE;
  
    s->next_link_available_time[s->direction+(2*s->source_dim)][vc] += ts;      

    e = tw_event_new( dst_lp, 
		    s->next_link_available_time[s->direction+(2*s->source_dim)][vc] 
		    - tw_now(lp), 
		    lp);

   m = tw_event_data(e);
   m->type = ARRIVAL;
   m->sender_lp = lp->gid;
 
   // Carry on the message info
   m->source_dim = msg->source_dim;
   m->source_direction = msg->source_direction;
   m->saved_vc = vc;
  }
  else
  {
   // reschedule the packet send event?
   bf->c2=0;
  // ts = RESCHEDULE_DELAY + tw_rand_exponential(lp->rng, (double)RESCHEDULE_DELAY/1000);
  // e = tw_event_new(lp->gid, ts, lp);
  // m = tw_event_data(e);
  // m->type = SEND;
  // if(msg->packet_ID == TRACK)
    s->waiting_messages[s->wait_indx]=msg; 
    s->wait_indx++;
  }
  // Carry on the message info
   for( i = 0; i < N_dims; i++ )
      m->dest[i] = msg->dest[i];
   m->dest_lp = msg->dest_lp;
   m->transmission_time = msg->transmission_time;
  
   m->packet_ID = msg->packet_ID;
   m->travel_start_time = msg->travel_start_time;
  
   m->my_N_hop = msg->my_N_hop;
   m->my_N_queue = msg->my_N_queue;
   m->queueing_times = msg->queueing_times;
  
   tw_event_send(e);
}


void
packet_process(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  int i;
  tw_event *e;
  nodes_message *m;

  bf->c3 = 1;
  // One packet leaves the queue
  int vc=0; 
  if(lp->gid > msg->dest_lp)
     vc=1;

  s->node_queue_length[msg->source_direction+(2*msg->source_dim)][vc]--;
  s->N_wait_to_be_processed--;
  
  if(lp->gid==msg->dest_lp)
    {
      // one packet arrives and dies
      bf->c3 = 0;
	N_finished++;
	int index = floor(N_COLLECT_POINTS*(tw_now(lp)/g_tw_ts_end));
	N_finished_storage[index]++;
	total_time += tw_now(lp) - msg->travel_start_time;
	if (max_latency<tw_now(lp) - msg->travel_start_time)
	  max_latency=tw_now(lp) - msg->travel_start_time;
	total_hops += msg->my_N_hop;
	total_queue_length += msg->my_N_queue;
	queueing_times_sum += msg->queueing_times;
    }
  else
    {
      e = tw_event_new(lp->gid, MEAN_PROCESS, lp);
      m = tw_event_data(e);
      m->type = SEND;
      
      // Carry on the message info
      for( i = 0; i < N_dims; i++ )
	m->dest[i] = msg->dest[i];
      m->dest_lp = msg->dest_lp;
      m->transmission_time = msg->transmission_time;
      m->sender_lp = msg->sender_lp;
      m->saved_vc = msg->saved_vc; 
      
      m->source_dim = msg->source_dim;
      m->source_direction = msg->source_direction;
      
      m->packet_ID = msg->packet_ID;	  
      m->travel_start_time = msg->travel_start_time;

      m->my_N_hop = msg->my_N_hop;
      m->my_N_queue = msg->my_N_queue;
      m->queueing_times = msg->queueing_times;

      tw_event_send(e);
    }
}

void 
packet_generate(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  //unsigned long long i;
  int i;
  tw_lpid dst_lp;
  tw_stime ts;
  tw_event *e;
  nodes_message *m;

  // Send the packet out
  e = tw_event_new(lp->gid, 0.1, lp);
  m = tw_event_data(e);
  m->type = SEND;
  m->transmission_time = PACKET_SIZE;
  
  // Set up random destination
  dst_lp = tw_rand_integer(lp->rng,0,N_nodes-1);
  //rand_total += dst_lp;
  int dim_N[N_dims];
  dim_N[0]=dst_lp;
   
  // find destination dimensions using destination LP ID 
  for (i=0; i<N_dims; i++)
    {
      m->dest[i] = dim_N[i]%dim_length[i];
      dim_N[i+1] = ( dim_N[i] - m->dest[i] )/dim_length[i];
    }

  // record start time
  m->travel_start_time = tw_now(lp);
  m->my_N_queue = 0;
  m->my_N_hop = 0;
  m->queueing_times = 0;
  
  // set up packet ID
  // each packet has a unique ID
  m->packet_ID = (lp->gid * g_tw_nlp) + s->packet_counter;
  
  m->dest_lp = dst_lp;
  tw_event_send(e);	    

  if(m->packet_ID == TRACK)
	  printf("\n (%lf) Packet %lld generated destination LP %lld ", tw_now(lp), m->packet_ID, dst_lp);
  // One more packet is generating 
  s->packet_counter++;
  int index = floor(N_COLLECT_POINTS*(tw_now(lp)/g_tw_ts_end));
  N_generated_storage[index]++;

  // schedule next GENERATE event
  ts = tw_rand_exponential(lp->rng, MEAN_INTERVAL);	
  //ts = MEAN_INTERVAL;
  e = tw_event_new(lp->gid, ts, lp);
  m = tw_event_data(e);
  m->type = GENERATE;
  m->saved_vc = -1;
  
  m->dest_lp = lp->gid;
  tw_event_send(e);

}

void packet_buffer_process(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
 printf("\n Processing packet");	
 if(s->buffer[msg->source_direction+(msg->source_dim*2)][msg->saved_vc] > 0)	
  s->buffer[msg->source_direction+(msg->source_dim*2)][msg->saved_vc]--;

 // Send the message waiting in the queue
 if(s->wait_indx!=0)
 { 
  packet_send(s, bf, s->waiting_messages[0], lp); 

  // Re-arrange the array
  int i; 
  for(i=0; i<s->wait_indx-1; i++)
   s->waiting_messages[i]=s->waiting_messages[i+1];
  
  s->wait_indx--;
 }
}

void
event_handler(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  switch(msg->type)
    {
    case GENERATE:
      packet_generate(s,bf,msg,lp);
      break;
    case ARRIVAL:
      packet_arrive(s,bf,msg,lp);
      break;
    case SEND:
      packet_send(s,bf,msg,lp);
      break;
    case PROCESS:
      packet_process(s,bf,msg,lp);
      break;
    case CREDIT:
      packet_buffer_process(s,bf,msg,lp);
      break;  
    }
}

void
rc_event_handler(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  int index = floor(N_COLLECT_POINTS*(tw_now(lp)/g_tw_ts_end));
  switch(msg->type)
    {
    case GENERATE:
      N_generated_storage[index]--;
      s->packet_counter--;
      tw_rand_reverse_unif(lp->rng);
      tw_rand_reverse_unif(lp->rng);
      break;
    case ARRIVAL:
      s->next_available_time = msg->saved_available_time;
      tw_rand_reverse_unif(lp->rng);
      msg->my_N_hop--;
      s->N_wait_to_be_processed--;
      
      msg->queueing_times -= s->N_wait_to_be_processed;

      //s->node_queue_length[msg->source_direction+(NUM_VC*msg->source_dim)]--;
      //msg->my_N_queue-=s->node_queue_length[msg->source_direction+(NUM_VC*msg->source_dim)];
      break;
    case SEND:      
      s->source_dim = msg->saved_source_dim;
      s->direction = msg->saved_direction;
      //s->next_link_available_time[s->direction][s->source_dim]=      
	//msg->saved_link_available_time[s->direction][s->source_dim];
      tw_rand_reverse_unif(lp->rng);
      break;
    case PROCESS:
     if ( bf->c3 == 0 )
	{
	  N_finished--;
	  N_finished_storage[index]--;
	  total_time -= tw_now(lp) - msg->travel_start_time;
	  total_hops -= msg->my_N_hop;
	  total_queue_length -= msg->my_N_queue;
	  //event_queue_length -= msg->queue_length;
	  queueing_times_sum -= msg->queueing_times;
	}
      //s->node_queue_length[msg->source_direction+(NUM_VC*msg->source_dim)]++;
      s->N_wait_to_be_processed++;
      break;
    }
}

void
final(nodes_state * s, tw_lp * lp)
{
}

tw_lptype nodes_lps[] =
{
	{
		(init_f) torus_init,
		(event_f) event_handler,
		(revent_f) rc_event_handler,
		(final_f) final,
		(map_f) mapping,
		sizeof(nodes_state),
	},
	{0},
};

const tw_optdef app_opt [] =
{
	TWOPT_GROUP("Nodes Model"),
	TWOPT_UINT("memory", opt_mem, "optimistic memory"),
	TWOPT_STIME("arrive_rate", MEAN_INTERVAL, "packet arrive rate"),
	TWOPT_END()
};


int
main(int argc, char **argv, char **env)
{
	int i;
	tw_opt_add(app_opt);
	tw_init(&argc, &argv);

	for (i=0; i<N_dims; i++)
	  N_nodes*=dim_length[i];

	MEAN_INTERVAL = 100;


	nlp_per_pe = N_nodes/tw_nnodes()/g_tw_npe;
	g_tw_events_per_pe = (32 * nlp_per_pe/g_tw_npe) + opt_mem;
	tw_define_lps(nlp_per_pe, sizeof(nodes_message), 0);

	for(i = 0; i < g_tw_nlp; i++)
	  tw_lp_settype(i, &nodes_lps[0]);

	tw_run();

	if(tw_ismaster())
	{
		printf("\nTorus Network Model Statistics:\n");

		printf("\t%-50s %11lld\n", "Number of nodes", 
			nlp_per_pe * g_tw_npe * tw_nnodes());
	}

	unsigned long long total_finished_storage[N_COLLECT_POINTS];
	unsigned long long total_dropped_storage[N_COLLECT_POINTS];
	unsigned long long total_generated_storage[N_COLLECT_POINTS];
	unsigned long long wait_length,event_length,N_total_finish,N_total_hop;
	tw_stime total_time_sum,g_max_latency;

	for( i=0; i<N_COLLECT_POINTS; i++ )
	  {
	    MPI_Reduce( &N_dropped_storage[i], &total_dropped_storage[i],1, 
			MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce( &N_finished_storage[i], &total_finished_storage[i],1,
                        MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce( &N_generated_storage[i], &total_generated_storage[i],1,
                        MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	  }

	MPI_Reduce( &queueing_times_sum, &event_length,1, 
		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &total_queue_length, &wait_length,1, 
		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &total_time, &total_time_sum,1, 
		    MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &N_finished, &N_total_finish,1, 
		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &total_hops, &N_total_hop,1, 
		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &max_latency, &g_max_latency,1, 
		    MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

//	unsigned long long total_rand_total;
//	MPI_Reduce( &rand_total, &total_rand_total,1, 
//		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	for( i=1; i<N_COLLECT_POINTS; i++ )
	  {
	    total_dropped_storage[i]+=total_dropped_storage[i-1];
	    total_finished_storage[i]+=total_finished_storage[i-1];
	    total_generated_storage[i]+=total_generated_storage[i-1];
	  }

	if(tw_ismaster())
	  {
	    printf("\n ****************** \n");
	    printf("\n Total drop:           %lld; \n",
		   total_dropped_storage[N_COLLECT_POINTS-1]);
	    printf("\n total finish:         %lld and %lld; \n",
		   total_finished_storage[N_COLLECT_POINTS-1],N_total_finish);
	    printf("\n total generate:       %lld; \n",
		   total_generated_storage[N_COLLECT_POINTS-1]);
	    printf("\n total hops:           %lf; \n",
		   (double)N_total_hop/total_finished_storage[N_COLLECT_POINTS-1]);
	    printf("\n total wait length:    %lf; \n",
		   (double)wait_length/total_finished_storage[N_COLLECT_POINTS-1]);
	    printf("\n total total queued:   %lf; \n",
		   (double)event_length/total_finished_storage[N_COLLECT_POINTS-1]);
	    printf("\n average travel time:  %lf; \n\n",
		   total_time_sum/total_finished_storage[N_COLLECT_POINTS-1]);
	    
	    for( i=0; i<N_COLLECT_POINTS; i++ )
	      {
		printf(" %d ",i*100/N_COLLECT_POINTS);
		printf("drop: %lld; finish: %lld; generate: %lld; alive: %lld\n",
		       total_dropped_storage[i],
		       total_finished_storage[i],
		       total_generated_storage[i],
		       total_generated_storage[i]-total_finished_storage[i]);
		       
	      }

	    // capture the steady state statistics
	    unsigned long long steady_sum=0;
	    for( i = N_COLLECT_POINTS/2; i<N_COLLECT_POINTS;i++)
	      steady_sum+=total_generated_storage[i]-total_finished_storage[i];
	    printf("\n Steady state, packet alive: %lld\n",
		   2*steady_sum/N_COLLECT_POINTS);
//	    printf("Aeverage is %lld\n",total_rand_total/total_generated_storage[N_COLLECT_POINTS-1]);
	    printf("\nMax latency is %lf\n\n",g_max_latency);

	  }

	tw_end();
	return 0;
}

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

  // calculate my torus co-ordinates
  for (i=0; i<N_dims; i++)
    {
      s->dim_position[i] = dim_N[i]%dim_length[i];
      dim_N[i+1] = ( dim_N[i] - s->dim_position[i] )/dim_length[i];

      half_length[i] = dim_length[i]/2;
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

  // calculate minus neighbour's lpID
  for (j=0; j<N_dims; j++)
    {
      temp_dim_pos[j] = (s->dim_position[j] -1 + dim_length[j])%dim_length[j];
      s->neighbour_minus_lpID[j]=0;
      for (i=0; i<N_dims; i++)
        s->neighbour_minus_lpID[j]+=factor[i]*temp_dim_pos[i];
      temp_dim_pos[j] = s->dim_position[j];
    }

  // calculate plus neighbour's lpID
  for (j=0; j<N_dims; j++)
    {
      temp_dim_pos[j] = (s->dim_position[j] + 1 + dim_length[j])%dim_length[j];
      s->neighbour_plus_lpID[j]=0;
      for (i=0; i<N_dims; i++)
        s->neighbour_plus_lpID[j]+=factor[i]*temp_dim_pos[i];
      temp_dim_pos[j] = s->dim_position[j];
    }

  // initialize each node's waiting list
     s->root = NULL;
  
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

  torus_setup(s, lp);

  //Start a GENERATE event on each LP
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

  // One more arrives and wait to be processed
  s->N_wait_to_be_processed++;

  msg->my_N_hop++;
  msg->saved_available_time = s->next_available_time;
  s->next_available_time = max(s->next_available_time, tw_now(lp));
    // consider 1% noise on packet header parsing
  ts = tw_rand_exponential(lp->rng, (double)MEAN_PROCESS/1000)+MEAN_PROCESS;

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

void credit_send(nodes_state * s, tw_bf * bf, tw_lp * lp, nodes_message * msg)
{
  tw_event * buf_e;
  nodes_message *m;
  tw_stime ts;

  ts = LINK_DELAY;
  
//  s->next_link_available_time[msg->source_direction+(msg->source_dim*2)][msg->saved_vc] =
//	       max(s->next_link_available_time[msg->source_direction+(msg->source_dim*2)][msg->saved_vc], tw_now(lp));
  
//  s->next_link_available_time[msg->source_direction+(msg->source_dim*2)][msg->saved_vc] += ts;

  buf_e = tw_event_new(msg->sender_lp, ts, lp);
  
  m = tw_event_data(buf_e);
  m->saved_vc = msg->saved_vc;
  m->source_direction=msg->source_direction;
  m->source_dim = msg->source_dim;

  m->type=CREDIT;
  tw_event_send(buf_e);
}

void update_waiting_list(nodes_state * s, nodes_message * msg, tw_lp * lp, int dir, int dim, int vc)
{
   waiting_list * tmp = malloc(sizeof(waiting_list));
   waiting_list * tmp2 = malloc(sizeof(waiting_list));

   tmp->dim = dim;
   tmp->dir = dir;
   tmp->packet = malloc(sizeof(nodes_message));
   tmp->packet = msg;
   tmp->vc = vc;
   tmp->next = NULL;

   if(s->root==NULL)
   {
      // Insert at the root
      s->root = tmp;
   }
  else
   {
     tmp2 = s->root;
     // Traverse down to the end of the list
     while(tmp2->next != NULL)
      tmp2 = tmp2->next;
     
     // Append at the end of the list
      tmp2->next = tmp;        
   }
}

void packet_send(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  int i;
  tw_lpid dst_lp;
  tw_stime ts;
  tw_event *e;
  nodes_message *m;
  
   // Routing in the torus, start from the first dimension
  if (lp->gid !=  msg->dest_lp )  
      dimension_order_routing(s,msg,&dst_lp);
  else
      dst_lp = lp->gid;
 
  int vc=0;
  if(lp->gid > dst_lp)
    vc=1;
  
 if(s->buffer[s->direction+(s->source_dim*2)][vc]<NUM_BUF_SLOTS)
  {
  if(msg->saved_vc != -1) 	  
    credit_send(s, bf, lp, msg);

  if(msg->packet_ID == TRACK)
    printf("\n (%lf) [LP %d] Packet %lld being sent to destination %lld source dim %d source dir %d ", 
			  tw_now(lp), (int)lp->gid, msg->packet_ID,dst_lp, s->source_dim, s->direction);

   s->buffer[s->direction+(s->source_dim*2)][vc]++;

   s->next_link_available_time[s->direction+(s->source_dim*2)][vc] = 
     max(s->next_link_available_time[s->direction+(s->source_dim*2)][vc], tw_now(lp));

   ts = tw_rand_exponential(lp->rng, (double)LINK_DELAY/1000)+
     LINK_DELAY+PACKET_SIZE;
  
    s->next_link_available_time[s->direction+(2*s->source_dim)][vc] += ts;      

    e = tw_event_new( dst_lp, s->next_link_available_time[s->direction+(2*s->source_dim)][vc] - tw_now(lp), lp);

   m = tw_event_data(e);
   m->type = ARRIVAL;
   m->sender_lp = lp->gid;
 
   // Carry on the message info
   m->source_dim = s->source_dim;
   m->source_direction = s->direction;
   m->saved_vc = vc;
   for( i = 0; i < N_dims; i++ )
      m->dest[i] = msg->dest[i];
   m->dest_lp = msg->dest_lp;
   m->transmission_time = msg->transmission_time;
  
   m->packet_ID = msg->packet_ID;
   m->travel_start_time = msg->travel_start_time;
  
   m->my_N_hop = msg->my_N_hop;
   tw_event_send(e);
  }
  else
  {
   // Update the information in the linked list
   update_waiting_list(s, msg, lp, s->direction, s->source_dim, vc);
  }
}


void packet_process(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  int i;
  tw_event *e;
  nodes_message *m;

  s->N_wait_to_be_processed--;
  
  if(lp->gid==msg->dest_lp)
    {
      // one packet arrives and dies
	N_finished++;
	int index = floor(N_COLLECT_POINTS*(tw_now(lp)/g_tw_ts_end));
	N_finished_storage[index]++;
	total_time += tw_now(lp) - msg->travel_start_time;
	if (max_latency<tw_now(lp) - msg->travel_start_time)
	  max_latency=tw_now(lp) - msg->travel_start_time;
	total_hops += msg->my_N_hop;
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
      tw_event_send(e);
    }
}

void packet_generate(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
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
  m->saved_vc = -1;
  
  // Set up random destination
  dst_lp = tw_rand_integer(lp->rng,0,N_nodes-1);
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
  m->my_N_hop = 0;
  
  // set up packet ID
  // each packet has a unique ID
  m->packet_ID = lp->gid + (s->packet_counter * g_tw_nlp);
  
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
  e = tw_event_new(lp->gid, ts, lp);
  m = tw_event_data(e);
  m->type = GENERATE;
  tw_event_send(e);
}

void schedule_waiting_msg(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
 if(s->root == NULL)
   return;

 waiting_list * current = s->root;
 waiting_list * head = s->root;
 waiting_list * prev;
 
 int pos=0;
 while(current != NULL)
 {
  if(current->dim == msg->source_dim && current->dir == msg->source_direction && current->vc == msg->saved_vc)
   {
    pos = 1;
    packet_send(s, bf, current->packet, lp); 
    if(current==head)
      s->root = head->next;
    else
     prev->next = current->next;
   free(current);
   break;
   }
   else
   {
     prev = current;
     current = current->next;	
   }
 } 
}

void packet_buffer_process(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
 if(s->buffer[msg->source_direction+(msg->source_dim*2)][msg->saved_vc] > 0)
 {	
  s->buffer[msg->source_direction+(msg->source_dim*2)][msg->saved_vc]--;
 }

 // Send the message waiting in the queue
  schedule_waiting_msg(s, bf, msg, lp);
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
final(nodes_state * s, tw_lp * lp)
{
}

tw_lptype nodes_lps[] =
{
	{
		(init_f) torus_init,
		(event_f) event_handler,
		(revent_f) NULL,
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

	MEAN_INTERVAL = 10;

	nlp_per_pe = N_nodes/tw_nnodes()/g_tw_npe;
	g_tw_events_per_pe = (512 * nlp_per_pe/g_tw_npe) + opt_mem;
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
	unsigned long long total_generated_storage[N_COLLECT_POINTS];
	unsigned long long wait_length,event_length,N_total_finish,N_total_hop;
	tw_stime total_time_sum,g_max_latency;

	for( i=0; i<N_COLLECT_POINTS; i++ )
	  {
	    MPI_Reduce( &N_finished_storage[i], &total_finished_storage[i],1,
                        MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	    MPI_Reduce( &N_generated_storage[i], &total_generated_storage[i],1,
                        MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	  }

	MPI_Reduce( &total_time, &total_time_sum,1, 
		    MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &N_finished, &N_total_finish,1, 
		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &total_hops, &N_total_hop,1, 
		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce( &max_latency, &g_max_latency,1, 
		    MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	for( i=1; i<N_COLLECT_POINTS; i++ )
	  {
	    total_finished_storage[i]+=total_finished_storage[i-1];
	    total_generated_storage[i]+=total_generated_storage[i-1];
	  }

	if(tw_ismaster())
	  {
	    printf("\n ****************** \n");
	    printf("\n total finish:         %lld and %lld; \n",
		   total_finished_storage[N_COLLECT_POINTS-1],N_total_finish);
	    printf("\n total generate:       %lld; \n",
		   total_generated_storage[N_COLLECT_POINTS-1]);
	    printf("\n total hops:           %lf; \n",
		   (double)N_total_hop/total_finished_storage[N_COLLECT_POINTS-1]);
	    printf("\n average travel time:  %lf; \n\n",
		   total_time_sum/total_finished_storage[N_COLLECT_POINTS-1]);
	    
	    for( i=0; i<N_COLLECT_POINTS; i++ )
	      {
		printf(" %d ",i*100/N_COLLECT_POINTS);
		printf("finish: %lld; generate: %lld; alive: %lld\n",
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
	    printf("\nMax latency is %lf\n\n",g_max_latency);

	  }

	tw_end();
	return 0;
}

#include "torus.h"

//If the number of LPs is not evenly divisible the the number of processors, there is a remainder
int
getRem()
{
     return 0;
}
/*Get the MPI process ID mapped to the torus node ID*/
int 
getProcID( tw_lpid lpid )
{
	  return lpid - N_nodes;
}

/*Takes a MPI LP id and a torus node LP ID, returns the process ID on which the lp is mapped */
tw_peid 
mapping( tw_lpid gid )
{
    int rank;
    int rem = 0;

    if(gid < N_nodes)   {
         rank = gid / nlp_nodes_per_pe;
      }
    else {
         rank = getProcID( gid ) / nlp_mpi_procs_per_pe;
      }
    return rank;
}

/*Initialize the torus model, this initialization part is borrowed from Ning's torus model */
void
torus_init( nodes_state * s, 
	   tw_lp * lp )
{
    int i, j;
    int dim_N[ N_dims + 1 ];
    dim_N[ 0 ]=( int )lp->gid;

  // calculate my torus co-ordinates
  for ( i=0; i < N_dims; i++ ) 
    {
      s->dim_position[ i ] = dim_N[ i ]%dim_length[ i ];
      dim_N[ i + 1 ] = ( dim_N[ i ] - s->dim_position[ i ] )/dim_length[ i ];

      half_length[ i ] = dim_length[ i ] / 2;
    }

  int factor[ N_dims ];
  factor[ 0 ] = 1;
  for ( i=1; i < N_dims; i++ )
    {
      factor[ i ] = 1;
      for ( j = 0; j < i; j++ )
        factor[ i ] *= dim_length[ j ];
    }

  int temp_dim_pos[ N_dims ];
  for ( i = 0; i < N_dims; i++ )
    temp_dim_pos[ i ] = s->dim_position[ i ];

//if( lp->gid == TRACK_LP )
// printf("\n LP %d assigned dimensions %d %d %d %d %d", (int)lp->gid, s->dim_position[ 0 ], s->dim_position[ 1 ], s->dim_position[ 2 ], s->dim_position[ 3 ], s->dim_position[ 4 ]);

  // calculate minus neighbour's lpID
  for ( j = 0; j < N_dims; j++ )
    {
      temp_dim_pos[ j ] = (s->dim_position[ j ] -1 + dim_length[ j ]) % dim_length[ j ];

      s->neighbour_minus_lpID[ j ] = 0;
      
      for ( i = 0; i < N_dims; i++ )
        s->neighbour_minus_lpID[ j ] += factor[ i ] * temp_dim_pos[ i ];
     
      temp_dim_pos[ j ] = s->dim_position[ j ];

    }
   if( lp->gid == TRACK_LP )
   {
     for ( j = 0; j < N_dims; j++ )
       printf( " Neighbor %d LP ID %d ", j, s->neighbour_minus_lpID[ j ] );
   }
  // calculate plus neighbour's lpID
  for ( j = 0; j < N_dims; j++ )
    {
      temp_dim_pos[ j ] = ( s->dim_position[ j ] + 1 + dim_length[ j ]) % dim_length[ j ];

      s->neighbour_plus_lpID[ j ] = 0;
      
      for ( i = 0; i < N_dims; i++ )
        s->neighbour_plus_lpID[ j ] += factor[ i ] * temp_dim_pos[ i ];

      temp_dim_pos[ j ] = s->dim_position[ j ];
    }

  if( lp->gid == TRACK_LP ) {
     for ( j = 0; j < N_dims; j++ )
       printf( " Neighbor %d LP ID %d ", j, s->neighbour_plus_lpID[ j ]);
   }

  for( j=0; j < 2 * N_dims; j++ )
   {
    for( i = 0; i < NUM_VC; i++ )
     {
       s->buffer[ j ][ i ] = 0; //NUM_BUF_SLOTS * num_chunks;
       s->next_link_available_time[ j ][ i ] = 0.0;
       s->next_credit_available_time[j][i] = 0.0;
     }
   }
  // record LP time
    s->packet_counter = 0;
    s->waiting_list = tw_calloc(TW_LOC, "waiting list", sizeof(struct waiting_packet), WAITING_PACK_COUNT);
    
    for (j = 0; j < WAITING_PACK_COUNT - 1; j++) {
    s->waiting_list[j].next = &s->waiting_list[j + 1];
    s->waiting_list[j].dim = -1;
    s->waiting_list[j].dir = -1;
    s->waiting_list[j].packet = NULL; 
    }
    
    s->waiting_list[j].next = NULL;
    
    s->head = &s->waiting_list[0];
    s->wait_count = 0;
//    s->next_available_time = 0;
//    s->N_wait_to_be_processed = 0;

}
// initialize MPI process LP
void 
mpi_init( mpi_process * s, 
	 tw_lp * lp)
{
    tw_event *e;
    tw_stime ts;
    nodes_message *m;
    s->message_counter = 0;
//    s->next_available_time = 0;

    s->row = getProcID(lp->gid) / NUM_ROWS;
    s->col = getProcID(lp->gid) % NUM_COLS;

    //Start a GENERATE event on each LP
    ts =  tw_rand_exponential(lp->rng, MEAN_INTERVAL);
    e = tw_event_new( lp->gid, ts, lp );
    m = tw_event_data( e );
    m->type = MPI_SEND;

    s->available_time = 0.0;
//    s->next_available_time = 0;
    s->zone_id = getProcID(lp->gid) / NUM_ZONE_NODES;

//    printf("\n Zone ID for proc %d is %d", getProcID(lp->gid), s->zone_id);
    tw_event_send( e );
}

/*Returns the next neighbor to which the packet should be routed by using DOR (Taken from Ning's code of the torus model)*/
void 
dimension_order_routing( nodes_state * s,
			     tw_lpid * dst_lp, 
			     int * dim, 
			     int * dir )
{
  int dim_N[ N_dims ], 
      dest[ N_dims ],
      i;

  dim_N[ 0 ] = *dst_lp;

  // find destination dimensions using destination LP ID 
  for ( i = 0; i < N_dims; i++ )
    {
      dest[ i ] = dim_N[ i ] % dim_length[ i ];
      dim_N[ i + 1 ] = ( dim_N[ i ] - dest[ i ] ) / dim_length[ i ];
    }

  for( i = 0; i < N_dims; i++ )
    {
      if ( s->dim_position[ i ] - dest[ i ] > half_length[ i ] )
	{
	  *dst_lp = s->neighbour_plus_lpID[ i ];
	  *dim = i;
	  *dir = 1;
	  break;
	}
      if ( s->dim_position[ i ] - dest[ i ] < -half_length[ i ] )
	{
	  *dst_lp = s->neighbour_minus_lpID[ i ];
	  *dim = i;
	  *dir = 0;
	  break;
	}
      if ( ( s->dim_position[ i ] - dest[ i ] <= half_length[ i ] ) && ( s->dim_position[ i ] - dest[ i ] > 0 ) )
	{
	  *dst_lp = s->neighbour_minus_lpID[ i ];
	  *dim = i;
	  *dir = 0;
	  break;
	}
      if (( s->dim_position[ i ] - dest[ i ] >= -half_length[ i ] ) && ( s->dim_position[ i ] - dest[ i ] < 0) )
	{
	  *dst_lp = s->neighbour_plus_lpID[ i ];
	  *dim = i;
	  *dir = 1;
	  break;
	}
    }
}
/*Generates a packet. If there are two buffer slots available, then the packet is 
injected in the network. Else, the packet is placed in the injection queue */
void 
packet_generate( nodes_state * s, 
		tw_bf * bf, 
		nodes_message * msg, 
		tw_lp * lp )
{
    int i, j, tmp_dir=-1, tmp_dim=-1;
    tw_stime ts;
    bf->c3 = 1;
    bf->c2 = 1;

//    event triggered when packet head is sent
    tw_event * e_h;
    nodes_message *m;

    if(TRAFFIC == NEAREST_NEIGHBOR)
	msg->dest_lp = s->neighbour_minus_lpID[0]; 
  
//    Nodes at the diagnol will not participate in the send/receive operation 
    if(TRAFFIC == TRANSPOSE && msg->dest_lp == -1)
      {
        bf->c2 = 0;
	return;
      }
	
    tw_lpid dst_lp = msg->dest_lp; //s->neighbour_plus_lpID[0];
    dimension_order_routing( s, &dst_lp, &tmp_dim, &tmp_dir );
    ts = 0.001 + tw_rand_exponential(lp->rng, MEAN_INTERVAL/10000);

  if(tmp_dir == -1 || tmp_dim == -1)
    printf("\n LP %d dest LP %d dim %d dir %d ", (int)lp->gid, (int)msg->dest_lp, tmp_dim, tmp_dir);
    for(j = 0; j < num_chunks; j++)
    { 
       e_h = tw_event_new( lp->gid, j + ts, lp );
       m = tw_event_data( e_h );
       m->next_stop = dst_lp;
       m->dest_lp = msg->dest_lp;
       m->count = msg->count;
       m->travel_start_time = msg->travel_start_time;
       m->packet_ID = msg->packet_ID;
       m->chunk_id = j;
       m->sender_lp = -1;

       int dim_N[ N_dims ];
       dim_N[ 0 ] = m->dest_lp;
   
      // find destination dimensions using destination LP ID 
       for (i=0; i < N_dims; i++)
        {
           m->dest[ i ] = dim_N[ i ] % dim_length[ i ];
           dim_N[ i + 1 ] = ( dim_N[ i ] - m->dest[ i ] ) / dim_length[ i ];
        }

//       For reverse computation  17-05 placed these two lindes out of the if condition below
       if(s->buffer[ tmp_dir + ( tmp_dim * 2 ) ][ 0 ] < NUM_BUF_SLOTS * num_chunks)
        {

	 m->my_N_hop = 0;
	 m->wait_type = -1;
	    // Send the packet out
	 m->type = SEND;
         m->source_direction = tmp_dir;
         m->source_dim = tmp_dim;
#if DEBUG
          if( msg->packet_ID == TRACK )
          {
         //printf("\n (%lf) msg ID %lld Generating packet dimension ", m->packet_ID, tw_now(lp));
 	   for( i = 0; i < N_dims; i++ )
	     printf(" %d ", s->dim_position[ i ]);
   
	   printf("\n");

	   for( i = 0; i < N_dims; i++ )
	     printf(" %d ", m->dest[ i ]);
         }
#endif
#if DEBUG
//if(lp->gid == TRACK_LP)
//   printf("\n (%lld) Packet generated %lld Buffer space %d tmp_dir %d tmp_dim %d num_chunks %d dest_lp %lld", lp->gid, m->packet_ID, s->buffer[ tmp_dir + ( tmp_dim * 2 ) ][ 0 ], tmp_dir, tmp_dim, num_chunks, msg->dest_lp );
#endif
        }
      else 
       {
#if DEBUG
//if(lp->gid == TRACK_LP)
//   printf("\n %d Packet queued in line, buffer space %d ", (int)lp->gid, s->buffer[ tmp_dir + ( tmp_dim * 2 ) ][ 0 ]);
#endif
            m->wait_type = GENERATE;
	    m->type = WAIT;		
	    m->wait_dir = tmp_dir;	
	    m->wait_dim = tmp_dim;
       }
       tw_event_send(e_h);
   }
}
/*Sends a 8-byte credit back to the torus node LP that sent the message */
void 
credit_send( nodes_state * s, 
	    tw_bf * bf, 
	    tw_lp * lp, 
	    nodes_message * msg)
{
#if DEBUG
//if(lp->gid == TRACK_LP)
//	printf("\n (%lf) sending credit tmp_dir %d tmp_dim %d %lf ", tw_now(lp), msg->source_direction, msg->source_dim, credit_delay );
#endif
    tw_event * buf_e;
    nodes_message *m;
    tw_stime ts;
//    tw_stime ts;
    ts =  tw_rand_exponential(lp->rng, credit_delay/1000) + credit_delay;
  
    //buf_e = tw_event_new( lp->gid, ts, lp );
 //   s->next_credit_available_time[2 * msg->source_dim + msg->source_direction][0] = max(s->next_credit_available_time[2 * msg->source_dim + msg->source_direction][0], tw_now(lp) );
 //   s->next_credit_available_time[2 * msg->source_dim + msg->source_direction][0] += 2.0;

//    buf_e = tw_event_new( msg->sender_lp, s->next_credit_available_time[2 * msg->source_dim + msg->source_direction][0] - tw_now(lp), lp );
    buf_e = tw_event_new( msg->sender_lp, ts, lp);

    m = tw_event_data(buf_e);
    m->source_direction = msg->source_direction;
    m->source_dim = msg->source_dim;

    m->type = CREDIT;
    tw_event_send( buf_e );
}
/*Inserts a packet in the injection queue which then waits to be sent over the network */
void 
update_waiting_list( nodes_state * s, 
		     nodes_message * msg, 
		     tw_lp * lp )
{
   int loc = s->wait_count;

   if(loc >= WAITING_PACK_COUNT)
	printf(" Reached maximum count of linked list %d ", s->wait_count);

   s->waiting_list[loc].dim = msg->wait_dim;
   s->waiting_list[loc].dir = msg->wait_dir;
   s->waiting_list[loc].packet = msg;

   s->wait_count++;
//   printf("\n Inserted message, wait count %d dim %d dir %d type %d ", s->wait_count, msg->source_dim, msg->source_direction, msg->wait_type);
}
// send a packet from one torus node to another torus node
// A packet can be up to 256 bytes on BG/L and BG/P and up to 512 bytes on BG/Q
void 
packet_send( nodes_state * s, 
	         tw_bf * bf, 
		 nodes_message * msg, 
		 tw_lp * lp )
{   
#if DEBUG
// if( lp->gid == TRACK_LP )
//	printf("\n Sending credit %lf ", tw_now(lp) );
#endif 

    int i, vc = 0, tmp_dir, tmp_dim;
    tw_stime ts;
    tw_event *e;
    nodes_message *m;
    tw_lpid dst_lp = msg->dest_lp;
    int tokens_min = 0;
  
    // For reverse computation
    bf->c3 = 1;
    bf->c1 = 1;
    bf->c2 = 1;

//   If the message is coming from waiting queue then it will have everything set     
//    if(msg->wait_type == GENERATE)
//	printf("\n Sending waiting message!");

    if( msg->next_stop == -1 )  
     {
 	 dimension_order_routing( s, &dst_lp, &tmp_dim, &tmp_dir );     
     }
    else
      {
	 dst_lp = msg->next_stop;  
         tmp_dim = msg->wait_dim;
         tmp_dir = msg->wait_dir;
      }
     if(msg->source_dim != tmp_dim)
      {
          tokens_min = (2 * PACKET_SIZE) / TOKEN_SIZE;
      } 
      else
      {
         tokens_min = PACKET_SIZE/TOKEN_SIZE;
      }


    if(s->buffer[ tmp_dir + ( tmp_dim * 2 ) ][ 0 ] >= NUM_BUF_SLOTS * num_chunks )
    {
         // re-schedule the message in the future
         bf->c3 = 0;
	 ts = 0.001 + tw_rand_exponential( lp->rng, MEAN_INTERVAL/10000);
	 e = tw_event_new( lp->gid, ts, lp );
	 m = tw_event_data( e );	
         m->wait_type = SEND;
	 m->type = WAIT;
  	 m->wait_dim = tmp_dim;
         m->next_stop = dst_lp;
	 m->chunk_id = msg->chunk_id;
         m->dest_lp = msg->dest_lp;
         m->count = msg->count;
         m->wait_dir = tmp_dir;
         
	 // added may 20
	 m->source_direction = msg->source_direction;
	 m->source_dim = msg->source_dim;

	 m->packet_ID = msg->packet_ID;
	 m->sender_lp = msg->sender_lp;

	 for (i=0; i < N_dims; i++)
           m->dest[i] = msg->dest[i];

	 tw_event_send(e);

         return;
     }

    msg->saved_src_dir = tmp_dir;
    msg->saved_src_dim = tmp_dim;
    ts = tw_rand_exponential( lp->rng, ( double )head_delay/10000 )+
                              head_delay;

//    For reverse computation 
    msg->saved_available_time = s->next_link_available_time[tmp_dir + ( tmp_dim * 2 )][0];

    s->next_link_available_time[tmp_dir + ( tmp_dim * 2 )][0] = max( s->next_link_available_time[ tmp_dir + ( tmp_dim * 2 )][0], tw_now(lp) );
    s->next_link_available_time[tmp_dir + ( tmp_dim * 2 )][0] += ts;
    
    e = tw_event_new( dst_lp, s->next_link_available_time[tmp_dir + ( tmp_dim * 2 )][0] - tw_now(lp), lp );
    
    if(msg->sender_lp != -1)
      {
	  bf->c2 = 0;	
          credit_send( s, bf, lp, msg );
      }
       
   if( msg->packet_ID == TRACK )
      printf("\n (%lf) Scheduling for next hop after %f ", tw_now( lp ), head_delay);

    //if(msg->packet_ID == TRACK)
   //printf("\n (%lf) [LP %d] Packet %lld being sent to destination %lld source dim %d source dir %d Link delay %f ", 
               //			  tw_now(lp), (int)lp->gid, msg->packet_ID,msg->dest_lp, s->source_dim, s->direction, link_delay);
    m = tw_event_data( e );
    m->type = ARRIVAL;
 
    //Carry on the message info
    m->source_dim = tmp_dim;
    m->source_direction = tmp_dir;
    m->next_stop = dst_lp;
//    m->origin_lp = msg->origin_lp;
    m->sender_lp = lp->gid;
    m->chunk_id = msg->chunk_id;
   
    for( i = 0; i < N_dims; i++ )
       m->dest[ i ] = msg->dest[ i ];
     
    m->dest_lp = msg->dest_lp;
  
    m->packet_ID = msg->packet_ID;
    m->travel_start_time = msg->travel_start_time;
    m->count = msg->count;
  
    m->my_N_hop = msg->my_N_hop;
    tw_event_send( e );

    s->buffer[ tmp_dir + ( tmp_dim * 2 ) ][ 0 ]++;
    
    if(msg->chunk_id == num_chunks - 1 && msg->sender_lp == -1)
    {
       bf->c1 = 0;
       int index = floor( N_COLLECT_POINTS * ( tw_now( lp ) / g_tw_ts_end ) );
       N_generated_storage[ index ]++;           
   }
}
/*Once a credit arrives at the node, this method picks a waiting packet in the injection queue and schedules it */
void
waiting_packet_free(nodes_state * s, int loc)
{
  int i, max_count = s->wait_count - 1;

  for(i = loc; i < max_count; i++) {
	s->waiting_list[i].dim = s->waiting_list[i + 1].dim;
	s->waiting_list[i].dir = s->waiting_list[i + 1].dir;
	s->waiting_list[i].packet = s->waiting_list[i + 1].packet;
    }
  s->waiting_list[max_count].dim = -1;
  s->waiting_list[max_count].dir = -1;
  s->waiting_list[max_count].packet = NULL; 
  
  s->wait_count++; 
  if(s->wait_count < 0)
	printf("\n After scheduling, the wait count is less than zero");
}

void 
copyData(nodes_message * from, nodes_message * to)
{
   int i;

   to->wait_loc = from->wait_loc;
   to->source_dim = from->source_dim;
   to->source_direction = from->source_direction;
   to->next_stop = from->next_stop;
   to->dest_lp = from->dest_lp;
   to->count = from->count;
   to->packet_ID = from->packet_ID;
   to->chunk_id = from->chunk_id;
   to->sender_lp = from->sender_lp;
      
   for (i=0; i < N_dims; i++)
     to->dest[i] = from->dest[i];

   to->my_N_hop = from->my_N_hop;
   to->wait_type = from->wait_type;
   to->type = from->type;
}

void 
schedule_waiting_msg( nodes_state * s, 
			   tw_bf * bf, 
			   nodes_message * msg, 
			   tw_lp * lp )
{
  bf->c3 = 1;

  if( s->wait_count <= 0 || s->head == NULL)
   {
     return;
   }

  tw_event * e_h;
  nodes_message * m;
  tw_stime ts;

  //int i;
  //for(i = 0; i < s->wait_count; i++)
  waiting_packet * current = s->head;
  int loc=s->wait_count, j=0;

  for(j = 0; j < loc; j++)
   {
    if( s->waiting_list[j].dim == msg->source_dim && s->waiting_list[j].dir == msg->source_direction)
     {
        bf->c3=0;
	ts = tw_rand_exponential(lp->rng, MEAN_INTERVAL/10000);
        e_h = tw_event_new( lp->gid, ts, lp );
        m = tw_event_data( e_h );
     	copyData(s->waiting_list[j].packet, m);
        
//        For reverse computation, also copy data to the msg 
   	copyData(s->waiting_list[j].packet, msg);
        msg->wait_loc = j;

	m->type = SEND;
//	printf("\n Source dim %d Source dir %d loc %d wait type %d", current->dim, current->dir, loc, current->packet->wait_type);
//        printf("\n Packet ID %d ", current->packet->packet_ID);
        tw_event_send(e_h);
        waiting_packet_free(s, j);
        return;
    }
  }
}
/*Processes the packet after it arrives on the from the neighboring torus node */
void packet_arrive( nodes_state * s, 
		    tw_bf * bf, 
		    nodes_message * msg, 
		    tw_lp * lp )
{
  int i, delay = head_delay;//HOP_DELAY + tw_rand_exponential(lp->rng, HOP_DELAY/1000);
  tw_event *e;
  tw_stime ts;
  nodes_message *m;

//  s->N_wait_to_be_processed--;
  
  msg->my_N_hop++;
  bf->c1 = 1;
  ts = tw_rand_exponential(lp->rng, MEAN_INTERVAL/10000);
  if( lp->gid == msg->dest_lp )
    {   
	bf->c1 = 0;
        credit_send( s, bf, lp, msg); 
    
        if( msg->chunk_id == num_chunks - 1 )    
        {
		//ts = MEAN_PROCESS + OVERHEADS;
		//ts += MEAN_PROCESS;
	        e = tw_event_new(lp->gid + N_nodes, ts, lp);
		m = tw_event_data(e);
	        m->type = MPI_RECV;
	        m->travel_start_time = msg->travel_start_time;
//		m->origin_lp = msg->origin_lp;
		m->my_N_hop = msg->my_N_hop;
		m->packet_ID = msg->packet_ID;
		m->count = msg->count;
		tw_event_send(e);
	//	printf("\n Last message arrived at %lf travel start time %lf ", tw_now( lp ), m->travel_start_time );
        }
    }
  else
    {
//      Additional hop delay for large messages
      //if( mpi_message_size > 1024 )
	//delay += (HOP_DELAY/12 * num_packets);
      e = tw_event_new(lp->gid, ts , lp);
//      e = tw_event_new(lp->gid, HOP_DELAY * (msg->count + 1), lp);
      m = tw_event_data( e );
      m->type = SEND;
 //     m->origin_lp = msg->origin_lp;
      
      // Carry on the message info
      for( i = 0; i < N_dims; i++ )
	m->dest[i] = msg->dest[i];

      m->dest_lp = msg->dest_lp;
      
      m->source_dim = msg->source_dim;
      m->source_direction = msg->source_direction;
      
      m->packet_ID = msg->packet_ID;	  
      m->travel_start_time = msg->travel_start_time;

      m->my_N_hop = msg->my_N_hop;
      m->count = msg->count;
      m->sender_lp = msg->sender_lp;
      m->chunk_id = msg->chunk_id;

      m->next_stop = -1;

      tw_event_send(e);
   }
}

/*Each MPI LP in this model generates a MPI message until a certain message count is reached.
This method 
          (i) keeps generating MPI messages, 
	 (ii) breaks a MPI message down to torus packets 
         (iii) sends those packets to the underlying torus node LP */
void mpi_msg_send(mpi_process * p, 
	          tw_bf * bf, 
		  nodes_message * msg, 
		  tw_lp * lp)
{
    tw_stime ts;
    tw_event *e;
    nodes_message *m;
    tw_lpid final_dst;
    int i;//, pack_size = mpi_message_size;

    bf->c4 = 1;

    if(p->message_counter >= injection_limit)
     {
	bf->c4 = 0;

	return;
     }

    bf->c3 = 1;
    bf->c2 = 1;

    //if( mpi_message_size < 32 )
    //   pack_size = 32;

    //if(num_packets > 1)
    //   pack_size = PACKET_SIZE;

    // Create a packet
   // If ping-pong is turned on then only MPI proc 0 sends the message
   // If ping-pong is turned off then all MPI processes participate in the test (bisection test)
    //if( p->message_counter <= num_mpi_msgs )
    // {
     // final_dst = lp->gid  - distance ;

    switch(TRAFFIC)
	{
	  case UNIFORM_RANDOM:
		{
                    bf->c3 = 0;

		    final_dst = tw_rand_integer( lp->rng, N_nodes, 2 * N_nodes - 1);
	
		    if( final_dst == lp->gid )
		      {
                        final_dst = N_nodes + ((lp->gid + 1) % N_nodes);
		      }
		}
	  break;

	 case DRAGONFLY_ZONES:
		{
		    bf->c3 = 0;

  		    final_dst = N_nodes + ((p->zone_id + 1) % num_zones * NUM_ZONE_NODES) + tw_rand_integer( lp->rng, 0, NUM_ZONE_NODES-1);
		
//		     Commented out on May 17, destination is in next zone so there is not a possibility of lp id being the same as destination
		    /*if( final_dst == lp->gid )
                      {
                        bf->c1 = 0;

                        final_dst = N_nodes + ((lp->gid + 1) % N_nodes);
                      }*/
		}
	 break;

	case TRANSPOSE:
		{
		   if( p->col == p->row )
		    {
		      bf->c2 = 0;

		      final_dst = -1;
                      return;
                    }

		   final_dst = N_nodes + p->col * NUM_ROWS + p->row;
		}
         break;

       case NEAREST_NEIGHBOR:
         {
           final_dst = -1;
        }
        break;
	
	}
      //if(final_dst < N_nodes)
      //   final_dst += N_nodes;

      //printf("\n LP %d final dest %d ", (int)lp->gid, (int)final_dst);
#if DEBUG
//if(lp->gid == TRACK_LP + N_nodes)
//  printf("\n MPI Rank %d sending message to rank %d ", getProcID(lp->gid), getProcID(final_dst));
#endif
      tw_stime base_time = MEAN_PROCESS;
	
      for( i=0; i < num_packets; i++ ) 
       {
	      // Send the packet out
	     ts = 1.0 + tw_rand_exponential(lp->rng, MEAN_INTERVAL/100000); 
             msg->saved_available_time = p->available_time;
	     p->available_time = max( p->available_time, tw_now(lp) );
	     p->available_time += ts;

	     e = tw_event_new( getProcID(lp->gid), p->available_time - tw_now(lp), lp );

	     //e = tw_event_new( getProcID(lp->gid), ts, lp );
	     m = tw_event_data( e );
	     m->type = GENERATE;
             m->packet_ID = packet_offset * ( lp->gid * num_mpi_msgs * num_packets ) + p->message_counter;

             p->message_counter++;
	     m->count = i;
	     m->travel_start_time = tw_now( lp ) + ts;

	     if(final_dst == -1)
		m->dest_lp = -1;
	     else
	        m->dest_lp = getProcID( final_dst );

 	     m->next_stop = -1; 
//	     m->origin_lp = lp->gid;
	     //available_time += pack_size;
            //if( lp->gid == N_nodes)
             // printf("\n Sending message to %d ", (int)m->dest_lp); 
	    //int dst_proc = tw_rand_integer(lp->rng, 0, N_nodes-1);
	    //int dst_proc = s->neighbour_minus_lpID[0];
	    //m->dest_lp = dst_proc;
	
            tw_event_send( e );
     } 
     ts = 0.1 + tw_rand_exponential( lp->rng, MEAN_INTERVAL );
     e = tw_event_new( lp->gid, ts, lp );
     m = tw_event_data( e );
     m->type = MPI_SEND;
     tw_event_send( e );
  //}
}
/*Computes final latencies after a message is received */
void mpi_msg_recv(mpi_process * p, 
		  tw_bf * bf, 
		  nodes_message * msg, 
		  tw_lp * lp)
{
 // Message arrives at final destination
//  if( msg->count == num_packets -1 )
//  {
   bf->c3 = 1; 
   N_finished_msgs++;
    
// For torus-dragonfly comparison only place the end time here
    N_finished_packets++;
    int index = floor(N_COLLECT_POINTS*(tw_now(lp)/g_tw_ts_end));
    N_finished_storage[index]++;

    total_time += tw_now( lp ) - msg->travel_start_time;
    total_hops += msg->my_N_hop;

    if (max_latency < tw_now( lp ) - msg->travel_start_time) {
	  bf->c3 = 0;
	  msg->saved_available_time = max_latency;
          max_latency=tw_now( lp ) - msg->travel_start_time;
     }

 /*  if(msg->origin_lp == N_nodes)
     {
     // printf( "\n Origin lp time %lf packet ID %lld ", tw_now( lp ) - msg->travel_start_time, msg->packet_ID );
       total_lp_time += tw_now( lp ) - msg->travel_start_time;
       lp_hops = msg->my_N_hop;
     }
//  }*/
}
void mpi_event_handler( mpi_process * p, 
		       tw_bf * bf, 
		       nodes_message * msg, 
		       tw_lp * lp )
{
  switch(msg->type)
  {
     case MPI_SEND:
   	      mpi_msg_send(p, bf, msg, lp);
     break;

     case MPI_RECV:
              mpi_msg_recv(p, bf, msg, lp);
      break; 
  }
}

void mpi_event_rc_handler( mpi_process * p,
                           tw_bf * bf,
                           nodes_message * msg,
                           tw_lp * lp)
{
   switch(msg->type)
  {
     case MPI_SEND:
                {
		  if(!bf->c2 || !bf->c4)
		    return;

		  if(!bf->c3)
		       tw_rand_reverse_unif(lp->rng);

//		   Commented out on May 17, add 1 in the lp id instead of re-running random number generator
		  //if(!bf->c1)
		    //   tw_rand_reverse_unif(lp->rng);

		    int i; 
		    for( i=0; i < num_packets; i++ )
		     {
			     tw_rand_reverse_unif(lp->rng);

			     p->available_time = msg->saved_available_time;
			     p->message_counter--;			
		     }
		    tw_rand_reverse_unif(lp->rng);
	     }
     break;

     case MPI_RECV:
		  {
			int index = floor( N_COLLECT_POINTS * ( tw_now( lp ) / g_tw_ts_end ) );
			N_finished_msgs--;
			N_finished_packets--;
			N_finished_storage[index]--;
			total_time -= tw_now( lp ) - msg->travel_start_time;
			total_hops-=msg->my_N_hop;

			if(!bf->c3)
			  max_latency= msg->saved_available_time;
		  }
      break;
  }
}

void
final( nodes_state * s, tw_lp * lp )
{ 

}

tw_lp * torus_mapping_to_lp( tw_lpid lpid )
{
    int index;

    if(lpid < N_nodes)
       index = lpid - g_tw_mynode * nlp_nodes_per_pe;
    else
       index = nlp_nodes_per_pe + (lpid - g_tw_mynode * nlp_mpi_procs_per_pe - N_nodes);

    return g_tw_lp[index];
}

void packet_buffer_process( nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp )
{
#if DEBUG
//if( lp->gid == TRACK_LP )
//    printf( "\n Buffer space updated %d dir %d dim %d ", s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][  0 ], msg->source_direction, msg->source_dim );
#endif
//  printf("\n Credit arrived at %lf", tw_now(lp));
 // if( s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][ 0 ] >= 0 && s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][ 0 ] <= NUM_BUF_SLOTS)
//   printf("\n Incorrect buffer space %d ", s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][ 0 ]);
 // {
    //bf->c3 = 1;

    s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][  0 ]--;
  
    //msg->saved_direction = msg->source_direction;
    //msg->saved_source_dim = msg->source_dim; 

    //if( s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][ 0 ] >= 2 * PACKET_SIZE / TOKEN_SIZE )
    // {
       //bf->c3 = 0;

       schedule_waiting_msg( s, bf, msg, lp );
    // }
//  }
  //else
  //    printf("\n Buffer size exceeded normal range Should be: %d is %d ", NUM_BUF_SLOTS, s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][ 0 ]);

  // Send the message waiting in the queu+1e
//   {
//  if( lp->gid == TRACK_LP )
//    printf( "\n SCHEDULING PACKET buffer space %d %d ", s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][ 0 ], msg->packet_size );
//   }
}
void 
node_rc_handler(nodes_state * s, tw_bf * bf, nodes_message * msg, tw_lp * lp)
{
  switch(msg->type)
    {
       case GENERATE:
		   {
		     if(!bf->c2)
			return;

		     tw_rand_reverse_unif(lp->rng);	

		   }
	break;
	
	case ARRIVAL:
		   {
		     msg->my_N_hop--;
		     if(!bf->c1)
  		       tw_rand_reverse_unif(lp->rng);
		     tw_rand_reverse_unif(lp->rng);
		   }
	break;	

	case SEND:
		 {
	            tw_rand_reverse_unif(lp->rng);
		    if(!bf->c3)
		     {
			return;
		     }
		     else
		      {
			if(!bf->c2)
   			   tw_rand_reverse_unif(lp->rng);
                        int next_dim = msg->saved_src_dim;
			int next_dir = msg->saved_src_dir;

			s->next_link_available_time[next_dir + ( next_dim * 2 )][0] = msg->saved_available_time;
			
		        //printf("\n Sent reverse handler called %d %d %d ", msg->saved_source_dim, msg->saved_direction, s->buffer[ next_dir + ( next_dim * 2 ) ][ 0 ]);
			s->buffer[ next_dir + ( next_dim * 2 ) ][ 0 ] --;
			
                        if(!bf->c1)
			  {
			    int index = floor( N_COLLECT_POINTS * ( tw_now( lp ) / g_tw_ts_end ) );
			    N_generated_storage[ index ]--;
			  }
		      }
		 }
	break;

        case WAIT:
		{
		   //printf("\n Wait reverse handler called ");
//		    Gave this error message with roll backs, probably event ties causing this error message?
		   if(s->wait_count <= 0)
		    	printf("\n Trying to decrement a zero valued variable %d ", s->wait_count);

		   s->wait_count--;
		   int loc = s->wait_count;
		   s->waiting_list[loc].dim = -1;
		   s->waiting_list[loc].dir = -1;
		   s->waiting_list[loc].packet = NULL;
		}
        break;

       case CREDIT:
		{
		  s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][  0 ]++;
		  if(!bf->c3)
		  {
		  //printf("\n credit reverse handler called %d %d %d %d ", msg->source_direction, msg->source_dim, s->buffer[ msg->source_direction + ( msg->source_dim * 2 ) ][  0 ], s->wait_count);
		     tw_rand_reverse_unif(lp->rng);
		     int loc = msg->wait_loc, i;
                     int max_count = s->wait_count;
		     if(s->wait_count == WAITING_PACK_COUNT)
			printf("\n Exceeded maximum count!!! ");
		     for(i = max_count; i > loc ; i--)  
                      {
		  	  s->waiting_list[i].dim = s->waiting_list[i-1].dim;
			  s->waiting_list[i].dir = s->waiting_list[i-1].dir;
			  s->waiting_list[i].packet = s->waiting_list[i-1].packet;
                      }
		     s->waiting_list[loc].dim = msg->source_dim;
		     s->waiting_list[loc].dir = msg->source_direction;
		     s->waiting_list[loc].packet = msg;		 
		     s->wait_count++;
		 }
              }
       break;
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
  case CREDIT:
    packet_buffer_process(s,bf,msg,lp);
   break;
  case WAIT:
    update_waiting_list(s, msg, lp);
   break;
 }
}
tw_lptype nodes_lps[] =
{
	{
		(init_f) torus_init,
		(event_f) event_handler,
		(revent_f) node_rc_handler,
		(final_f) final,
		(map_f) mapping,
		sizeof(nodes_state),
	},
	{
               (init_f) mpi_init,
	       (event_f) mpi_event_handler,
	       (revent_f) mpi_event_rc_handler,
	       (final_f) final,
	       (map_f) mapping,
	       sizeof(mpi_process),
	},
	{0},
};

void torus_mapping(void)
{
  tw_lpid kpid;
  tw_pe * pe;
  int nkp_per_pe=16;

  for(kpid = 0; kpid < nkp_per_pe; kpid++)
   tw_kp_onpe(kpid, g_tw_pe[0]);

  int i;
  for(i = 0; i < nlp_nodes_per_pe; i++)
   {
     kpid = i % g_tw_nkp;
     pe = tw_getpe(kpid % g_tw_npe);
     tw_lp_onpe(i, pe, g_tw_mynode * nlp_nodes_per_pe + i + getRem() );
     tw_lp_onkp(g_tw_lp[i], g_tw_kp[kpid]);
     tw_lp_settype(i, &nodes_lps[0]);
   }
  for(i = 0; i < nlp_mpi_procs_per_pe; i++)
   {
     kpid = i % g_tw_nkp;
     pe = tw_getpe(kpid % g_tw_npe);
     tw_lp_onpe(nlp_nodes_per_pe+i, pe, N_nodes + g_tw_mynode * nlp_mpi_procs_per_pe + i + getRem() );
     tw_lp_onkp(g_tw_lp[nlp_nodes_per_pe + i], g_tw_kp[kpid]);
     tw_lp_settype(nlp_nodes_per_pe + i, &nodes_lps[1]);
   }
}
const tw_optdef app_opt [] =
{
	TWOPT_GROUP("Nodes Model"),
	TWOPT_UINT("memory", opt_mem, "optimistic memory"),
	TWOPT_ULONG("mpi-message-size", mpi_message_size, "mpi-message-size"),
//	TWOPT_UINT("num-mpi-msgs", num_mpi_msgs, "num-mpi-msgs"),
//	TWOPT_UINT("distance", distance, "distance"),
	TWOPT_UINT("mem_factor", mem_factor, "mem_factor"),
	TWOPT_UINT("traffic", TRAFFIC, "UNIFORM RANDOM=1, DRAGONFLY ZONES=2, TRANSPOSE=3, NEAREST NEIGHBOR=4"), 
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
        {
	  N_nodes*=dim_length[i];
	  N_mpi_procs*=dim_length[i];
	}
	nlp_nodes_per_pe = N_nodes/tw_nnodes()/g_tw_npe;
	nlp_mpi_procs_per_pe = N_mpi_procs/tw_nnodes()/g_tw_npe;

	total_lps = g_tw_nlp * tw_nnodes();
	node_rem = N_nodes % (tw_nnodes()/g_tw_npe);

	 /*if(g_tw_mynode < node_rem)
	   {
		nlp_nodes_per_pe++;
		nlp_mpi_procs_per_pe++;
	   }*/
        num_packets=1;
        num_chunks = PACKET_SIZE/CHUNK_SIZE;
        num_zones = N_nodes/NUM_ZONE_NODES;

	if( N_nodes % NUM_ZONE_NODES != 0)
	   num_zones++;
      
        if( mpi_message_size > PACKET_SIZE)
         {
          num_packets = mpi_message_size / PACKET_SIZE;  
          
	  if(mpi_message_size % PACKET_SIZE != 0 )
	    num_packets++;
         }

	g_tw_mapping=CUSTOM;
     	g_tw_custom_initial_mapping=&torus_mapping;
        g_tw_custom_lp_global_to_local_map=&torus_mapping_to_lp;

	g_tw_events_per_pe = mem_factor * 1024 * (nlp_nodes_per_pe/g_tw_npe + nlp_mpi_procs_per_pe/g_tw_npe) + opt_mem;
	tw_define_lps(nlp_nodes_per_pe + nlp_mpi_procs_per_pe, sizeof(nodes_message), 0);

	head_delay = (1 / BANDWIDTH) * CHUNK_SIZE;
	
        // BG/L torus network paper: Tokens are 32 byte chunks that is why the credit delay is adjusted according to bandwidth * 32
	credit_delay = (1 / BANDWIDTH) * 8;
	packet_offset = (g_tw_ts_end/MEAN_INTERVAL) * num_packets; 
	
	injection_limit = 20000/MEAN_INTERVAL;

	//printf("\n nlp_nodes_per_pe %d g_tw_nlp %d ", nlp_nodes_per_pe, (int)g_tw_nlp);

	if(tw_ismaster())
	{
		printf("\nTorus Network Model Statistics:\n");
		printf("\t%-50s %11d\n", "Number of nodes", N_nodes);
		printf(" Injection limit %d \n", injection_limit);
	}

	tw_run();
	unsigned long long total_finished_storage[N_COLLECT_POINTS];
	unsigned long long total_generated_storage[N_COLLECT_POINTS];
	unsigned long long wait_length,event_length,N_total_packets_finish, N_total_msgs_finish, N_total_hop;
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
	MPI_Reduce( &N_finished_packets, &N_total_packets_finish,1, 
		    MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce( &N_finished_msgs, &N_total_msgs_finish,1, 
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

//	if( total_lp_time > 0 )
//	   total_lp_time /= num_mpi_msgs;

	if(tw_ismaster())
	  {
	    printf("\n ****************** \n");
	    printf("\n total packets finished:         %lld and %lld; \n",
		   total_finished_storage[N_COLLECT_POINTS-1],N_total_packets_finish);
	    printf("\n total MPI messages finished:         %lld; \n",
                   N_total_msgs_finish);
	    printf("\n total generate:       %lld; \n",
		   total_generated_storage[N_COLLECT_POINTS-1]);
	    printf("\n total hops:           %lf; \n",
		   (double)N_total_hop/total_finished_storage[N_COLLECT_POINTS-1]);
	    printf("\n average travel time:  %lf; \n\n",
		   total_time_sum/N_total_msgs_finish);
	    
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

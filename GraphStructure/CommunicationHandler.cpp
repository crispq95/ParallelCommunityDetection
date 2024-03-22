#include "CommunicationHandler.h"
#include "DistributedGraph.h"


/*
 *    Class: CommunicationHandler  
 * Function: CommunicationHandler
 * --------------------
 * CommunicationHandler class constructor
 * 
 * -:-
 * 
 * returns: -
 */
CommunicationHandler::CommunicationHandler(){}

/*
 *    Class: CommunicationHandler  
 * Function: ~CommunicationHandler
 * --------------------
 * CommunicationHandler class destructor
 * 
 * -:-
 * 
 * returns: -
 */
CommunicationHandler::~CommunicationHandler(){}

/*
 *    Class: CommunicationHandler  
 * Function: send_data
 * --------------------
 * Sends the data from the send_buffer to neighbor PEs 
 * 
 * ghost_vertices: pointer to the vector of ghost_vertices of this PE
 * 
 * returns: -
 */
void CommunicationHandler::init_communications(std::vector<GhostNode> *ghost_vertices){
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    ID_T g_idx = 0; 
    no_of_neighbor_PEs = 0; 

    // get the neighbor PEs
    while ( g_idx < (*ghost_vertices).size() && no_of_neighbor_PEs < world_size){
        if(neighborPEs.find((*ghost_vertices)[g_idx].pe_id) == neighborPEs.end()){
            neighborPEs[(*ghost_vertices)[g_idx].pe_id] = no_of_neighbor_PEs;
            no_of_neighbor_PEs++; 
        }
        g_idx++; 
    }
    s_buffer.resize(no_of_neighbor_PEs); 
    send_request.resize(no_of_neighbor_PEs*2); // one for the msg size another for the msg itself 
    recv_request.resize(no_of_neighbor_PEs*2); 
    rcv_buffer.resize(no_of_neighbor_PEs); 

    // DPC _ LPA 
    rcv_buffer_weight.resize(no_of_neighbor_PEs); 
    s_buffer_weight.resize(no_of_neighbor_PEs); 
    
}

/*
 *    Class: CommunicationHandler  
 * Function: clear_buffers
 * --------------------
 * Clears the buffers of the CommunicationHandler 
 * 
 * -:-
 * 
 * returns: -
 */
void CommunicationHandler::clear_buffers()
{
    for( int i = 0 ; i < s_buffer.size() ; i++ )
        s_buffer[i].clear(); 
    for( int i = 0; i< rcv_buffer.size(); i++)
        rcv_buffer[i].clear();
}

/*
 *    Class: CommunicationHandler  
 * Function: add_all_to_send
 * --------------------
 * Adds a local vertex to the send_buffer of neighborPEs
 * 
 * pe_ids: global_IDs of the neighbor PEs of the local vertex
 * global_id: global_ID of the local vertex to be sent 
 * label: label to be sent 
 * 
 * returns: -
 */
void CommunicationHandler::add_all_to_send(std::unordered_set<int> * pe_ids, ID_T global_id, ID_T label)
{
    for ( auto pe_id : (*pe_ids) ){
        int pe_idx = neighborPEs[pe_id];
        s_buffer[pe_idx].push_back(global_id);
        s_buffer[pe_idx].push_back(label);  
    }
}

/*
 *    Class: CommunicationHandler  
 * Function: add_label_to_send
 * --------------------
 * Adds the LABEL of the local vertex to the send_buffer of neighborPEs
 * To use this function we must first order the vtx on each PE to work properly on the send/rcv 
 * 
 * pe_ids: global_IDs of the neighbor PEs of the local vertex
 * label: label to be sent 
 * 
 * returns: -
 */
void CommunicationHandler::add_label_to_send(std::unordered_set<int> * pe_ids, ID_T label){
    for ( auto pe_id : (*pe_ids) ) s_buffer[neighborPEs[pe_id]].push_back(label);  
}


/*
 *    Class: CommunicationHandler  
 * Function: send_data
 * --------------------
 * Sends the data from the send_buffer to neighbor PEs 
 * 
 * -:-
 * 
 * returns: -
 */
void CommunicationHandler::send_data(){
    MPI_Request request;
    // send msg sizes 
    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        unsigned long int msg_size = s_buffer[pe_idx].size();

        MPI_Isend(&s_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, &request); 
    }
}

/*
 *    Class: CommunicationHandler  
 * Function: recv_data
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void CommunicationHandler::recv_data(){
   for ( auto n_PE : neighborPEs ){
        MPI_Status status; 
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        int msg_size; 

        MPI_Probe(pe_rank, 11, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &msg_size);

        // std::vector<ID_T> temp_rcv_data;
        // temp_rcv_data.resize(msg_size); 
        rcv_buffer[pe_idx].resize(msg_size); 

        // if(msg_size > 0)
        MPI_Recv(&rcv_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
   }
}

/*
 *    Class: CommunicationHandler  
 * Function: wait_requests
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void CommunicationHandler::wait_requests(){
    clear_buffers();
}

/*
 *    Class: CommunicationHandler  
 * Function: wait_requests
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void CommunicationHandler::send_rcv_inactive(DistributedGraph* graph){
    int no_local_vtx = graph->get_local_vertices()->size();

    int num_deg_1 = 0; 

    // Send and recieve data from the boundary inactive vtx 
    for ( auto vtx : (*graph->get_local_vertices()) ) {
        if( !vtx.is_boundary || vtx.active )
            continue; 
        
        std::unordered_set<int> boundary_neighbor_PEs;

        for( auto neigh : (*vtx.edges) ){
            if( neigh.target < no_local_vtx )
                continue; 

            int local_idx =  graph->from_local_ghost_to_index(neigh.target);
            int ghost_pe_id = graph->get_ghost_vertices()->at(local_idx).pe_id;

            // if the neighbor is a ghost -> keep track of the PE it belongs to 
            if(boundary_neighbor_PEs.find(ghost_pe_id) == boundary_neighbor_PEs.end()){
                boundary_neighbor_PEs.insert(ghost_pe_id); 
                if(boundary_neighbor_PEs.size() >= no_of_neighbor_PEs)
                    break; 
            }
        }

        ID_T vtx_global_id = graph->from_local_to_global(vtx.id);
        for ( auto i : boundary_neighbor_PEs ) {
            int pe_idx = neighborPEs[i];
            s_buffer[pe_idx].push_back(vtx_global_id);
            num_deg_1++; 
        }
    }

    // std::cout << "On rank " << my_rank << " I found " << num_deg_1 << " boundary vtx of degree 1 " << std::endl; 

    // send the data to neighbor PEs 
    MPI_Request request;
    // send msg sizes 
    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        int msg_size = s_buffer[pe_idx].size();

        MPI_Isend(&s_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, &request); 
    }

    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        int msg_size;
        MPI_Status status; 

        MPI_Probe(pe_rank, 21, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &msg_size);

        rcv_buffer[pe_idx].resize(msg_size); 

        // if(msg_size > 0)
        MPI_Recv(&rcv_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
        // if(msg_size > 0){
        //     std::cout << "On rank " << my_rank << " I GOT " << rcv_buffer[pe_idx].size() << " boundary vtx of degree 1 from " << pe_rank << " : "; 
        //     for(auto r : rcv_buffer[pe_idx])
        //         std::cout << r << " " << std::endl; 
        //     std::cout << std::endl; 
        // }
    }

}

// test 
void CommunicationHandler::send_recv_data(){
    MPI_Request request;

    /* Sending data to neighbors */
    // For all neighbor PEs 
    // for ( auto n_PE : neighborPEs ){
    //     // n_PE : first = RANK / second = id 
    //     int pe_rank = n_PE.first, pe_idx = n_PE.second;
        
    //     // this message sends the amount of data to be sent to this process 
    //     unsigned long int msg_size = s_buffer[pe_idx].size();
    //     MPI_Isend(&msg_size, 1, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, &request); // tag of nº = 11 

    //     if(msg_size != 0){
    //         MPI_Isend(&s_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, &request); // tag of nº = 21 
    //     }
    // }

    // send msg sizes 
    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        unsigned long int msg_size = s_buffer[pe_idx].size();

        MPI_Isend(&s_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, &request); 
    }


    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        int msg_size;
        MPI_Status status; 

        MPI_Probe(pe_rank, 11, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &msg_size);

        rcv_buffer[pe_idx].resize(msg_size); 

        if(msg_size > 0)
            MPI_Recv(&rcv_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // for ( auto n_PE : neighborPEs ){
    //     int pe_rank = n_PE.first, pe_idx = n_PE.second;
    //     unsigned long int msg_size; 
        
    //     // get how many ghosts I'm recvng 
    //     MPI_Recv(&msg_size, 1, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //     // rcv the data from process pe_rank 
    //     if(msg_size != 0){
    //         MPI_Recv(&rcv_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         // for( int a=0 ; a < temp_rcv_data.size() ; a+=2 ){
    //         //     // T local_g_indx = g->ghost_global_ids[temp_rcv_data[a]];
    //         //     T local_g_indx = g->from_ghost_global_to_index(temp_rcv_data[a]);
    //         //     LABEL_T g_label = temp_rcv_data[a+1];

    //         //     // update label 
    //         //     // (*g->ghost_vertices)[local_g_indx].current_label = g_label;
    //         //     (*g->get_ghost_vertices())[local_g_indx].current_label = g_label;
    //         // }
    //     }
    // }
    // clear_buffers();
} 


/*
 *    Class: CommunicationHandler  
 * Function: order_ghosts
 * --------------------
 * Creates an ordered index array to decide the order of recieved ghosts (if sending only labels)
 * 
 * ghost_vertices : Pointer to the vector of ghost vertices of this process
 * 
 * returns: -
 */
void CommunicationHandler::order_ghosts(std::vector<GhostNode> *ghost_vertices){
    sort_indexes_by_label<GhostNode>((*ghost_vertices), ordered_ghost_indices);
} 

/* Code for DPC LPA */
/* Nothing here is tested : */
// adds seeds to be sent, its label and weight included for each seed candidate
void CommunicationHandler::add_seeds_to_send(const std::unordered_set<int>& pe_ids, ID_T global_id, double weight){
    for ( auto pe_id : pe_ids ){
        int pe_idx = neighborPEs[pe_id];

        s_buffer[pe_idx].push_back(global_id);
        s_buffer_weight[pe_idx].push_back(weight);  
    }
}

// sends and recieves info of seed candidates 
void CommunicationHandler::send_recv_candidate_data(){
    MPI_Request request;

    // send msg sizes 
    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        unsigned long int msg_size = s_buffer[pe_idx].size();

        MPI_Isend(&s_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, &request); 
        MPI_Isend(&s_buffer_weight[pe_idx][0], msg_size, MPI_DOUBLE, pe_rank, 21, MPI_COMM_WORLD, &request); 
    }


    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        int msg_size;
        MPI_Status status; 

        MPI_Probe(pe_rank, 11, MPI_COMM_WORLD, &status);
        MPI_Probe(pe_rank, 21, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_UNSIGNED_LONG, &msg_size);
        
        rcv_buffer[pe_idx].resize(msg_size); 
        rcv_buffer_weight[pe_idx].resize(msg_size); 

        if(msg_size > 0){
            MPI_Recv(&rcv_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rcv_buffer_weight[pe_idx][0], msg_size, MPI_DOUBLE, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
} 

// addes candidate that changed into queue -> must be changed later to send only <label_id:list_ids, label_id2:list_ids, etc.> 
void CommunicationHandler::add_candidate_to_send(const std::unordered_set<int>& pe_ids, ID_T global_id, ID_T label, double weight){
    for ( auto pe_id : pe_ids ){
        int pe_idx = neighborPEs[pe_id];

        s_buffer[pe_idx].push_back(global_id);
        s_buffer[pe_idx].push_back(label);
        s_buffer_weight[pe_idx].push_back(weight);  
    }
}
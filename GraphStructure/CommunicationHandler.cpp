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
    for( int i = 0 ; i < s_buffer.size() ; i++ ){
        s_buffer[i].clear(); 
        
    }

    for( int i = 0; i< rcv_buffer.size(); i++){
        rcv_buffer[i].clear();
    }
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

/* TO DO */
// una funcion que guarde los ids en orden de los ghosts por PE 
// void CommunicationHandler::order_neighborPE_ghosts(std::vector<GhostNode> *gn){
//     for( auto g : (*gn)){

//     }
// }

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
    
        // if(my_rank == 0 || my_rank == 1)
        // {
        //     std::cout << "sending data from " << my_rank << " to " << pe_rank << " : ";
        //     for(int i = 0; i < s_buffer[pe_idx].size() ; i+=2){
        //         std::cout << "[" << s_buffer[pe_idx][i] << ", " << s_buffer[pe_idx][i+1] << "]" << " ";
        //     }
        //     std::cout << std::endl;
        // }
    }
}

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

        // if(my_rank == 2){
        //     std::cout << "Recv data from " << pe_rank << " : ";
        //     for(int i = 0; i < rcv_buffer[pe_idx].size() ; i+=2){
        //         std::cout << "[" << rcv_buffer[pe_idx][i] << ", " << rcv_buffer[pe_idx][i+1] << "]" << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // if(msg_size != 0){
        //     for( int a=0 ; a < rcv_buffer[pe_idx].size() ; a+=2 ){
        //         ID_T local_g_indx = g->from_ghost_global_to_index(rcv_buffer[pe_idx][a]);

        //         // update label 
        //         (*g->get_ghost_vertices())[local_g_indx].current_label = rcv_buffer[pe_idx][a+1];
        //     }
        // }
   }
}

/* WORK IN PROGRESS */
void CommunicationHandler::send_recv_data(DistributedGraph* g){
    MPI_Request request;
    
    // if(my_rank == 2){
    //     for( int i=0 ; i< s_buffer[1].size(); i+=2){
    //         std::cout << "Value of vtx sending to rank 1 id : " << s_buffer[1][i] << " val :" << s_buffer[1][i+1] << "to rank 1"   << std::endl;
    //     }
    // }

    /* Sending data to neighbors */
    // For all neighbor PEs 
    for ( auto n_PE : neighborPEs ){
        // n_PE : first = RANK / second = id 
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        
        // this message sends the amount of data to be sent to this process 
        unsigned long int msg_size = s_buffer[pe_idx].size();
        MPI_Isend(&msg_size, 1, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, &request); // tag of nº = 11 

        if(msg_size != 0){
            MPI_Isend(&s_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, &request); // tag of nº = 21 
        }
    }


    for ( auto n_PE : neighborPEs ){
        int pe_rank = n_PE.first, pe_idx = n_PE.second;
        unsigned long int msg_size; 
        
        // get how many ghosts I'm recvng 
        MPI_Recv(&msg_size, 1, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        std::vector<ID_T> temp_rcv_data;
        temp_rcv_data.resize(msg_size); 

        rcv_buffer[pe_idx].resize(msg_size); 

        // rcv the data from process pe_rank 
        if(msg_size != 0){
            // MPI_Recv(&temp_rcv_data[0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&rcv_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // if(msg_size > 0)
            // for( int a=0 ; a < temp_rcv_data.size() ; a+=2 ){
            //     // T local_g_indx = g->ghost_global_ids[temp_rcv_data[a]];
            //     ID_T local_g_indx = g->from_ghost_global_to_index(temp_rcv_data[a]);
            //     LABEL_T g_label = temp_rcv_data[a+1];

            //     if(my_rank == 1 && pe_rank == 2)
            //         std::cout << "Recv : " << temp_rcv_data[a] << " with label : " << g_label << std::endl; 
            //     // update label 
            //     // (*g->ghost_vertices)[local_g_indx].current_label = g_label;
            //     // (*g->get_ghost_vertices())[local_g_indx].current_label = g_label;
            // }
        }
    }
    // clear_buffers();
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


// template <typename T>
// std::vector<size_t> sort_indexes(const std::vector<T> &v) {
//   // initialize original index locations
//   std::vector<size_t> idx(v.size());
//   iota(idx.begin(), idx.end(), 0);

//   // sort indexes based on comparing values in v
//   // using std::stable_sort instead of std::sort
//   // to avoid unnecessary index re-orderings
//   // when v contains elements of equal values 
//   stable_sort(idx.begin(), idx.end(),
//        [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

//   return idx;
// }

/* TO DO */
void CommunicationHandler::order_ghosts(std::vector<GhostNode> *ghost_vertices){
    // std::vector<std::vector<ID_T>> ghost_ids_vector;
    // // std::vector<int> pe_g_cnt;
    // ghost_ids_vector.reserve(no_of_neighbor_PEs);
    // ordered_ghost_indices.reserve(no_of_neighbor_PEs);
    // // pe_g_cnt.reserve(no_of_neighbor_PEs); 
    // // std::fill(pe_g_cnt.begin(), pe_g_cnt.end(), 0); // init to 0s 

    // // We want, for all boundary vtx -> get the 
    // for( ID_T i = 0 ; i < (*ghost_vertices).size() ; i++ ){
    //     // get its global id and the PE it belongs to 
    //     ID_T global_id = (*ghost_vertices)[i].current_label;
    //     int loc_pe = neighborPEs[(*ghost_vertices)[i].pe_id];

    //     std::cout << "loc pe : " << loc_pe << std::endl; 
    //     ghost_ids_vector[0].push_back(global_id);
    // //     // ordered_ghost_indices[loc_pe].push_back(i);
    // }

    // //once we have all neighbor global_ids and indices
    // // for( int i = 0 ; i < no_of_neighbor_PEs ; i++ ){
    // //     ghost_ids_vector[i].reserve(ghost_ids_vector[i].size());
    // //     for( auto a : sort_indexes<ID_T>((ghost_ids_vector[i]))){
    // //         ordered_ghost_indices[i].push_back(a);  
    // //     }
    // // }

    // // if(my_rank != 0){
    // //     std::cout << "On rank " << my_rank << " we have the following indexes : ";
    // //     for ( auto a : ordered_ghost_indices[0] ){
    // //         std::cout << a << " ";
    // //     }
    // //     std::cout << std::endl; 
    // // }
} 
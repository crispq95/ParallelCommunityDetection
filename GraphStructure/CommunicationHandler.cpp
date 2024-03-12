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
void CommunicationHandler::send_rcv_inactive(std::vector<LocalNode> *local_vertices){
    // Send and recieve data from the boundary inactive vtx 
    // for ( auto vtx : (*local_vertices) ) {
    //     if( !vtx.is_boundary || vtx.active )
    //         continue; 
        
    //     std::unordered_map<int, int> boundary_neighbor_PEs;

    //     for( auto neigh : (*vtx.edges) ){
    //         if( neigh->id < no_of_vtx_PE)
    //             continue; 
    //         int ghost_pe_id = graph->get_ghost_vertices()->at(neigh->id).pe_id;

    //         // if the neighbor is a ghost -> keep track of the PE it belongs to 
    //         if(boundary_neighbor_PEs.find(ghost_pe_id) == boundary_neighbor_PEs.end())
    //             boundary_neighbor_PEs.insert(ghost_pe_id); 

    //         if(boundary_neighbor_PEs.size() == no_of_neighbor_PEs)
    //             break; 
    //     }

    //     for ( auto i : boundary_neighbor_PEs ) {
    //         int pe_idx = neighborPEs[i];
    //         ID_T vtx_global_id = graph->from_local_to_global(vtx.id);
    //         s_buffer[pe_idx].push_back(vtx_global_id);
    //     }
    // }

    // // send the data to neighbor PEs 
    // MPI_Request request;
    // // send msg sizes 
    // for ( auto n_PE : neighborPEs ){
    //     int pe_rank = n_PE.first, pe_idx = n_PE.second;
    //     unsigned long int msg_size = s_buffer[pe_idx].size();

    //     MPI_Isend(&s_buffer[pe_idx][0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, &request); 
    // }

    // for ( auto n_PE : neighborPEs ){
    //     int pe_rank = n_PE.first, pe_idx = n_PE.second;
    //     unsigned long int msg_size;
    //     MPI_Status status; 

    //     MPI_Probe(pe_rank, 11, MPI_COMM_WORLD, &status);
    //     MPI_Get_count(&status, MPI_UNSIGNED_LONG, &msg_size);

    //     rcv_buffer[pe_idx].resize(msg_size); 

    //     if(msg_size > 0)
    //         MPI_Recv(&rcv_buffer[pe_idx][0], );
    // }
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
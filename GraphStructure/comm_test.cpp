#include "comm.h"

CommunicationHandler::CommunicationHandler(){}
CommunicationHandler::~CommunicationHandler(){}

/* Init the required info for communications */
void CommunicationHandler::init_communications(std::vector<GhostNode> *ghost_nodes){
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    T g_idx = 0; 
    no_of_neighbor_PEs = 0; 

    // get the neighbor PEs
    while ( g_idx < (*ghost_nodes).size() && no_of_neighbor_PEs < world_size){
        if(neighborPEs.find((*ghost_nodes)[g_idx].pe_id) == neighborPEs.end()){
            neighborPEs[(*ghost_nodes)[g_idx].pe_id] = no_of_neighbor_PEs;
            no_of_neighbor_PEs++; 
        }
        g_idx++; 
    }
    s_buffer.resize(no_of_neighbor_PEs); 
}

void CommunicationHandler::clearBuffers()
{
    for( int i = 0 ; i < s_buffer.size() ; i++ )
        s_buffer[i].clear(); // rmv elements for each pe ? 
}

// pseudocode to attach data to each send buffer
void CommunicationHandler::addToSend(DistributedGraph *g, LocalNode node)
{
    if(node.is_boundary){

        T label = (*g->local_nodes)[node.id].next_label; //  current_label
        T node_id = (*g->local_nodes)[node.id].id;

        const std::vector<Edge>* neighbors = g->get_neighbors(node.id);
        for (auto n : (*neighbors)){
            if(g->is_ghost(n.target)){
                T indx_of_ghost = g->from_local_ghost_to_index(n.target);

                int PEid = (*g->ghost_nodes)[indx_of_ghost].pe_id; // not really like this tbh 
                int pe_idx = neighborPEs[PEid];

                T global_id = g->from_local_to_global(node_id); // change type 

                // if this node is set to be sent to this PE ? nothing : else send 
                // change this for array to check if something is up to be sent ? 
                bool already_in = false; 
                for ( T i = 0; i < s_buffer[pe_idx].size() ; i+=2 )
                {
                    if(global_id == s_buffer[pe_idx][i]){
                        already_in = true; 
                        break;
                    }
                }

                // needs to use global ids to send info : 
                if(!already_in){
                    s_buffer[pe_idx].push_back(global_id);
                    s_buffer[pe_idx].push_back(label);    
                }
            }
        }
    }
}

/* WORK IN PROGRESS */
void CommunicationHandler::send_recv_data(DistributedGraph* g){
    MPI_Request request;

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

        std::vector<T> temp_rcv_data;
        temp_rcv_data.resize(msg_size); 

        // rcv the data from process pe_rank 
        if(msg_size != 0){
            MPI_Recv(&temp_rcv_data[0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for( int a=0 ; a < temp_rcv_data.size() ; a+=2 ){
                T local_g_indx = g->ghost_global_ids[temp_rcv_data[a]];
                LABEL_T g_label = temp_rcv_data[a+1];

                // update label 
                (*g->ghost_nodes)[local_g_indx].current_label = g_label;
            }
        }
    }
    clearBuffers();
}
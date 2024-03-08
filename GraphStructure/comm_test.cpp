#include "comm.h"
// // Handle the required Information for communications 

CommunicationHandler::CommunicationHandler(){}
CommunicationHandler::~CommunicationHandler(){}

/* Init the required info for communications */
void CommunicationHandler::init_communications(std::vector<GhostNode> *ghost_nodes){
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    T g_idx = 0; 
    no_of_neighbor_PEs = 0; 

    // const int test = 1 ; 
    // if(std::find(neighborPEs.begin(), neighborPEs.end(), test) != neighborPEs.end()){
    //     std::cout << "Found : " <<  std::endl; 
    // }

    // std::cout << " I have " << (*ghost_nodes).size() << " ghosts" << std::endl; 
    // get the neighbor PEs
    while ( g_idx < (*ghost_nodes).size() && no_of_neighbor_PEs < world_size){
        if(neighborPEs.find((*ghost_nodes)[g_idx].pe_id) == neighborPEs.end()){
            neighborPEs[(*ghost_nodes)[g_idx].pe_id] = no_of_neighbor_PEs;
            // n_ranks.push_back((*ghost_nodes)[g_idx].pe_id); 
            no_of_neighbor_PEs++; 
            
            // if(my_rank == 1)
            //     std::cout << (*ghost_nodes)[g_idx].pe_id << " - " << neighborPEs[(*ghost_nodes)[g_idx].pe_id] << " " << std::endl; 
        }
       
        // if(std::find(neighborPEs.begin(), neighborPEs.end(), (*ghost_nodes)[g_idx].pe_id) == neighborPEs.end()){
        //     neighborPEs.push_back((*ghost_nodes)[g_idx].pe_id);
        //     no_of_neighbor_PEs++;

        //     // std::cout << "MIAU : " << neighborPEs[(*ghost_nodes)[g_idx].pe_id] << std::endl;
        // }
        // else{
        //     // std::cout << "HOLA : " << neighborPEs[(*ghost_nodes)[g_idx].pe_id] << std::endl;
        // }
        g_idx++; 
    }

    

    // if(my_rank == 0){
    //     std::cout << "Found " << no_of_neighbor_PEs << " neighbors." <<  std::endl; 
    //     for( auto a : neighborPEs){
    //         std::cout << a.first << " -> " << a.second << std::endl; 
    //     }
    // }

    s_buffer.resize(no_of_neighbor_PEs); 
    rcv_buffer.resize(no_of_neighbor_PEs); 

    // if(my_rank == 0){
    //     std::cout << "Found " << no_of_neighbor_PEs << " neighbor PEs." << std::endl; 

    //     for(auto i : neighborPEs){
    //         std::cout << i.first << " - " << i.second << " " << std::endl; 
    //     }
    // }


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

        // get the data itself 
        // rcv_buffer[i].resize(msg_size);

        std::vector<T> temp_rcv_data;
        temp_rcv_data.resize(msg_size); 

        // std::cout << "I'm on " << my_rank << " and rcv " << msg_size/2 << " vertices from " << pe_rank << std::endl; 
        // rcv the data from process pe_rank 
        if(msg_size != 0){
            // std::cout << "Aiming to recv [" << msg_size << "] on " << my_rank << " from " << pe_rank << std::endl; 
            // MPI_Recv(&temp_rcv_data[0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&temp_rcv_data[0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // std::cout << "Temp data on " << my_rank << " from : " << pe_idx << std::endl;  
            for( int a=0 ; a < temp_rcv_data.size() ; a+=2 ){
                // get local_id 
                T local_g_id = g->ghost_global_ids[temp_rcv_data[a]];

                // get local_index 
                T local_g_idx = g->from_local_ghost_to_index(local_g_id); 

                // get label 
                LABEL_T g_label = temp_rcv_data[a+1];

                // update label 
                // std::cout << "Global id : " << temp_rcv_data[a] << "[LOCAL: " << g->ghost_global_ids[temp_rcv_data[a]] << "] , label " << temp_rcv_data[a+1] << " changes from " << (*g->ghost_nodes)[local_g_idx].current_label << std::endl; 
                (*g->ghost_nodes)[local_g_id].current_label = g_label;
            }
        }
    }
    clearBuffers();
}
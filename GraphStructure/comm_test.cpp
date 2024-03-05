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
            n_ranks.push_back((*ghost_nodes)[g_idx].pe_id); 
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
    // for( int i = 0 ; i < (*s_buffer).size() ; i++ )
    //     s_buffer[i].resize(0); // rmv elements for each pe ? 
}

// pseudocode to attach data to each send buffer
void CommunicationHandler::addToSend(DistributedGraph *g, LocalNode node)
{
    // if (my_rank == 1 )
    //     std::cout << "I'm adding a vtx " << node.id ;
    // if the node is boundary :
    if(node.is_boundary){
        // if (my_rank == 1 )
        //     std::cout << " and it is BOUNDARY" << std::endl;
        nodePackage vtx_to_send;
        // vtx_to_send.label = node.next_label; //  current_label
        // vtx_to_send.node_id = node.id;

        vtx_to_send.label = (*g->local_nodes)[node.id].next_label; //  current_label
        vtx_to_send.node_id = (*g->local_nodes)[node.id].id;

        const std::vector<Edge>* neighbors = g->get_neighbors(node.id);
        // if (my_rank == 1 )
        //     std::cout << "It has : " << (*neighbors).size() << " neighbors ";

            // int cnt = 0;
    //     // if(my_rank==0)
    //     //     std::cout << "Im a "<< node.is_boundary << " b node " << node.id << " at rank " << my_rank << " and I have " << neighbors->size() << " neighbors." << std::endl;
    //     // for all its neighbors 
        for (auto n : (*neighbors)){
            // if(my_rank == 1)
            //     std::cout << n.target << " " << std::endl; 
            // if the neighbor is a ghost 
            if(g->is_ghost(n.target)){
                T indx_of_ghost = g->from_local_ghost_to_index(n.target);

                int PEid = (*g->ghost_nodes)[indx_of_ghost].pe_id; // not really like this tbh 
                int pe_idx = neighborPEs[PEid];

                // if(my_rank == 1)
                //     std::cout << n.target << " is a ghost with index " << indx_of_ghost << " found on PE " << PEid <<  std::endl;

                int global_id = g->from_local_to_global(vtx_to_send.node_id); // change type 

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

                if(!already_in){
                    // needs to use global ids to send info : 

                    // if(my_rank == 1)
                    //     std::cout << "From local " << vtx_to_send.node_id << " to global " << global_id << std::endl; 
                    s_buffer[pe_idx].push_back(global_id);
                    s_buffer[pe_idx].push_back(vtx_to_send.label);    
                }
            }
        }
    }

    // if(my_rank==0){
    // // // //     std::cout << "MIAU" << std::endl; 
    // //     std::cout << "Sending to process " << neighborPEs[0] << " labeled as " << 0 << "- " << s_buffer[0].size()/2 << "vtx" << std::endl; 
    // //     for ( int i = 0 ; i < s_buffer[0].size() ; i+=2){
    // //         std::cout << s_buffer[0][i] << ", " << s_buffer[0][i+1] << " | " << std::endl;
    // //     }

    //     for(int i = 0 ; i < no_of_neighbor_PEs ; i++ ){
    //         int n_PE = neighborPEs[i];
    //         std::cout << "Detected from " << my_rank << " to process " << n_PE << " vtx : " << s_buffer[n_PE].size() << std::endl; 
    //     }
    // }
}

/* WORK IN PROGRESS */
void CommunicationHandler::send_recv_data(DistributedGraph* g){
    MPI_Request request;

    // if(my_rank == 0){
    //     std::cout << "printing neighbors at rank " << my_rank << std::endl;
    //     for(auto a : neighborPEs){
    //         // if(n_PE == a.second){
    //             // pe_rank = a.first; 
    //             std::cout << "Rank " << a.first << " neighbor index : " << a.second << std::endl;
    //             // break; 
    //         // }
    //     }
    // }

    /* Sending data to neighbors */
    // For all neighbor PEs 
    for(int i = 0 ; i < no_of_neighbor_PEs ; i++ ){
        // get the rank of the neighbor PE & get the number of vtx to send (*2)
        int n_PE = neighborPEs[i], msg_size; 
        int pe_rank = n_ranks[i], pe_idx;

        for(auto a : neighborPEs){
            if(pe_rank == a.first){
                pe_idx = a.second; 
                break; 
            }
        }

        msg_size = s_buffer[pe_idx].size();

        // if(my_rank == 0)
            // std::cout << "Sending " << msg_size << " vtx from " << my_rank << " to process with rank " << pe_rank << std::endl; 
        // this message sends the amount of data to be sent to this process 
        // MPI_Send(&size, 1, MPI_UNSIGNED_LONG, n_PE, 11, MPI_COMM_WORLD);
        MPI_Isend(&msg_size, 1, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, &request); // tag of nº = 11 
        

        // send the data iself 
        std::vector<T> temp_data; // array needs to be flattened for send to work ? 
        for( int i = 0; i< s_buffer[pe_idx].size(); i++){
            temp_data.push_back(s_buffer[pe_idx][i]);
        }

        // if(my_rank == 1)
        if(temp_data.size() != 0){
            MPI_Isend(&temp_data[0], temp_data.size(), MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, &request); // tag of nº = 21 
            
            // std::cout << "sending from " << my_rank << " to " << pe_rank << " a buffer of " << msg_size << " stored at "  << pe_idx << "[ Confirmed by : " << neighborPEs[pe_rank] << "]" <<  std::endl;

        }
        // else{
        //     std::cout << "Not sending anything from " << my_rank << " to " << pe_rank << std::endl; 

        // }
            // MPI_Isend(&temp_data[0], 1, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, &request); // tag of nº = 21 
    }

    /* Reciving data from neighbors */
    for( int i = 0 ; i < no_of_neighbor_PEs ; i++ ){
        // get the rank of the neighbor PE 
        int n_PE = neighborPEs[i], msg_size; 
        int pe_rank = n_ranks[i], pe_idx;

        for(auto a : neighborPEs){
            if(pe_rank == a.first){
                pe_idx = a.second; 
                break; 
            }
        }


        // get how many ghosts I'm recvng 
        MPI_Recv(&msg_size, 1, MPI_UNSIGNED_LONG, pe_rank, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // get the data itself 
        rcv_buffer[i].resize(msg_size);
        // rcv_buffer[pe_idx].push_back(1);
        
        std::vector<T> temp_rcv_data;
        temp_rcv_data.resize(msg_size); 


        // MPI_Recv(temp_rcv_data.data(), 1, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(msg_size != 0){
            // std::cout << "Aiming to recv [" << msg_size << "] on " << my_rank << " from " << pe_rank << std::endl; 
            MPI_Recv(&temp_rcv_data[0], msg_size, MPI_UNSIGNED_LONG, pe_rank, 21, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // if(my_rank == 0){
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
            // }
        }

        // std::cout << "Gonna store on rcv_buffer[" << pe_idx << "] the " << msg_size << " vertices arriving on " << my_rank << " from " << pe_rank << std::endl;
        // if(my_rank != 0)
        // std::cout << "Reciving " << temp_rcv_data.size() << " vtx in " << my_rank << " from process with rank " << pe_rank << std::endl; 

    }

    // std::cout << my_rank << " -- IS DONE !!________ " << std::endl; 
    

    // determine the amount of data to be sent to this PE 
    // for( int i = 0 ; i < no_of_neighbor_PEs ; i++ ){
    //     int n_PE = neighborPEs[i], size = (*s_buffer)[n_PE].size(); 
    //     std::cout << "Buffer size : " << (*s_buffer)[n_PE].size()<< std::endl; 
        
    //     // send the size of the label msg 
    //     MPI_Isend(&size, 1, MPI_LONG_INT, n_PE, 11, MPI_COMM_WORLD, &request); // tag of nº = 11 

    //     // send the msg with {id:label} entry for each vtx 
    //     MPI_Isend(&((*s_buffer)[n_PE]), (*s_buffer)[n_PE].size(), MPI_LONG_INT, n_PE, 0, MPI_COMM_WORLD, &request); 
    // }
    
//     MPI_Request requests[no_of_neighbor_PEs];
//     MPI_Status statuses[no_of_neighbor_PEs];
    
//     for( int i = 0 ; i < no_of_neighbor_PEs ; i++ ){
//         int n_PE = neighborPEs[i];  
//         long int size[neighborPEs.size()]; 
//         // rcv the number of elements that get new label per PE 
//         MPI_Irecv(&size[n_PE], 1, MPI_LONG_INT, my_rank, 0, MPI_COMM_WORLD, &request+i);
        
//         // max size to recv = nº ghost nodes ? 
//         MPI_Irecv(&(*rcv_buffer)[n_PE], size[n_PE], MPI_LONG_INT, my_rank, 0, MPI_COMM_WORLD, &request+i);
//     }

//     MPI_Waitall(no_of_neighbor_PEs, requests, statuses);
}

void CommunicationHandler::update_data(DistributedGraph *g){
    // how to update the data ? 
    // we recv info at rcv_buffer so 
    // for( int i = 0; i < no_of_neighbor_PEs ; i++ ){
    //     int pe_id = neighborPEs[i]; 
    //     long int cnt = 0; 
    //     // for each element rcvd from a rank 
    //     for( auto n : rcv_buffer[pe_id] ){
    //         // get the global id -> turn into local index 
    //         T local_indx = g->from_ghost_global_to_index(n[cnt]);

    //         // update this vtx 
    //         (*g->ghost_nodes)[local_indx].current_label = n[cnt+1];
    //         cnt+=2; 
    //     }
    // }

    // restart buffers 
    // this->clearBuffers(); 
}
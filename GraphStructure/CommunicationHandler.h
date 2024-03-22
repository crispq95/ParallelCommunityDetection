#ifndef CommunicationHandler_H
#define CommunicationHandler_H

// #include "DistributedGraph.h"
#include "define.h"

class DistributedGraph; 

/* Class to handle the communications of the distributed graph */
class CommunicationHandler{
    // stores the neighbor PEs and their local IDs 
    std::unordered_map<int, int> neighborPEs;   // key : Neighbor process RANK  |  value : local id for s/rcv buffers 

    // Stores the information regarding which vtx are we sending to each PE 
    std::vector<std::vector<ID_T>> s_buffer, rcv_buffer; 

    // Stores the MPI requests to verify it worked properly
    std::vector<MPI_Request> send_request, recv_request;

    // rank / world_size / number of neighbor PEs of this process 
    int my_rank, world_size, no_of_neighbor_PEs;
    
    //array to keep the ghost vtx ordered 
    std::vector<int> ordered_ghost_indices;

    public: 
        CommunicationHandler();
        ~CommunicationHandler();    

        std::vector<std::vector<ID_T>> get_recv_buffer(){ return rcv_buffer; } ;
        int get_neigh_PEs(){ return no_of_neighbor_PEs; };


        void init_communications(std::vector<GhostNode> *ghost_vertices);
        // void addToSend(DistributedGraph* g, LocalNode node); 
        void add_all_to_send(std::unordered_set<int> * pe_ids, ID_T global_id, ID_T label);
        void add_label_to_send(std::unordered_set<int> * pe_ids, ID_T label); 
        void clear_buffers();
        void wait_requests(); 

        void send_recv_data(); 

        void send_data();
        void recv_data();  
        void recv_labels_data();  

        /* TO DO : Code to order ghosts / boundary nodes at the beginning -- not tested*/
        void order_ghosts(std::vector<GhostNode> *ghost_vertices); 
        void send_rcv_inactive( DistributedGraph* graph);

};
#endif
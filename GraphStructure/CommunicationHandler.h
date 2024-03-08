#ifndef CommunicationHandler_H
#define CommunicationHandler_H

#include "DistributedGraph.h"
#include "define.h"

/* Class to handle the communications of the distributed graph */
class CommunicationHandler{
    // stores the neighbor PEs and their local IDs 
    std::unordered_map<int, int> neighborPEs;   // key : Neighbor process RANK  |  value : local id for s/rcv buffers 

    // Marks a vtx to be sent within a LPA step 
    std::unordered_map<T, bool> to_be_sent; 

    // Stores the information regarding which vtx are we sending to each PE 
    std::vector<std::vector<T>> s_buffer; 

    // Stores the MPI requests to verify it worked properly
    std::vector<MPI_Request> recv_request;

    // rank / world_size / number of neighbor PEs of this process 
    int my_rank, world_size, no_of_neighbor_PEs;


    public: 
        CommunicationHandler();
        ~CommunicationHandler();

        void init_communications(std::vector<GhostNode> *ghost_vertices);
        // void addToSend(DistributedGraph* g, LocalNode node); 
        void add_to_send(std::unordered_set<int> * pe_ids, T global_id, T label);
        void clearBuffers();
        void send_recv_data(DistributedGraph* g);
};
#endif
#ifndef COMM_H
#define COMM_h

// #include "define.h"
// #include "graph_structure.h"
#include <vector>
#include <unordered_map>

#include "graph_structure.h"
#include "define.h"
#include <algorithm>

// Handle the required Information for communications 

// what we need to know ? 
    // IDs of the neighboring PEs (the PEs any of our local vtx are connected to)
    // current rank ? 
    // buffer to send data to each PE ? buff[0:NºPE-1] 
        // and to mark if a vtx is sent ? is_present_buff[0:NºPE-1] <bool> ?
    // class that allows us to store this info and more / methods for :
        // add a local vertex to a buffer to be sent 
        // send comm buffers (isend?) + rcv last step buff ? 
            // how to send data ? id:label ? label:id1,id2,..,idN
        // store recv data into ghost nodes 
            // translate id 
            // store new label 
struct nodePackage{
    LABEL_T label;
    long int node_id; 

    // bool operator==(nodePackage vtx1, nodePackage vtx2) const {
    //     return (vtx1.node_id == vtx2.node_id);
    // }
};


class CommunicationHandler{
    int my_rank, world_size;
    // std::vector<int> neighborPEs; 
    std::unordered_map<T, bool> to_be_sent; 
    std::unordered_map<int, int> neighborPEs;
    std::vector<MPI_Request> recv_request;
    int no_of_neighbor_PEs;

    std::vector<std::vector<nodePackage>> *send_buffer; // , *rcv_buffer; //, *is_sent;
    std::vector<std::vector<T>> s_buffer, rcv_buffer; 
    std::vector<int> n_ranks; 

    public: 
        // init n
        CommunicationHandler();
        ~CommunicationHandler();

        void init_communications(std::vector<GhostNode> *ghost_nodes);
        void addToSend(DistributedGraph* g, LocalNode node); 
        void clearBuffers();
        void send_recv_data(DistributedGraph* g);
        void update_data(DistributedGraph *g);
        // void recvData(); 
};
#endif
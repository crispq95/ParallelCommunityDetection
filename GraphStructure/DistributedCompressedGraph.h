#ifndef DistributedCompressedGraph_H 
#define DistributedCompressedGraph_H 

#include "define.h"
#include "DistributedGraph.h"

// define a class ID_To handle ID_The graph
class DistributedCompressedGraph : public DistributedGraph {
    private: 
        std::vector<ID_T> degrees;
    public: 
        DistributedCompressedGraph();
        ~DistributedCompressedGraph();

        std::vector<ID_T>* get_degrees() { return &degrees; };

        // overriden class methods
        void create_graph_from_METIS(std::string filename);

        void update_local_labels();
        void update_inactive_ghosts(std::vector<std::vector<ID_T>> recv_buffer);
};

#endif 
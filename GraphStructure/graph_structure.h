
#include <mpi.h> 
#include <vector.h>

typedef unsigned long int T;

struct Edge {
    T target; 
    T edge_weight;
};

struct Node{
    T id;
    T node_weight;
    T current_label; 
    // bool active; // maybe we can mark unactive nodes so they dont update anymore 
};

struct Node : LocalNode{
    T next_label; 
    std::vector<Edge> *edges; 
    bool is_boundary;
};

struct Node : GhostNode{
    int pe_id; // id of the PE it belongs to 
}; 

// define a class to handle the graph
class DistributedGraph(){
    DistributedGraph();
    ~DistributedGraph();
    public: 
        std::vector<LocalNode> local_nodes; 
        std::vector<GhostNode> ghost_nodes; 

        T no_local_vtx; 

        void set_local_vtx( T no_vtx );
        T get_local_vtx( return no_local_vtx; ); 

        // class methods
        void create_graph_from_METIS();
};
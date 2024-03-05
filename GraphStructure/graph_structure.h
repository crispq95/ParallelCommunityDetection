#ifndef GRAPH_STRUCTURE_H 
#define GRAPH_STRUCTURE_H

#include <mpi.h> 
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h> 
#include <iomanip>
#include <vector>
#include <unordered_map>

// #include "define.h"
// #include "comm.h"

typedef unsigned long int T; // unsigned long int ? 
typedef long int LABEL_T;

struct Edge {
    T target; 
    T edge_weight;
};

struct Node{
    T id;
    T node_weight;
    LABEL_T current_label; 
    // bool active; // maybe we can mark unactive nodes so they dont update anymore 
};

struct LocalNode : Node {
    LABEL_T next_label; 
    std::vector<Edge> *edges; 
    bool is_boundary;
};

struct GhostNode : Node {
    int pe_id; // id of the PE it belongs to 
}; 


// define a class to handle the graph
class DistributedGraph{
    private: 

        T no_local_vtx, no_total_edg, no_total_vtx;
        T vtx_begin, vtx_end; 

        // CommunicationHandler commH; 
    public: 
        // check where to place this 
        std::vector<LocalNode> *local_nodes; 
        std::vector<GhostNode> *ghost_nodes; 
        std::unordered_map<T,T> ghost_global_ids; 

        DistributedGraph();
        ~DistributedGraph();

        // getters and setters <- no hacen falta para esas var ?
        void set_local_vtx( T vtx ){ no_local_vtx = vtx; };
        void set_total_vtx( T vtx ){ no_total_vtx = vtx; };
        void set_total_edges( T edges ){ no_total_edg = edges; };
        void set_next_label( LABEL_T new_label, T n_id ){ 
            // int rank; 
            // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            // if(rank == 1)
            //     std::cout << "Changing label to local vtx " << (*local_nodes)[n_id].id << " to " << new_label << std::endl; 
            (*local_nodes)[n_id].next_label = new_label; 
        
        };

        T get_local_vtx(){      return no_local_vtx; }; 
        T get_total_vtx(){      return no_total_vtx; }; 
        T get_total_edges(){    return no_total_edg; }; 
        T get_vtx_begin(){      return vtx_begin;}; 
        T get_vtx_end(){        return vtx_end;  }; 

        // methods to transform IDs from local <-> global  
        T from_global_to_local(T global_id);
        T from_local_to_global(T local_id); 
        T from_ghost_global_to_index(T ghost_global_id);
        T from_local_ghost_to_index(T local_ghost_index);
        LABEL_T get_label(int index) { return (*local_nodes)[index].current_label; };
        T get_local_vtx_id(int index) { return (*local_nodes)[index].id; };

        // class methods
        void create_graph_from_METIS(std::string filename);
        const std::vector<Edge>* get_neighbors(T local_id); 
        bool is_ghost( T n_index ); 

};

#endif 
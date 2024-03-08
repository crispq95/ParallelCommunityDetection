#ifndef DistributedGraph_H 
#define DistributedGraph_H

#include "define.h"

// define a class to handle the graph
class DistributedGraph{
    private: 
        T no_local_vtx, no_total_edg, no_total_vtx;
        T vtx_begin, vtx_end; 

        std::unordered_map<T,T> ghost_global_ids; 
        std::vector<GhostNode> *ghost_vertices;     // may be better as arrays ?
        std::vector<LocalNode> *local_vertices; 
    public: 
        DistributedGraph();
        ~DistributedGraph();

        void set_local_vtx( T vtx ){ no_local_vtx = vtx; };
        void set_total_vtx( T vtx ){ no_total_vtx = vtx; };
        void set_total_edges( T edges ){ no_total_edg = edges; };
        void set_next_label( LABEL_T new_label, T n_id ){ (*local_vertices)[n_id].next_label = new_label; };

        T get_local_vtx(){      return no_local_vtx; }; 
        T get_total_vtx(){      return no_total_vtx; }; 
        T get_total_edges(){    return no_total_edg; }; 
        T get_vtx_begin(){      return vtx_begin;}; 
        T get_vtx_end(){        return vtx_end;  }; 
        std::vector<GhostNode>* get_ghost_vertices() { return ghost_vertices; }
        std::vector<LocalNode>* get_local_vertices() { return local_vertices; }
        
        LABEL_T get_label(int index) { return (*local_vertices)[index].current_label; };
        T get_local_vtx_id(int index) { return (*local_vertices)[index].id; };
        LocalNode get_local_vtx(int index) { return (*local_vertices)[index]; };


        // methods to transform IDs from local <-> global  
        T from_global_to_local(T global_id);
        T from_local_to_global(T local_id); 
        T from_ghost_global_to_index(T ghost_global_id);
        T from_local_ghost_to_index(T local_ghost_index);

        // class methods
        void create_graph_from_METIS(std::string filename);
        const std::vector<Edge>* get_neighbors(T local_id); 
        bool is_ghost( T n_index ); 
        void update_local_labels(); 
};

#endif 
#ifndef DistributedGraph_H 
#define DistributedGraph_H

#include "define.h"

// define a class ID_To handle ID_The graph
class DistributedGraph{
    private: 
        ID_T no_local_vtx, no_total_edg, no_total_vtx;
        ID_T vtx_begin, vtx_end; 

        std::unordered_map<ID_T,ID_T> ghost_global_ids; 
        std::vector<GhostNode> *ghost_vertices;     // may be better as arrays ?
        std::vector<LocalNode> *local_vertices; 
    public: 
        DistributedGraph();
        ~DistributedGraph();

        void set_local_vtx( ID_T vtx ){ no_local_vtx = vtx; };
        void set_total_vtx( ID_T vtx ){ no_total_vtx = vtx; };
        void set_total_edges( ID_T edges ){ no_total_edg = edges; };
        void set_next_label( LABEL_T new_label, ID_T n_id ){ (*local_vertices)[n_id].next_label = new_label; };

        ID_T get_local_vtx(){      return no_local_vtx; }; 
        ID_T get_total_vtx(){      return no_total_vtx; }; 
        ID_T get_total_edges(){    return no_total_edg; }; 
        ID_T get_vtx_begin(){      return vtx_begin;}; 
        ID_T get_vtx_end(){        return vtx_end;  }; 
        std::vector<GhostNode>* get_ghost_vertices() { return ghost_vertices; }
        std::vector<LocalNode>* get_local_vertices() { return local_vertices; }
        
        LABEL_T get_label(int index) { return (*local_vertices)[index].current_label; };
        ID_T get_local_vtx_id(int index) { return (*local_vertices)[index].id; };
        LocalNode get_local_vtx(int index) { return (*local_vertices)[index]; };


        // methods ID_To ID_Transform IDs from local <-> global  
        ID_T from_global_to_local(ID_T global_id);
        ID_T from_local_to_global(ID_T local_id); 
        ID_T from_ghost_global_to_index(ID_T ghost_global_id);
        ID_T from_local_ghost_to_index(ID_T local_ghost_id);

        // class methods
        void create_graph_from_METIS(std::string filename);
        const std::vector<Edge>* get_neighbors(ID_T local_id); 
        bool is_ghost( ID_T n_index ); 
        void update_local_labels(); 
        void update_ghost_labels(std::vector<std::vector<ID_T>> recv_buffer); 
};

#endif 
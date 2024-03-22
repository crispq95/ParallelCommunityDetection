#include "DistributedGraph.h"

/*
 *    Class: DistributedGraph  
 * Function: DistributedGraph
 * --------------------
 * DistributedGraph class constructor
 * 
 * -:-
 * 
 * returns: -
 */
DistributedGraph::DistributedGraph(){
    local_vertices = new std::vector<LocalNode>();
    ghost_vertices = new std::vector<GhostNode>(); 
}

/*
 *    Class: DistributedGraph  
 * Function: ~DistributedGraph
 * --------------------
 * DistributedGraph class destructor
 * 
 * -:-
 * 
 * returns: -
 */
DistributedGraph::~DistributedGraph(){
    for(auto n : (*local_vertices))
        delete n.edges;
    delete local_vertices;
    delete ghost_vertices;
}

/*
 *    Class: DistributedGraph  
 * Function: from_ghost_global_to_index
 * --------------------
 * Translates a global ghost id into a local index of the ghost_vertex array
 * 
 * ghost_global_id: Global ghost id to be translated
 * 
 * returns: the translated local index of the ghost vertex
 */
ID_T DistributedGraph::from_ghost_global_to_index( ID_T ghost_global_id){
    return ghost_global_ids[ghost_global_id]; 
}

/*
 *    Class: DistributedGraph  
 * Function: from_local_ghost_to_index
 * --------------------
 * Translates a local ghost index into a local index of the ghost_vertex array
 * 
 * local_ghost_id: Local ghost id to be translated
 * 
 * returns: the translated index of the ghost vertex
 */
ID_T DistributedGraph::from_local_ghost_to_index( ID_T local_ghost_id){
    return local_ghost_id - no_local_vtx;
}

/*
 *    Class: DistributedGraph  
 * Function: from_global_to_local
 * --------------------
 * Translates a global LOCAL id into a local index of the local_vertex array
 * 
 * global_id: Global id to be translated (MUST BE FROM A LOCAL VTX) 
 * 
 * returns: the local index of a local vertex 
 */
ID_T DistributedGraph::from_global_to_local( ID_T global_id){
    return global_id - vtx_begin;
}

/*
 *    Class: DistributedGraph  
 * Function: from_local_to_global
 * --------------------
 * Translates a local id of a local vertex to a global id
 * 
 * local_id: local id of a LOCAL vertex. The vtx must be local.
 * 
 * returns: the translated global id of the local vertex 
 */
ID_T DistributedGraph::from_local_to_global( ID_T local_id){
    return local_id + vtx_begin;
}

/*
 *    Class: DistributedGraph  
 * Function: get_neighbors
 * --------------------
 * Gets the neighbors of a local vertex 
 * 
 * local_id: local id of a local vertex
 * 
 * returns: the neighbors of a local vertex if sucessful, else NULL
 */
const std::vector<Edge>* DistributedGraph::get_neighbors( ID_T local_id){
    if( vtx_begin <= local_id < vtx_end){
        return (*local_vertices)[local_id].edges; 
    }
    return NULL;
}

/*
 *    Class: DistributedGraph  
 * Function: create_graph_from_METIS
 * --------------------
 * Loads from a file the information of a graph. 
 * 
 * filename: the path to the file to be used to load the graph
 * 
 * returns: -
 */
void DistributedGraph::create_graph_from_METIS(std::string filename){
    std::ifstream myFile(filename);
    std::string line, vtx, edgs; 
    
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (!myFile.is_open())
        exit(0);

    // get first line of graph file
    myFile >> vtx >> edgs;
    std::getline(myFile, line); // skip this line ? 

    // store the number of vertices and edges of the graph  
    no_total_vtx = stoul(vtx);
    no_total_edg = stoul(edgs); 

    ID_T my_part = no_total_vtx / world_size;
	ID_T remainder = no_total_vtx % world_size;

    no_local_vtx = ((rank == world_size-1) ? my_part+remainder : my_part); //: (remainder <= rank) ? my_part : my_part+1
    local_vertices->reserve(no_local_vtx); 

    unsigned long int total = 0, n_boundary = 0, cnt=0;

    // this is not the best way to do this -- CHANGE 
    vtx_begin = rank*my_part; 
    vtx_end=vtx_begin+no_local_vtx;
    cnt=0;
    ID_T ghost_index = 0; 
    
	// iterate over the rest of the file and create remaining CSR vectors 
    while (std::getline(myFile, line)) {
		if( total<vtx_begin){ // use a better way to skip lines ? 
			total++;
            continue;
		}
        if( total>=vtx_end)
			break; 
        
		std::istringstream iss(line);
		std::string word;

        LocalNode ln; 
        ln.is_boundary = false; 
        

        ln.id = from_global_to_local(total); 
        ln.node_weight = 1;    // default = 1 
        ln.current_label = total;
        ln.next_label = -1;
        ln.edges = new std::vector<Edge>(); // maybe store memory for edges 
        int num_neigh = 0;
        // int num_neigh = ceil(iss.rdbuf()->in_avail()/2);
        // (*ln.edges).reserve(num_neigh);

		// iterate over the words on the file 
		// the rest are its neighbours 
		while (iss >> word) {
            Edge e; 
            ID_T n_global_id = stoul(word);
            e.target = from_global_to_local(n_global_id); 
            e.edge_weight = 1; // here we need to set weight if available 
                        
            // check if current node is boundary (have neighbors outside)
            if( n_global_id < vtx_begin || n_global_id >= vtx_end ){
                ln.is_boundary = true; 
                ID_T ghost_id; // = ghost_global_ids.size()+(no_local_vtx);
                auto it = ghost_global_ids.find(n_global_id);

                if (it != ghost_global_ids.end()) {
                    // Neighbor already exists in ghost vertices
                    ghost_id = (*ghost_vertices)[it->second].id;
                } else {
                    // Create a new ghost 
                    ghost_id = ghost_global_ids.size()+no_local_vtx;
                    
                    GhostNode gn; 
                    gn.id = ghost_id;
                    gn.node_weight = 1; // for now, gotta change this 
                    gn.current_label = n_global_id;
                    gn.pe_id =  std::min(static_cast<unsigned int>(floor(static_cast<double>(n_global_id) / my_part)), static_cast<unsigned int>(world_size - 1));
                    ghost_global_ids[n_global_id] = ghost_index;
                    ghost_vertices->push_back(gn);
                    ghost_index++; 
                }
                e.target = ghost_id; 
            }
            ln.edges->push_back(e);
            num_neigh++;
		}
        // int num_neigh_test = ceil(word.in_avail()/2);
        // std::cout << "Vtx with global id " << total << " has " << num_neigh << " neighbors thus it is " << ln.active << std::endl;

        local_vertices->push_back(ln); 
        cnt++; 
        total++;
    }

    myFile.close();
}


/*  */
/*
 *    Class: DistributedGraph  
 * Function: is_ghost
 * --------------------
 * Evaluates wether a local index belongs to a ghost
 * 
 * n_index: Index of a vertex to be evaluated
 * 
 * returns: Returns true if a vtx id belongs to a ghost
 */
bool DistributedGraph::is_ghost( ID_T n_index ){
    return ( n_index < no_local_vtx ) ?  false : true ;
}

/*
 *    Class: DistributedGraph  
 * Function: update_local_labels
 * --------------------
 * Updates local labels from current to next label
 * 
 * -: -
 * 
 * returns: -
 */
void DistributedGraph::update_local_labels(){
    for (ID_T i = 0 ; i < (*local_vertices).size() ; i++ )
        (*local_vertices)[i].current_label = (*local_vertices)[i].next_label; 
}

/*
 *    Class: DistributedGraph  
 * Function: update_ghost_labels
 * --------------------
 * Updates ghost labels with the recv_buffer data 
 * 
 * recv_buffer: buffer that holds the data recieved from other PEs, each entry on recv buffer matches {global_id, label} pattern 
 * 
 * returns: -
 */
void DistributedGraph::update_ghost_labels(std::vector<std::vector<ID_T>> recv_buffer){
    for(int i = 0; i < recv_buffer.size(); i++){
        if(recv_buffer[i].size() > 0){
            for( int a=0 ; a < recv_buffer[i].size() ; a+=2 ){
                // std::cout << "RCV id : " << recv_buffer[i][a] << std::endl; 
                ID_T local_g_indx = from_ghost_global_to_index(recv_buffer[i][a]);
                
                // std::cout << "Global id : " << recv_buffer[i][a] << " label : " << recv_buffer[i][a+1] << std::endl; 

                // update label 
                (*ghost_vertices)[local_g_indx].current_label = recv_buffer[i][a+1];
            }
        }
    }
}

// NOT TESTED !!! - TO DO 
/*
 *    Class: DistributedGraph  
 * Function: update_ghost_labels_from_labels
 * --------------------
 * -
 * 
 * -: -
 * 
 * returns: -
 */
void DistributedGraph::update_ghost_labels_from_labels(std::vector<std::vector<ID_T>> recv_buffer, 
                                std::vector<int> id_order, std::unordered_map<int,int> neighborPEs){
    std::vector<int> cnt(recv_buffer.size(), 0);
    for ( auto i : id_order ){
        // get the PE rank of each element (in order) and transform to local 
        int local_pe_id = neighborPEs[(*ghost_vertices)[i].pe_id]; 

        // now go to the recv buffer of this pe and pop an element and store into ith ghost
        (*ghost_vertices)[i].current_label = recv_buffer[local_pe_id][cnt[local_pe_id]];
        cnt[local_pe_id]++; 
    }
}


    
int DistributedGraph::count_neighbor_labels( const LocalNode& local_vtx, CommunicationHandler& cm ){
    std::unordered_map<ID_T,ID_T> label_cnt; 
    LABEL_T max_label_value = 0; 
    std::vector<ID_T> max_labels; 
    int updated = 0;
    std::unordered_set<int> boundary_neighbor_PEs;

    label_cnt[local_vtx.current_label] = 1; 

    for(auto neigh : (*local_vtx.edges)){
        ID_T n_idx = neigh.target; 
        LABEL_T n_label; 

        if(local_vtx.is_boundary)
            n_label = (n_idx < no_local_vtx) ? (*local_vertices)[n_idx].current_label : (*ghost_vertices)[from_local_ghost_to_index(n_idx)].current_label;
        else 
            n_label = (*local_vertices)[n_idx].current_label;  
        
        if( n_idx > no_local_vtx ){
            int ghost_pe_id = (*ghost_vertices)[from_local_ghost_to_index(n_idx)].pe_id;
            if(boundary_neighbor_PEs.find(ghost_pe_id) == boundary_neighbor_PEs.end())
                boundary_neighbor_PEs.insert(ghost_pe_id); 
        }

        /* add value to label counter */
        label_cnt[n_label] += 1; 

        // if new label = maximal 
        if(max_label_value < label_cnt[n_label]){
            max_label_value = label_cnt[n_label]; 
            max_labels.clear(); // resize to 0 ? 
            max_labels.push_back(n_label); // append new max label 
        // if new label is the same value another one found previously
        }else if ( label_cnt[n_label] == max_label_value ){
            max_labels.push_back(n_label);
        }
    } 

    if(max_labels.size() >= 1){
        int rng_label = std::rand() % max_labels.size();   // pick randomly one label
        if(max_labels[rng_label] != local_vtx.current_label){
            set_next_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
            if(local_vtx.is_boundary)
                cm.add_all_to_send(&boundary_neighbor_PEs, from_local_to_global(local_vtx.id), max_labels[rng_label]);
            updated++;
        }
    }
    
    return updated;   
}
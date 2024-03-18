#include "DistributedCompressedGraph.h"

/*
 *    Class: DistributedCompressedGraph  
 * Function: DistributedCompressedGraph
 * --------------------
 * DistributedCompressedGraph class constructor
 * 
 * -:-
 * 
 * returns: -
 */
DistributedCompressedGraph::DistributedCompressedGraph(){}

/*
 *    Class: DistributedCompressedGraph  
 * Function: ~DistributedCompressedGraph
 * --------------------
 * DistributedCompressedGraph class destructor
 * 
 * -:-
 * 
 * returns: -
 */
DistributedCompressedGraph::~DistributedCompressedGraph(){}


/*
 *    Class: DistributedCompressedGraph  
 * Function: create_graph_from_METIS
 * --------------------
 * Loads from a file the information of a graph. 
 * 
 * filename: the path to the file to be used to load the graph
 * 
 * returns: -
 */
void DistributedCompressedGraph::create_graph_from_METIS(std::string filename){
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
                    gn.active = true; 
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
        degrees.insert(degrees.begin()+ln.id, num_neigh);
        ln.active = num_neigh < 2 ? false : true;
        // std::cout << "Vtx with global id " << total << " has " << num_neigh << " neighbors thus it is " << ln.active << std::endl;

        local_vertices->push_back(ln); 
        cnt++; 
        total++;
    }

    myFile.close();
}

/*
 *    Class: DistributedCompressedGraph  
 * Function: update_local_labels
 * --------------------
 * Updates local labels from current to next label
 * 
 * -: -
 * 
 * returns: -
 */
void DistributedCompressedGraph::update_local_labels(){
    for (ID_T i = 0 ; i < (*local_vertices).size() ; i++ )
        if((*local_vertices)[i].active)
            (*local_vertices)[i].current_label = (*local_vertices)[i].next_label; 
}

/*
 *    Class: DistributedGraph  
 * Function: update_inactive_ghosts
 * --------------------
 * Updates active state on ghosts from recieved data 
 * 
 * -: -
 * 
 * returns: -
 */
void DistributedCompressedGraph::update_inactive_ghosts(std::vector<std::vector<ID_T>> recv_buffer){
    for( int pe_idx = 0 ; pe_idx < recv_buffer.size(); pe_idx++ ){
        for ( auto ghost_global_id : recv_buffer[pe_idx]){
            int local_gh_indx = from_ghost_global_to_index(ghost_global_id);

            // std::cout << "Local vtx : " << local_gh_indx << " for ghost global id : " << ghost_global_id << std::endl; 
            (*ghost_vertices)[local_gh_indx].active = false;   // set them as inactive 
        } 
    }
}

#include "DistributedGraph.h"

DistributedGraph::DistributedGraph(){
    local_vertices = new std::vector<LocalNode>();
    ghost_vertices = new std::vector<GhostNode>(); 
}

DistributedGraph::~DistributedGraph(){
    for(auto n : (*local_vertices))
        delete n.edges;
    delete local_vertices;
    delete ghost_vertices;
}

/* This function only works if global_id belongs to this rank */
T DistributedGraph::from_ghost_global_to_index(T ghost_global_id){
    return ghost_global_ids[ghost_global_id]; // (ghost_global_ids.find(ghost_global_id) != ghost_global_ids.end()) ? 
}

/* This function only works if local_id belongs to this rank */
T DistributedGraph::from_local_ghost_to_index(T local_ghost_index){
    return local_ghost_index - no_local_vtx;
}

/* This function only works if global_id belongs to this rank */
T DistributedGraph::from_global_to_local(T global_id){
    return global_id - vtx_begin;
}

/* This function only works if local_id belongs to this rank */
T DistributedGraph::from_local_to_global(T local_id){
    return local_id + vtx_begin;
}

const std::vector<Edge>* DistributedGraph::get_neighbors(T local_id){
    if( vtx_begin <= local_id < vtx_end){
        return (*local_vertices)[local_id].edges; 
    }
    return NULL;
}

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

    T my_part = no_total_vtx / world_size;
	T remainder = no_total_vtx % world_size;

    no_local_vtx = ((rank == world_size-1) ? my_part+remainder : my_part); //: (remainder <= rank) ? my_part : my_part+1
    local_vertices->reserve(no_local_vtx); 

    unsigned long int total = 0, n_boundary = 0, cnt=0;

    // this is not the best way to do this -- CHANGE 
    vtx_begin = rank*my_part; 
    vtx_end=vtx_begin+no_local_vtx;
    cnt=0;
    T ghost_index = 0; 
    
	// iterate over the rest of the file and create remaining CSR vectors 
    while (std::getline(myFile, line)) {
		if(total<vtx_begin){ // use a better way to skip lines ? 
			total++;
            continue;
		}
        if(total>=vtx_end)
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
        (*ln.edges).reserve(ceil(iss.rdbuf()->in_avail()/2));

		// iterate over the words on the file 
		// the rest are its neighbours 
		while (iss >> word) {
            Edge e; 
            T n_global_id = stoul(word);
            e.target = from_global_to_local(n_global_id); 
            e.edge_weight = 1; // here we need to set weight if available 
                        
            // check if current node is boundary (have neighbors outside)
            if( n_global_id < vtx_begin || n_global_id >= vtx_end ){
                ln.is_boundary = true; 
                T ghost_id; // = ghost_global_ids.size()+(no_local_vtx);
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
		}

        local_vertices->push_back(ln); 
        cnt++; 
        total++;
    }

    myFile.close();
}

/* Returns true if a vtx id belongs to a ghost */
bool DistributedGraph::is_ghost( T n_index ){
    return ( n_index < no_local_vtx ) ?  false : true ;
}

/* Updates local labels from current to next label */
void DistributedGraph::update_local_labels(){
    for( T i = 0 ; i < (*local_vertices).size() ; i++ )
        (*local_vertices)[i].current_label = (*local_vertices)[i].next_label; 
}

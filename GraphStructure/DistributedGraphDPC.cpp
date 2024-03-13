#include "DistributedGraph.h"
#include "DistributedGraphDPC.h"

/*
 *    Class: DistributedGraphDPC
 * Function: DistributedGraphDPC
 * --------------------
 * DistributedGraph class constructor
 * 
 * -:-
 * 
 * returns: -
 */
DistributedGraphDPC::DistributedGraphDPC(){}

/*
 *    Class: DistributedGraphDPC
 * Function: ~DistributedGraphDPC
 * --------------------
 * DistributedGraphDPC class destructor
 * 
 * -:-
 * 
 * returns: -
 */
DistributedGraphDPC::~DistributedGraphDPC(){}

// test funct 
void DistributedGraphDPC::create_graph_for_DCP_from_METIS(std::string filename){
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
    degree.reserve(no_local_vtx);
    
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
        ln.active = num_neigh < 2 ? false : true;
        degree[ln.id] = num_neigh; 
        // std::cout << "Vtx with global id " << total << " has " << num_neigh << " neighbors thus it is " << ln.active << std::endl;

        local_vertices->push_back(ln); 
        cnt++; 
        total++;
    }

    myFile.close();
}

// quality = weight_vtx -> our case = 1 
void DistributedGraphDPC::calculate_quality(){
    // get max value for degree (of this PE):
    // max_weight = std::max_element(degrees.begin(), degrees.end()) 

    // get max GLOBAL value 
    // MPI_Allreduce(MPI_IN_PLACE, &max_weight, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 

    for ( int i = 0; i < no_local_vtx ; i++ ) 
        quality[i] = 1; // pruebas 
}

// density (normalized) = deg/max_deg 
void DistributedGraphDPC::calculate_density(){
    // get max value for degree (of this PE):
    max_degree = std::max_element(degrees.begin(), degrees.end()) 

    // get max GLOBAL value 
    MPI_Allreduce(MPI_IN_PLACE, &max_degree, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 

    // density.reserve();
    for ( auto deg : degrees ) {
        density.push_back(deg/max_deg);
    }
}

// centrality index = qual * dens -> vertex with large centrality index = comm_seeds 
void DistributedGraphDPC::calculate_centrality_index(){
    for ( int i = 0; i < no_local_vtx ; i++ ) {
        centrality_index.push_back(density[i]*quality[i]); 
    }
}

void DistributedGraphDPC::construct_second_order_diff_decreasing_sequence(){
    sort_indexes_by_label(&centrality_index, &decreasing_sequence);

    for( int i = 1; i < no_local_vtx ; i++ ){
        second_order_difference[i] = fabs( ( centrality_index[decreasing_sequence[i]] - centrality_index[decreasing_sequence[i+1]]  ) - ( centrality_index[decreasing_sequence[i+1]] - centrality_index[decreasing_sequence[i+2]] ) );
    }
}

void DistributedGraphDPC::draw_potential_seeds(){
    

    
}

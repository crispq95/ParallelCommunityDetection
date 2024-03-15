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
    degrees.reserve(no_local_vtx);
    
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
        degrees.insert(degrees.begin()+ln.id, num_neigh); 
        // std::cout << "Vtx with global id " << total << " has " << num_neigh << " neighbors thus it is " << ln.active << std::endl;

        local_vertices->push_back(ln); 
        cnt++; 
        total++;
    }

    myFile.close();
}

/*
 *    Class: DistributedGraphDPC
 * Function: calculate_quality
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
// quality = weight_vtx -> our case = 1 
void DistributedGraphDPC::calculate_quality(){
    // get max value for degree (of this PE):
    // max_weight = std::max_element(degrees.begin(), degrees.end()) 

    // get max GLOBAL value 
    // MPI_Allreduce(MPI_IN_PLACE, &max_weight, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    quality.reserve(no_local_vtx); 
    for ( int i = 0; i < no_local_vtx ; i++ ) 
        quality.push_back(1); // pruebas 
}

/*
 *    Class: DistributedGraphDPC
 * Function: calculate_density
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
// density (normalized) = deg/max_deg 
void DistributedGraphDPC::calculate_density(){
    // get max value for degree (of this PE):
    max_degree = (*std::max_element(degrees.begin(), degrees.end()));

    // get max GLOBAL value 
    MPI_Allreduce(MPI_IN_PLACE, &max_degree, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 

    density.reserve(no_local_vtx); 
    for ( auto deg : degrees ) {
        density.push_back(deg/(double)max_degree); 
    }
}

/*
 *    Class: DistributedGraphDPC
 * Function: calculate_centrality_index
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
// centrality index = qual * dens -> vertex with large centrality index = comm_seeds 
void DistributedGraphDPC::calculate_centrality_index(){
    centrality_index.reserve(no_local_vtx);
    for ( int i = 0; i < no_local_vtx ; i++ ) {
        centrality_index.push_back(density[i]*quality[i]); 
    }
}

/*
 *    Class: DistributedGraphDPC
 * Function: construct_second_order_diff_decreasing_sequence
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void DistributedGraphDPC::construct_second_order_diff_decreasing_sequence(){
    decreasing_sequence.reserve(no_local_vtx);
    second_order_difference.reserve(no_local_vtx);  // -2 ? 
 
    sort_decreasing_indexes(centrality_index, decreasing_sequence);

    for( int i = 0; i < no_local_vtx-2 ; i++ ){
        second_order_difference.push_back(fabs( ( centrality_index[decreasing_sequence[i]] - centrality_index[decreasing_sequence[i+1]]  ) - ( centrality_index[decreasing_sequence[i+1]] - centrality_index[decreasing_sequence[i+2]] ) ));
        // printf("%f for vtx %ld \n", second_order_difference[i], from_local_to_global(decreasing_sequence[i]));
    }
}

/*
 *    Class: DistributedGraphDPC
 * Function: get_number_seeds
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void DistributedGraphDPC::get_number_seeds(){
    // get the number of seeds : argmax h <- first max h 
    /* OPTION 1 : */
    num_seeds = 0;
    double max = 0;
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    for( int i = 0; i < no_local_vtx ; i ++ ){
        double max_candidate = second_order_difference[i];
        if( max < max_candidate ){
            max = max_candidate;
            num_seeds = i;

            // if(rank == 2 )
            //     std::cout << " Rank  " << rank << " sending value " << max << " with i : " << from_local_to_global(decreasing_sequence[num_seeds]) << std::endl;
        }
    } 

    // get the first maximum (lower arg_max among PEs)
        // ta mal -> hay que buscar el max de las que hay y su indice -> MPI_MAXLOC ? 
    // create datatype maxloc
    struct { 
        double val; 
        ID_T   i; 
    } num_max; 

    num_max.val = max;
    num_max.i = from_local_to_global(decreasing_sequence[num_seeds]);



    // std::cout << " Rank  " << rank << " sending value " << num_max.val << " with i : " << num_max.i << std::endl;
    MPI_Allreduce(MPI_IN_PLACE, &num_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); 
    // std::cout << " Rank  " << rank << " RCV value " << num_max.val << " with i : " << num_max.i << std::endl;

    num_seeds=0;
    // this is not exactly what it should be -> we should send a portion of the vtx and then decide which ones are 
    // the seeds
    //locates current seeds  
    for( auto a : decreasing_sequence ){
        if ( std::round(second_order_difference[a] * 100) >= std::round(num_max.val * 100)){
            seed_candidates.push_back(from_local_to_global((*local_vertices)[decreasing_sequence[a]].id)); 
            num_seeds++;         
            std::cout << "At rank " << rank << " my vtx " << from_local_to_global((*local_vertices)[decreasing_sequence[a]].id) << " is a seed [" << centrality_index[a] << "]." << std::endl; 
        }
    }

    // on processes where there is a candidate seed send it to others 
    // MPI_AllGather(MPI_IN_PLACE, &num_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); 
}

/*
 *    Class: DistributedGraphDPC
 * Function: get_seeds
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void DistributedGraphDPC::get_seeds(){
    /* For testing */
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /* END - For testing */

    // get normalized quality and density; 
    calculate_quality();
    calculate_density();

    // get centrality index and construct 2nd order diff decreasing sequence 
    calculate_centrality_index();
    construct_second_order_diff_decreasing_sequence(); 

    // // get the number of seeds and reserve space for them in this PE 
    get_number_seeds();
    // community_seeds.reserve(num_seeds); 

    // // get any seed below this number ? -> should order seeds first somehow ?
    // // in sequential -> put on queue your vtx in order of g
    // //                  add a vtx to seeds if its not neighbor of a current seed 
    // //                  remove from queue its neighbors ?
    // std::vector<bool> seeds_found(no_local_vtx, 0);
    // int i = 0, found_seeds; 
    
    // // do similar ? order, decide which ones, send to other PEs candidates and draw results 
    // while( found_seeds < num_seeds || i >= no_local_vtx ){
    //     if( seeds_found[decreasing_sequence[i]] == 0 ){
    //         community_seeds.push_back(decreasing_sequence[i]); 
    //         for(auto neigh : (*(*local_vertices)[decreasing_sequence[i]].edges)){
    //             if( neigh.target < no_local_vtx){ // if its local 
    //                 seeds_found[neigh.target] = 1;
    //                 found_seeds++;
    //             }
    //         }  
    //     } 
    //     i++;
    // }

    // if(rank == 0)
    //     std::cout << " Found " << found_seeds << " of " << num_seeds << std::endl;
    // std::cout << " Candidates on rank " << rank << ": ";
    // for( auto a : community_seeds )
    //     std::cout << a << " " ;
    // std::cout << std::endl; 

    // // here we have local vtx that are seed candidates -> send to all and decide for candidate_smaller seeds 
    //     // should prob send Global_ID + centrality_index of candidates
    // // MPI_Bcast(&community_seeds[0], num_seeds, MPI_UNSIGNED_LONG, rank, MPI_COMM_WORLD); 
    // std::vector<ID_T> tmp_seed_buff(num_seeds*rank,0);
    // MPI_Allgather(&community_seeds[0], num_seeds, MPI_UNSIGNED_LONG, &tmp_seed_buff[0], num_seeds, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    // // order the recv array by centrality_index and decide final seeds <- all proc do this 
    
    //     // so they all get the same seeds 
    // community_seeds.fill(0);
    // seeds_found.resize(num_seeds*rank);
    // seeds_found.fill(0);

    // // while( found_seeds >= num_seeds || i >= no_local_vtx ){
    // //     if( seeds_found[decreasing_sequence[i]] == 0 ){
    // //         community_seeds.push_back(decreasing_sequence[i]); 
    // //         for(auto neigh : (*(*local_vertices)[decreasing_sequence[i]].edges)){
    // //             if( neigh.target < no_local_vtx) // if its local 
    // //                 seeds_found[neigh.target] = 1;
    // //         }  
    // //     } 
    // //     i++;
    // // }

}

/*
 *    Class: DistributedGraphDPC
 * Function: find_cores
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void DistributedGraphDPC::find_cores(){
    // for now we ignore this, no need to find cores ?
    // for ( auto cs : community_seeds ){
    //     // get neighbors of a seed : 
    //     for ( auto neigh : (*cs.edges) ){
    //         // for each neighbor if not neighb to other seed : add to CORE of this seed [change label to the seeds]
    //         // if its local 
    //         if( neigh.target < no_vtx_local ){
    //             // get the id and label 
    //             ID_T neigh_id = neigh.target;
    //             LABEL_T neigh_lab = (*local_vertices)[neigh.target].current_label; 

    //             // no inicializada 
    //             if(neigh_lab == -1){
    //                 if(find(community_seeds.begin(), community_seeds.end(), (*local_vertices)[neigh.target].id) == 1){
    //                     // add to CORE 
    //                     (*local_vertices)[neigh.target].current_label = cs.label;
    //                 }
    //             }
    //         }
    //     }
    // }

    // Just run a single LPA step blocking seeds and if a draw happens -> dont pick a label 

}

/*
 *    Class: DistributedGraphDPC
 * Function: run_DPC_LPA_step
 * --------------------
 * -
 * 
 * -:-
 * 
 * returns: -
 */
void DistributedGraphDPC::run_DPC_LPA_step(){

}
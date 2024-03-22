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
        ln.current_label = -1;
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
                    // gn.current_label = n_global_id;
                    gn.current_label = -1; // label is now -1 -> no node is initialized 
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
    // max_degree = (*std::max_element(degrees.begin(), degrees.end()));

    // std::cout << "MAX DEG : " << max_degree << std::endl;
    // get max GLOBAL value 
    // MPI_Allreduce(MPI_IN_PLACE, &max_degree, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 

    // calculate avg degree 
    // double avg_deg = 0; 
    // density.reserve(no_local_vtx); 
    // for ( auto deg : degrees ) {
    //     density.push_back(deg/(double)max_degree); 
    //     avg_deg += deg; 
    // }
    
    // MPI_Allreduce(MPI_IN_PLACE, &avg_deg, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    // avg_deg = avg_deg/no_total_vtx;

    // int cnt_higer_avg=0;
    // for ( auto deg : degrees ) {
    //     if(avg_deg*2 < deg)
    //         cnt_higer_avg++;    
    // }
    // std::cout << "[from a total of : "<< no_local_vtx <<" /"<< no_total_vtx << "] Nº of higher than avg nodes :  " << cnt_higer_avg << std::endl;
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
    }
}

// /*
//  *    Class: DistributedGraphDPC
//  * Function: get_number_seeds
//  * --------------------
//  * - USE LATER 
//  * 
//  * -:-
//  * 
//  * returns: -
//  */
// void DistributedGraphDPC::get_number_seeds(){
//     // // get the number of seeds : argmax h <- first max h 
//     // /* OPTION 1 : */
//     // num_seeds = 0;
//     // double max = 0;
//     // int world_size, rank;
//     // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//     // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     // for( int i = 0; i < no_local_vtx-2 ; i ++ ){
//     //     double max_candidate = second_order_difference[i];
//     //     if( max < max_candidate ){
//     //         max = max_candidate;
//     //         num_seeds = i;
//     //     }
//     // } 

//     // // get the first maximum (lower arg_max among PEs)
//     //     // ta mal -> hay que buscar el max de las que hay y su indice -> MPI_MAXLOC ? 
//     // // create datatype maxloc
//     // struct { 
//     //     double val; 
//     //     ID_T   i; 
//     // } num_max; 

//     // num_max.val = max;
//     // num_max.i = from_local_to_global(decreasing_sequence[num_seeds]);

//     // // std::cout << " Rank  " << rank << " sending value " << num_max.val << " with i : " << num_seeds << std::endl;
//     // MPI_Allreduce(MPI_IN_PLACE, &num_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD); 
//     // // std::cout << " Rank  " << rank << " RCV value " << num_max.val << " with i : " << num_seeds << std::endl;

//     // // num_seeds=0;
//     // // // this is not exactly what it should be -> we should send a portion of the vtx and then decide which ones are 
//     // // // the seeds
//     // // //locates current seeds  
//     // // for( auto a : decreasing_sequence ){
//     // //     if ( std::round(second_order_difference[a] * 1000) >= std::round(num_max.val * 1000)){
//     // //         // seed_candidates.push_back(from_local_to_global((*local_vertices)[decreasing_sequence[a]].id)); 
//     // //         num_seeds++;         
//     // //         // std::cout << "At rank " << rank << " my vtx " << from_local_to_global((*local_vertices)[decreasing_sequence[a]].id) << " is a seed [" << centrality_index[a] << "]." << std::endl; 
//     // //     }
//     // // }

//     // // on processes where there is a candidate seed send it to others 
// }

// /*
//  *    Class: DistributedGraphDPC
//  * Function: get_seeds
//  * --------------------
//  * Provisional -- consigue seeds ordenando las disponibles + eligiendo mejores de entre el nº que se necesitan
//  * 
//  * -:-
//  * 
//  * returns: -
//  */
 void DistributedGraphDPC::get_seeds( int num_seeds, CommunicationHandler* cm ){
    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    calculate_density(); 
    // community_seeds.resize(num_seeds); 
    
    std::vector<ID_T> community_seed_degrees(num_seeds); 
    int size_send = num_seeds/world_size, rcv_size=size_send; 

    // get candidate seeds of this PE 
    decreasing_sequence.reserve(no_local_vtx);
    ck_distance_weights.resize(no_local_vtx, -1);
    // std::fill(ck_distance_weights.begin(), ck_distance_weights.end(), -1);
    sort_decreasing_indexes(degrees, decreasing_sequence);


    // for( int i = 0; i < degrees.size(); i++){
    //     std::cout << degrees[decreasing_sequence[i]] << ", " << from_local_to_global(decreasing_sequence[i]) << " -- " << std::endl;
    // }

    /* set labels to seed candidates on each PE   */
    std::vector<bool> is_candidate(no_local_vtx, true);

    // rmv candidates that are neighbors -> set them as cores : 
    int candidate = 0, found_communities=0;
    while(found_communities < num_seeds && candidate < no_local_vtx){
        ID_T candidate_seed = decreasing_sequence[candidate];
        candidate++;
        bool boundary_candidate = (*local_vertices)[candidate_seed].is_boundary;
        // std::cout  << rank << " ["<< from_local_to_global(candidate_seed) << " from " << no_local_vtx << "] .. found : " << found_communities << std::endl;

        if(!is_candidate[candidate_seed])
            continue; 

        (*local_vertices)[candidate_seed].next_label = from_local_to_global(candidate_seed); 
        // ck_distance_weights.insert(ck_distance_weights.begin()+candidate_seed, (DEFAULT_WEIGTH + SEED_WEIGTH + 1));
        ck_distance_weights[candidate_seed] = double(DEFAULT_WEIGTH + SEED_WEIGTH + 1);  
        candidate_queue.push(candidate_seed);

        // if it is boundary -> add to send [label updated!]
        found_communities++; 
        std::unordered_set<int> boundary_neighbor_PEs;

        // remove neighbor candidate from beign a possible seed 
        for( auto neigh : (*(*local_vertices)[candidate_seed].edges) ){
            // if( vtx_begin <= neigh.target && neigh.target < vtx_end ){
            if( !is_ghost(neigh.target) ){
                is_candidate[neigh.target] = false; 
            }else if (boundary_candidate){
                // add PE of neighbor ghost to list :  
                int ghost_pe_id = (*ghost_vertices)[from_local_ghost_to_index(neigh.target)].pe_id;
                
                if(boundary_neighbor_PEs.find(ghost_pe_id) == boundary_neighbor_PEs.end())
                    boundary_neighbor_PEs.insert(ghost_pe_id);
            }
        }
        
        if(boundary_candidate){   // only necessary to send candidate seed id (label == id) <- modify later
            cm->add_all_to_send(&boundary_neighbor_PEs, from_local_to_global(candidate_seed), from_local_to_global(candidate_seed));
        }
    }
    
    cm->send_data(); 
    update_local_labels(); 
    cm->recv_data(); 
    update_ghost_labels(cm->get_recv_buffer());


    // std::cout << "Printing ghost labels on "<< rank << " : " << std::endl;
    // for(auto gh : ghost_global_ids ){
    //     ID_T g_id = gh.first, l_idx = gh.second; 

    //     std::cout << g_id << "," << (*ghost_vertices)[l_idx].current_label << " , " << l_idx << "] || "; 
    // }   


    // std::cout << std::endl; 

    // update all labels to next 

    // std::cout<< "Out from " << rank << std::endl; 

    // printf("Candidates on rank %d are : ", rank);
    // while(!candidate_queue.empty()){
    //     ID_T q = candidate_queue.front(); 
    //     candidate_queue.pop(); 
    //     printf("%ld [deg %ld] - ",  from_local_to_global((*local_vertices)[q].id), degrees[(*local_vertices)[q].id]);
    // }

    // for( int i = 0; i < no_local_vtx; i++ ){
    //     if((*local_vertices)[i].next_label != -1)
    //         printf("%ld [deg %ld] - ",  from_local_to_global((*local_vertices)[i].id), degrees[(*local_vertices)[i].id]);
    // }

    // printf("\n");

    // std::cout << "Nº elem : " << ck_distance_weights.size() << std::endl; 

    // std::vector<ID_T> tmp_degree, tmp_ids; 
    // tmp_degree.reserve(no_local_vtx);  
    // tmp_ids.reserve(no_local_vtx); 

    // for( int i = 0; i < size_send; i++ ){
    //     tmp_degree.push_back(degrees[decreasing_sequence[i]]);
    //     tmp_ids.push_back(from_local_to_global(decreasing_sequence[i])); 
    // }

    // // send-recv candidates from other PEs
    // MPI_Allgather(&tmp_degree[0] , size_send, MPI_UNSIGNED_LONG, &community_seed_degrees[0], rcv_size, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    // MPI_Allgather(&tmp_ids[0], size_send, MPI_UNSIGNED_LONG, &community_seeds[0], rcv_size, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    // // order arrays by degree
    // decreasing_sequence.clear(); 
    // decreasing_sequence.resize(num_seeds); 
    // sort_decreasing_indexes(tmp_degree, decreasing_sequence);

    // // decide which ones are seeds 
    // int found_communities = 0, candidate = 0; 
    // while(found_communities < num_seeds && candidate < num_seeds){
    //     // get current candidate 
    //     ID_T candidate_seed = community_seeds[decreasing_sequence[candidate]];
    //     candidate++;
    //     // if candidate cannot be a seed -> skip this candidate 
    //     if(label_counters[candidate_seed] == 1){
    //         continue; // go next 
    //     }
    //     if( candidate_seed < vtx_begin || tmp_ids >= vtx_end) {
    //         (*local_vertices)[candidate_seed].next_label = from_global_to_local(candidate_seed); 

    //     }else {
    //         (*ghost_vertices)[candidate_seed].current_label = from_local_ghost_to_index(ghost_global_ids[candidate_seed]); 
    //     }
    //     // else -> candidate is a seed 
    //     // mark it as seed (membership = candidate)
    //     ck_distance_weights.push_back(DEFAULT_WEIGTH + SEED_WEIGTH + 1);
    //     seed_nodes.push_back(candidate_seed);

    //     for(int i=0; i<igraph_vector_size(&seed_neighbors); i++){
    //         int v = VECTOR(seed_neighbors)[i];
    //         // remove all neighbors from candidate list 
    //         VECTOR(label_counters)[v] += 1;   
    //     }

    //     found_communities++;
    // }


//   printf("FOUND COMMUNITIES : %d\n", found_communities); 
    
    // init the weight of seeds of this PE 
 }
 

// /*
//  *    Class: DistributedGraphDPC
//  * Function: find_cores
//  * --------------------
//  * -
//  * 
//  * -:-
//  * 
//  * returns: -
//  */
void DistributedGraphDPC::find_cores(){
    // set the ghost nodes that are seeds : 
    std::unordered_set<ID_T> ghost_candidate_seeds;
    
    for( auto g : (*ghost_vertices)){
        if( g.current_label != -1 )
            ghost_candidate_seeds.insert(g.id);
    }

    while(!candidate_queue.empty()){
        ID_T q = from_global_to_local(candidate_queue.front());
        candidate_queue.pop();
        bool found = false; // if not found -> the seed adds its neighbors to the queue 

        // for the neighbors of the 
        for( auto neigh : (*(*local_vertices)[q].edges)){
            ID_T neighbo_local_id = neigh.target;
            // they may be a candidate only if the neighbor is a ghost 
            if(is_ghost(neighbo_local_id)){
                //if it is a ghost, check if it is on the set of candidates 
                if( ghost_candidate_seeds[neighbo_local_id] != ghost_candidate_seeds.end() ){
                    // keep the label of the vertex with more OR keep in candidates until further notice 
                    found = true; 
                }
            }else{
                local_vertices[neighbo_local_id].next_label = local_vertices[q].current_label;
            }
        }

        // if it is neighbor to a single seed : 
        if(!found){
            /* add all neighbors to community CORE  */
            // set the label of the neighbors 
            
        // if it is neighbor to more than one seed : 
        }else {

        }
    }

}


// /*
//  *    Class: DistributedGraphDPC
//  * Function: run_DPC_LPA_step
//  * --------------------
//  * -
//  * 
//  * -:-
//  * 
//  * returns: -
// //  */
// // void DistributedGraphDPC::run_DPC_LPA_step( int nstep ){
// //     // For all the candidates 
// //     for( auto local_vtx : current_candidates ){
// //         std::unordered_map<ID_T,ID_T> label_cnt; 
// //         LABEL_T max_label_value = 0; 
// //         std::vector<ID_T> max_labels; 
// //         int updated = 0;
// //         std::unordered_set<int> boundary_neighbor_PEs;

// //         label_cnt[local_vtx.current_label] = distance_weight[local_vtx.idx];  

// //         // for neighbors of candidate 
// //         for( auto neigh : (*local_vertices)[local_vtx.id] ){
// //             ID_T n_idx = neigh.target; 
// //             LABEL_T n_label = -1; 

// //             if(local_vtx.is_boundary)
// //                 n_label = (n_idx < no_local_vtx) ? (*local_vertices)[n_idx].current_label : (*ghost_vertices)[from_local_ghost_to_index(n_idx)].current_label;
// //             else 
// //                 n_label = (*local_vertices)[n_idx].current_label;  
            
// //             if( n_label == -1 ) // if the vtx has no label assigned : 
// //                 continue; 

// //             if( n_idx > no_local_vtx ){
// //                 int ghost_pe_id = (*ghost_vertices)[from_local_ghost_to_index(n_idx)].pe_id;
// //                 if(boundary_neighbor_PEs.find(ghost_pe_id) == boundary_neighbor_PEs.end())
// //                     boundary_neighbor_PEs.insert(ghost_pe_id); 
// //             }

// //             /* add value to label counter */
// //             label_cnt[n_label] += distance_weight[neigh.idx]; // sum the distance weight to the label_counter  

// //             // if new label = maximal 
// //             if(max_label_value < label_cnt[n_label]){
// //                 max_label_value = label_cnt[n_label]; 
// //                 max_labels.clear(); // resize to 0 ? 
// //                 max_labels.push_back(n_label); // append new max label 
// //             // if new label is the same value another one found previously
// //             }else if ( label_cnt[n_label] == max_label_value ){
// //                 max_labels.push_back(n_label);
// //             }
// //         }

// //         // calculate weight of the assigned candidate :   
// //         if(max_labels.size() >= 1){
// //             int rng_label = std::rand() % max_labels.size();   // pick randomly one label
// //             if(max_labels[rng_label] != local_vtx.current_label){
// //                 set_next_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
// //                 if(local_vtx.is_boundary)
// //                     cm.add_all_to_send(&boundary_neighbor_PEs, from_local_to_global(local_vtx.id), max_labels[rng_label]);

// //                 /* Add the neighbors of the labeled vtx */

// //                 /* Add the distance weight to the labeled vtx */

// //                 updated++;
// //             }
// //         }
// //     }
// // }
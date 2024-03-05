
#include "graph_structure.h"

DistributedGraph::DistributedGraph(){
    local_nodes = new std::vector<LocalNode>();
    ghost_nodes = new std::vector<GhostNode>(); 
}

DistributedGraph::~DistributedGraph(){
    for(auto n : (*local_nodes))
        delete n.edges;
    delete local_nodes;
    delete ghost_nodes;
}

/* TO DO */
/* This function only works if global_id belongs to this rank */
T DistributedGraph::from_ghost_global_to_index(T ghost_global_id){
    return ghost_global_ids[ghost_global_id]; // (ghost_global_ids.find(ghost_global_id) != ghost_global_ids.end()) ? 
}

/* This function only works if local_id belongs to this rank */
T DistributedGraph::from_local_ghost_to_index(T local_ghost_index){
    return local_ghost_index - no_local_vtx;
}
/* END TO DO */

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
        return (*local_nodes)[local_id].edges; 
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

    // std::cout << "I'm on process " << rank << " with " << no_total_vtx << " total vtx, my part is " << my_part << " and the remainder is " << remainder << std::endl; 

    // revisar 
    no_local_vtx = ((rank == world_size-1) ? my_part+remainder : my_part); //: (remainder <= rank) ? my_part : my_part+1
    local_nodes->reserve(no_local_vtx); 

    unsigned long int total = 0, n_boundary = 0, cnt=0;

    // this is not the best way to do this -- CHANGE 
    vtx_begin = rank*my_part; 
    vtx_end=vtx_begin+no_local_vtx;
    cnt=0;
    T ghost_index = 0; 
    
    // std::cout <<  "I'm on process " << rank << " begining at " << vtx_begin << " and ending at " << vtx_end <<  std::endl;
    // std::cout << "I'm on process " << rank << " with " << no_local_vtx << " assigned, starting at " << vtx_begin << " and going until " << vtx_end << std::endl; 
	
	// iterate over the rest of the file and create remaining CSR vectors 
    while (std::getline(myFile, line)) {
		std::istringstream iss(line);
		std::string word;
		if(total<vtx_begin){ // use a better way to skip lines ? 
			total++;
            continue;
		}if(total>=vtx_end)
			break; 
        
        LocalNode ln; 
        ln.is_boundary = false; 

        // ln.id = total; // id global 
        ln.id = from_global_to_local(total); 
        ln.node_weight = 1;    // default = 1 
        // ln.current_label = ln.id;  // Labeled to the id of the vertex
        ln.current_label = total;
        ln.next_label = -1;
        ln.is_boundary = false; 
        ln.edges = new std::vector<Edge>(); // maybe store memory for edges 

		// iterate over the words on the file 
		// the rest are its neighbours 
		while (iss >> word) {
            Edge e; 
            // e.target = stoul(word); // store globals
            e.target = from_global_to_local(stoul(word)); 
            e.edge_weight = 1; // here we need to set weight if available 
                        
            // check if current node is boundary (have neighbors outside)
            if( stoul(word) < vtx_begin || stoul(word) >= vtx_end ){
                ln.is_boundary = true; 
                // T ghost_id = ghost_global_ids.size()+vtx_end;
                T ghost_id = ghost_global_ids.size()+(no_local_vtx);
                
                // std::cout << "I'm getting this number : " << ghost_index << " from " << (*ghost_nodes).size() << " + " << no_local_vtx << std::endl;
                if (ghost_global_ids.find(stoul(word)) != ghost_global_ids.end()){
                    // get its id and set the neighbor of current node 
                    ghost_id = (*ghost_nodes)[ghost_global_ids[stoul(word)]].id;
                    // if(rank == 1)
                        // std::cout << "I'm getting this number : " << ghost_id << " from " << stoul(word) << std::endl;
                }else{
                    GhostNode gn; 
                    gn.id = ghost_id;
                    gn.node_weight = 1; // for now, gotta change this 
                    // gn.current_label = e.target; // conserva su label original 
                    gn.current_label = stoul(word);
                    // gn.pe_id = ceil((double)(stoul(word))/(double)my_part)-1;    // calculate pe id 
                    unsigned int rank_neigh_PE = floor((double)(stoul(word))/(double)my_part) >= world_size ? (world_size-1) : floor((double)(stoul(word))/(double)my_part);
              
                    // if this rank is the local index instead of the rank everything is easier ?
                        // other than that, need a conversion from global rank to local rank index 
                    gn.pe_id =  rank_neigh_PE;   // calculate pe id 


                    // if(rank == 0)
                    //     std::cout << "[" << world_size << "] Ghost vtx with id " << gn.id << " belongs to PE " << gn.pe_id << " -- " << stoul(word) << "/" << my_part << ": [" << floor((double)(stoul(word))/(double)my_part) << "]" << std::endl;
                    // if(rank == 1)
                    //     std::cout << "[g_id :" << ghost_id << "]" << stoul(word) <<" / " << my_part << " = " << (double)(stoul(word))/(double)my_part << " redondeado a : " << floor((double)(stoul(word))/(double)my_part)<< std::endl;
                    ghost_global_ids[stoul(word)] = ghost_index;
                    ghost_nodes->push_back(gn);
                    ghost_index++; 
                }
                e.target = ghost_id; 
            }

            ln.edges->push_back(e);

		}

        local_nodes->push_back(ln); 
        cnt++; 
        total++;
    }

    // // TESTING : 
    // if (rank == 1){
    //     for (auto n : (*ghost_nodes))
    //         std::cout << "ID : " << n.id << " - PE_origen : " << n.pe_id << std::endl; 
    //         // std::cout << "in " << rank << " vertex with: ID : " << n.id << ",label "<< n.next_label << ",boundary " << n.is_boundary << ",nº of edges " << n.edges->size() << " : ";
    //         // for(auto e : (*n.edges)){
    //         //     std::cout << e.target << " " ;
    //         // }
    //         std::cout << std::endl;

    //         // std::cout << "These are locally named : " << std::endl; 
            
    //         // for (auto i = ghost_global_ids.begin(); i != ghost_global_ids.end(); i++) {
    //         //     T index = i->second;

    //         //     std::cout << i->first << " \t\t\t" << i->second << std::endl; 
    //         //     std::cout << "[" << (*ghost_nodes).size() << "] Ghost with hash index " << index << " has value : " << (*ghost_nodes)[index].id << std::endl; 

    //         //     // index = 19; 
    //         //     // T local_index = from_local_ghost_to_index(index);
    //         //     // std::cout << "Local index for ghost : " << index << "translates to ghost :" << (*ghost_nodes)[local_index].id << std::endl;
    //         // }
                

    //         // std::cout << std::endl;


    //     }

    myFile.close();
}

/* Returns true if a vtx id belongs to a ghost */
bool DistributedGraph::is_ghost( T n_index ){
    return ( n_index < no_local_vtx ) ?  false : true ;
}

void DistributedGraph::update_labels(){
    for( T i = 0 ; i < (*local_nodes).size() ; i++ )
        (*local_nodes)[i].current_label = (*local_nodes)[i].next_label; 
}
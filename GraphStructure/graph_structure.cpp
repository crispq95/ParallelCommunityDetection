
// DistributedGraph::DistributedGraph(){

// }

// DistributedGraph::~DistributedGraph(){}; 

void DistributedGraph::create_graph_from_METIS(std::string filename, int rank, int world_size){
    std::ifstream myFile(filename);
    std::string line, vtx, edgs; 
    
    if (!myFile.is_open())
        exit(0);

    unsigned long int n_vertex, n_edges; 

    // get first line of graph file
    myFile >> vtx >> edgs;
    std::getline(myFile, line); // skip this line ? 

    // store the number of vertices and edges of the graph  
    n_vertex = stoul(vtx);
    n_edges = stoul(edgs); 

    unsigned long int my_part = n_vertex / world_size;
	unsigned long int remainder = n_vertex % world_size;

    // revisar 
    no_local_vtx = (remainder <= rank) ? my_part : my_part+1; //: ((rank == world_size-1) ? my_part+remainder : my_part)
    
    local_nodes.reserve(no_local_vtx); 

    unsigned long int total = 0, n_boundary = 0;
	cnt=0;

	// iterate over the rest of the file and create remaining CSR vectors 
    while (std::getline(myFile, line)) {
		std::istringstream iss(line);
		std::string word;
		if(total<vtx_begin){ // use a better way to skip lines ? 
			total++;
            continue;
		}if(total>vtx_end)
			break; 
        
        local_nodes.push_back(new LocalNode); 
        local_nodes.node_weight = 1;    // default = 1 
        local_nodes.current_label = -1;  // unlabeled 
        local_nodes.next_label = -1;
        // local_nodes.edges = ; // maybe store memory for edges 
        local_nodes.is_boundary = false; 

		// iterate over the words on the file 
		//the rest are its neighbours 
		while (iss >> word) {
            local_nodes.edges.push_back(stoi(word));
            // check if current node is boundary (have neighbors outside)
            if( stoi(word) < vtx_begin || stoi(word) > vtx_end ){
                local_nodes.is_boundary = true; 
            
            }

			cnt++; 
		}
        total++;
    }
    myFile.close();
}

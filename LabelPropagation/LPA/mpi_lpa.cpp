#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <vector>
     

#define KARATE_CLUB 0

// gets the label of a node 
long int get_local_id(int global_vtx_id, int rank, std::vector<long int> *vtx_dist){
    // we have a global id [0 to nº of vtx of graph] and want a local id [0 to vtx in this PE]
    return global_vtx_id-(*vtx_dist)[rank]; //(-1 bc it is for graphs starting at 1 ->> CHANGE THIS PLS) 

    // return vtx_dist[rank]+local_vertex; // get global from local 
}

int get_target_PE(int global_vtx_id, int rank, std::vector<long int> *vtx_dist){
    for(int i=1; i<vtx_dist->size(); i++)
        if(global_vtx_id < (*vtx_dist)[i]) return i; 
    return -1; // Not found  
}

// returns the neighbors of a node 
void get_neighbors(int rank, long int global_vtx_id, std::vector<long int> *vtx_dist, 
                std::vector<long int> *xadj, std::vector<long int> *adjncy, std::vector<long int> *neighbors){

    if(global_vtx_id < (*vtx_dist)[rank] || global_vtx_id > (*vtx_dist)[rank+1]) //+1 only for graphs starting at 1 -> CHANGE  
        return; 

    long int local_vtx_id = get_local_id(global_vtx_id, rank, vtx_dist), neigh_begin=(*xadj)[local_vtx_id], neigh_end=(*xadj)[local_vtx_id+1]; 
    long int size = neigh_end - neigh_begin;
    neighbors->reserve(size); // save memory for the neighbors of the vertex
    std::vector<long int>::const_iterator first = adjncy->begin() + neigh_begin;
    std::vector<long int>::const_iterator last = adjncy->begin() + neigh_end;
    std::copy(first, last, std::back_inserter(*neighbors)); // get neighbors
}


// this function loads karate graph manually
void load_graph(std::string filename, int rank, int world_size, std::vector<long int> *vtx_dist, 
                std::vector<long int> *xadj, std::vector<long int> *adjncy, std::vector<long int> *global_ids, std::vector<bool> *is_boundary){
    std::ifstream myFile(filename);
    std::string line, vtx, edgs; 
    
    if (!myFile.is_open())
        exit(0); // terminate if we cannot open the file 

    // The input file is in METIS format 
        // - first line has 3 values: nº of --, nº of --, nº of --
        // - other lines contain an adjacency list for each vertex of the graph   

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
    unsigned long int no_of_local_vtx = (remainder <= rank) ? my_part : my_part+1; //: ((rank == world_size-1) ? my_part+remainder : my_part)

    // construct CSR format vectors   
    // reserve memory
    vtx_dist->reserve(world_size+1); 
    xadj->reserve(no_of_local_vtx+1); 
    global_ids->reserve(no_of_local_vtx); 
    is_boundary->reserve(no_of_local_vtx); 

    //1r elem vtx_dist = 0
    (*vtx_dist).push_back(0);
    unsigned long int cnt = 0; 
    
    for(unsigned long int i=0; i<world_size; i++){
        // set the number of vtx on each process
        cnt += (remainder <= i) ? my_part : my_part+1; //: ((i == world_size-1) ? my_part+remainder : my_part)

        (*vtx_dist).push_back(cnt);
    }
    
    //1r elem xadj = 0
	(*xadj).push_back(0);
	unsigned long int vtx_begin = (*vtx_dist)[rank], vtx_end=(*vtx_dist)[rank+1];
   
	//vecinos de cada uno de los agentes 
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
        
        is_boundary->insert(is_boundary->begin() + n_boundary, false); 
		// iterate over the words on the file 
		//the rest are its neighbours 
		while (iss >> word) {
			(*adjncy).push_back(stoi(word));
            // check if current node is boundary (have neighbors outside)
            if( stoi(word) < vtx_begin || stoi(word) > vtx_end )
                (*is_boundary)[n_boundary] = true;
			cnt++; 
		}

		(*xadj).push_back(cnt);
        (*global_ids).push_back(total); 
		
        n_boundary++; 
        total++;
    }

    myFile.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Init variables */
    long int local_vertices, local_edges, boundary_vertices;
    std::vector<long int> vtx_dist, xadj, adjncy, global_ids;
    std::vector<bool> is_boundary;

    std::string filename = "/home/crispq/parallelCommunityDetection/GraphExamples/0_karate_club_metis.txt";

    /* Load graph from file & Construct graph into CSR format */
    // Each process loads its portion of the vertices of the graph 
    load_graph(filename, rank, world_size, &vtx_dist, &xadj, &adjncy, &global_ids, &is_boundary);

    ///////////////////////////////////////////////////////////
    ////////////////////////// DEBUG //////////////////////////
    ///////////////////////////////////////////////////////////
    
    // tested with a graph starting from 1 -> should work with one starting from 0 tho ? 
    // if(rank==0){
    //     std::cout << "Nº of vertices for each proc: " << std::endl;
    //     for (long int i=0; i<vtx_dist.size(); i++)
    //         std::cout << vtx_dist[i] << ' ';
    //     std::cout << std::endl; 

    //     // std::cout << "Nº of vertices belonging to "<< rank <<" : " << xadj.size() << std::endl;
    //     // for (long int i=0; i<xadj.size(); i++)
    //     //     std::cout << xadj[i] << ' ';
    //     // std::cout << std::endl; 

    //     std::cout << "Edge list of "<< rank <<" size : " << adjncy.size() << std::endl;
    //     for (long int i=0; i<adjncy.size(); i++)
    //         std::cout << adjncy[i] << ' ';
    //     std::cout << std::endl; 

    //     // std::cout << "Global IDs of "<< rank <<" size : " << global_ids.size() << std::endl;
    //     // for (long int i=0; i<global_ids.size(); i++)
    //     //     std::cout << global_ids[i] << ' ';
    //     // std::cout << std::endl; 

    //     std::cout << "Boundary nodes of "<< rank <<" size : " << is_boundary.size() << std::endl;
    //     for (long int i=0; i<is_boundary.size(); i++)
    //         std::cout << is_boundary[i] << ' ';
    //     std::cout << std::endl; 
    // }

    ///////////////////////////////////////////////////////////
    /////////////////////// END DEBUG /////////////////////////
    ///////////////////////////////////////////////////////////
    
    /* TO DO - Set up ghosts */
    // Each Process determines its boundary vertices and the required ghosts 
    // (+including the processes owning copies of boundary nodes)
    

    /* Initialization of the graph labels -> label = the id of the node */
    std::vector<long int> local_labels, local_next_labels, ghost_labels; //= new std::vector<long int>(xadj.size()-1); // size = nº of vtx 
    local_labels.reserve(global_ids.size()); 
    local_next_labels.reserve(global_ids.size()); 

    // IDs to initial labels for local vertices  
    // std::cout << "Copying IDs... ["<< rank << "]" << std::endl; 
    std::copy(global_ids.begin(), global_ids.end(), std::back_inserter(local_labels));
    // std::cout << "COPY - DONE ["<< rank << "]" << std::endl; 
    
    ///////////////////////////////////////////////////////////
    ////////////////////////// DEBUG //////////////////////////
    ///////////////////////////////////////////////////////////
    // if(rank==1){
    //     std::cout << "Global IDs of "<< rank <<" size : " << global_ids.size() << std::endl;
    //     for (long int i=0; i<global_ids.size(); i++)
    //         std::cout << global_ids[i] << ' ';
    //     std::cout << std::endl; 

    //     // std::cout << "Labels of "<< rank <<" size : " << labels.size() << std::endl;
    //     // for (long int i=0; i<labels.size(); i++)
    //     //     std::cout << labels[i] << ' ';
    //     // std::cout << std::endl; 
    // }
    ///////////////////////////////////////////////////////////
    /////////////////////// END DEBUG /////////////////////////
    ///////////////////////////////////////////////////////////

    /* TO DO - Init communications -> get which boundary nodes are connected to each PE */
        // maybe I dont need to ? 
            // when a node changes label and is boundary -> for each of its neighbors attach the node to mssg buffer they belong to 
            // when sending message -> check if boundary and belongs to each PE -> send to whomever needs them

    
    /* TO DO - Each process runs LPA until convergence */
        /* TO DO - Synchronization */
        // Boundaries / Ghost must be sent / rcv 
    bool convergence = false;
    while(!convergence){
        // Run LPA for each local node 
        std::vector<long int> label_count;
        std::vector<bool> is_maximal; 

        label_count.resize(); 

        // for each local node 
        for (long int i = 0 ; i < xadj.size() ; i++ ){
            

            // for each of the neighbors of a node
            for(long int j = 0 ; j < adjncy.size() ; j++ ){
                // look for the maximum label(s) among its neighbors
                // get neighbor label 
                if(){

                }
                
                if(label_count[]){

                }                
            
            }
        }

        // send data to next step 

        convergence=true; 
    }
    
   /* TO DO - Calculate Modularity */


    MPI_Finalize();
}

// TESTING NEIGHBOR FUNCTION 
    // std::vector<long int> neighbors; 
    // long int global_id = 12;
    // std::cout << "Getting neighbors... ["<< rank << "]" << std::endl; 
    // get_neighbors(rank, global_id, &vtx_dist, &xadj, &adjncy, &neighbors);
    // std::cout << "Getting neighbors - DONE ["<< rank << "]" << std::endl; 

    // // ///////////////////////////////////////////////////////////
    // // ////////////////////////// DEBUG //////////////////////////
    // // ///////////////////////////////////////////////////////////
    // // if(rank==0){
    //     std::cout << "Neighbors of "<< global_id << " at rank : " << rank << " [total : " << neighbors.size() << "]" << std::endl;
    //     for (long int i=0; i<neighbors.size(); i++)
    //         std::cout << neighbors[i] << ' ';
    //     std::cout << std::endl; 
    // }
    // ///////////////////////////////////////////////////////////
    // /////////////////////// END DEBUG /////////////////////////
    // ///////////////////////////////////////////////////////////
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <vector>

// #include "define.h"
#include "comm.h"
#include "graph_structure.h"

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

    // if(rank==0)
    //     std::cout << "Here I am" << std::endl;

    /* Init variables */
    long int local_vertices, local_edges, boundary_vertices;
    std::vector<long int> vtx_dist, xadj, adjncy, global_ids;
    std::vector<bool> is_boundary;

    // std::string filename = "/home/crispq/label_propagation/ParallelCommunityDetection/GraphExamples/0_karate_club_metis.txt";
    // std::string filename = "/home/crispq/parallelCommunityDetection/GraphExamples/small_test.txt";
    std::string filename = "/home/crispq/parallelCommunityDetection/GraphExamples/0_karate_club_metis.txt";

    /* Load graph from file & Construct graph into CSR format */
    // Each process loads its portion of the vertices of the graph 
    // load_graph(filename, rank, world_size, &vtx_dist, &xadj, &adjncy, &global_ids, &is_boundary);
    
    // if(rank==0)
    //     std::cout << "1st load done" << std::endl;

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

    /* TO DO - Test Graph Structure */
    // create graph obj 
    DistributedGraph graph; 
    CommunicationHandler cm; 

    // if(rank==0){
    //     std::cout << "Loading a graph..." << std::endl; 
    //     }
    // std::cout << "Loading the graph" << std::endl;
    graph.create_graph_from_METIS(filename);
    // std::cout << "LoadED !!" << std::endl;

    cm.init_communications(graph.ghost_nodes); 


    // std::cout << "Nº of local nodes in rank "<< rank << " - " << graph.local_nodes->size() << std::endl;
    // if(rank==0){
    //     std::cout << "Loading a graph - DONE !" << std::endl; 
    //     }

    /* TO DO - Init communications -> get which boundary nodes are connected to each PE */
        // maybe I dont need to ? 
            // when a node changes label and is boundary -> for each of its neighbors attach the node to mssg buffer they belong to 
            // when sending message -> check if boundary and belongs to each PE -> send to whomever needs them


    // 

    /* TO DO - RUN LPA */
    /* DEBUG __ Checking local labels */
    // int check_r = 2;
    // if(rank == check_r){
    //     // std::cout << "Ghost on rank " << rank << " has id " << 2 << " and global id is " << graph.ghost_global_ids[2] << " local id of the ghost : " << (*graph.ghost_nodes)[graph.from_ghost_global_to_index(3)].id << std::endl; 
    //     // std::cout << "Ghost on rank " << rank << " has id " << 5 << " and global id is " << graph.ghost_global_ids[5] << " local id of the ghost : " << (*graph.ghost_nodes)[graph.from_ghost_global_to_index(5)].id << std::endl; 
    //     // std::cout << "Ghost on rank " << rank << " has id " << 8 << " and global id is " << graph.from_ghost_global_to_index(8) << " local id of the ghost : " << (*graph.ghost_nodes)[graph.from_ghost_global_to_index(8)].id << std::endl; 


    //     // int cnt_b = 0;

    //     std::cout << "BEFORE Local vtx on " << rank << " [id:label]: "; 
    //     for(auto local_vtx : (*graph.local_nodes)){
    //         // if(local_vtx.is_boundary){
    //         //     std::cout << "vtx w local id : " << local_vtx.id << " [global:" << graph.from_local_to_global(local_vtx.id) << "]" << " is boundary" << std::endl;
    //         //     cnt_b++;
    //         // }
    //         std::cout << local_vtx.id << " - " << local_vtx.current_label << std::endl; 
    //     }

    //     std::cout << "BEFORE Ghost vtx on " << rank << " [id:label]: "; 
    //     for(auto local_vtx : (*graph.ghost_nodes)){
    //         // if(local_vtx.is_boundary){
    //         //     std::cout << "vtx w local id : " << local_vtx.id << " [global:" << graph.from_local_to_global(local_vtx.id) << "]" << " is boundary" << std::endl;
    //         //     cnt_b++;
    //         // }
    //         std::cout << local_vtx.id << " - " << local_vtx.current_label << std::endl; 
    //     }
    //     // std::cout << "nº boundary " << cnt_b << std::endl;
    // }

    // run this until convergence 
    // for all local vtx : 
    int nsteps = 0; 
        while(nsteps < 20){
        for(auto local_vtx : (*graph.local_nodes)){
            // init label counter & maximal labels 
            std::unordered_map<T,T> label_cnt;
            std::vector<T> max_labels; 
            LABEL_T max_label_value = 0;

            // if(rank == 1)
            //     if(local_vtx.id == 1) // vtx 8 
            //         std::cout << "I'm vtx 4 at step " << nsteps << " printing the label of my neighbors "; 

            // for all neighbors of the local vtx 
            for(auto neigh : (*local_vtx.edges)){
                T n_idx = neigh.target; 
                LABEL_T n_label = -1; 
                
                if(n_idx < graph.get_local_vtx()){
                    n_label = (*graph.local_nodes)[n_idx].current_label; 
                }else{
                    int local_idx =  graph.from_local_ghost_to_index(n_idx);
                    n_label = (*graph.ghost_nodes)[local_idx].current_label; 
                }

                if(rank == 1)
                if(local_vtx.id == 1) 
                        std::cout << n_label << " "; 

                /* add value to label counter */
                // if label_cnt does not have this label yet -> add it 
                if(label_cnt.find(n_label) == label_cnt.end()){
                    label_cnt[n_label] = 1;
                }else{
                    label_cnt[n_label] += 1; 
                }

                // if new label = maximal 
                if(max_label_value < label_cnt[n_label]){
                    max_label_value = label_cnt[n_label]; 
                    max_labels.clear(); // resize to 0 ? 
                    max_labels.push_back(n_label); // append new max label 

                // if new label is the same value another one found previously
                }else if (label_cnt[n_label] == max_label_value){
                    max_labels.push_back(n_label);
                }
            }

            // if there is at least one max label -> change label 
            if(max_labels.size() >= 1){
                // if(rank == 1)
                //     if(local_vtx.id == 1) 
                //     {  
                //         std::cout << std::endl; 
                //         std::cout << "Found " << max_labels.size() << " maximal labels for this vtx : ";
                //         for( int i = 0; i < max_labels.size() ; i++){
                //                 std::cout << max_labels[i] << " ";  
                //         }
                //         std::cout << std::endl;
                //     }

                int rng_label = std::rand() % max_labels.size();   // pick randomly one label

                if(max_labels[rng_label] != local_vtx.current_label){
                    graph.set_next_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                    // append changed vtx to queue to be sent ? 
                    cm.addToSend(&graph, local_vtx); 

                    // if(rank == 1)
                    //     if(local_vtx.id == 1) 
                    //         std::cout << "Adding next label as " << max_labels[rng_label];
                }
            }
        }
        // if(rank == 2)
        //     std::cout << std::endl; 
        cm.send_recv_data(&graph); 
        graph.update_labels(); 

        nsteps++; //number of LPA steps, change later 
    }

    std::cout << "DONE WITH LPA" << std::endl; 
    // communications are to be done here -> send any modified boundary vtx, 
                                          // rcv any ghost sent by other PEs 
    /* HOW TO ? */
    // each PE = slice of the cake -> write in order in file 
    std::string output_filename = "/home/crispq/parallelCommunityDetection/output_small_test.txt";
    std::ofstream myOutputFile;

    int current_writer = 0; 
    bool done = false; 
    while( !done ){
        // std::cout << "On rank " << rank << "the current_writer is " << current_writer << std::endl;
        if(rank == current_writer){
            // open the file 
            // std::cout << "Opening the file on " << rank << std::endl;
            if(rank == 0)
                myOutputFile.open(output_filename);
            else
                myOutputFile.open(output_filename, std::ios::app);

            // std::cout << "Writting file on " << rank << std::endl;
            // write into the output file its portion of the graph + labels 
            for( auto local_vtx : (*graph.local_nodes) ){
                T global_id = graph.from_local_to_global(local_vtx.id); 
                myOutputFile << global_id << " " << local_vtx.current_label << std::endl; 
            }

            //close file 
            // std::cout << "Closing the file on " << rank << std::endl;
            myOutputFile.close(); 

            // send the next writter a message 
            current_writer++;
            // std::cout << "Sending the msg from " << rank << " to " << rank+1 << " with current_writer " << current_writer << std::endl;
            if(rank < world_size-1)
                MPI_Send(&current_writer, 1, MPI_INT, rank+1, 30, MPI_COMM_WORLD);

            // this one is done 
            done = true;  
        }else{
            MPI_Recv(&current_writer, 1, MPI_INT, rank-1, 30, MPI_COMM_WORLD, MPI_STATUS_IGNORE);     // wait til getting message 
        }
    }

    // std::cout << "Done with my writting from " << rank << std::endl;


    MPI_Barrier(MPI_COMM_WORLD);

    // calculate modularity ? 

    // if(rank == check_r){
    //     // std::cout << "Ghost on rank " << rank << " has id " << 2 << " and global id is " << graph.ghost_global_ids[2] << " local id of the ghost : " << (*graph.ghost_nodes)[graph.from_ghost_global_to_index(3)].id << std::endl; 
    //     // std::cout << "Ghost on rank " << rank << " has id " << 5 << " and global id is " << graph.ghost_global_ids[5] << " local id of the ghost : " << (*graph.ghost_nodes)[graph.from_ghost_global_to_index(5)].id << std::endl; 
    //     // std::cout << "Ghost on rank " << rank << " has id " << 8 << " and global id is " << graph.from_ghost_global_to_index(8) << " local id of the ghost : " << (*graph.ghost_nodes)[graph.from_ghost_global_to_index(8)].id << std::endl; 


    //     // int cnt_b = 0;

    //     std::cout << "AFTER Local vtx on " << rank << " [id:label]: "; 
    //     for(auto local_vtx : (*graph.local_nodes)){
    //         // if(local_vtx.is_boundary){
    //         //     std::cout << "vtx w local id : " << local_vtx.id << " [global:" << graph.from_local_to_global(local_vtx.id) << "]" << " is boundary" << std::endl;
    //         //     cnt_b++;
    //         // }
    //         std::cout << local_vtx.id << " - " << local_vtx.current_label << std::endl; 
    //     }

    //     std::cout << "AFTER Ghost vtx on " << rank << " [id:label]: "; 
    //     for(auto local_vtx : (*graph.ghost_nodes)){
    //         // if(local_vtx.is_boundary){
    //         //     std::cout << "vtx w local id : " << local_vtx.id << " [global:" << graph.from_local_to_global(local_vtx.id) << "]" << " is boundary" << std::endl;
    //         //     cnt_b++;
    //         // }
    //         std::cout << local_vtx.id << " - " << local_vtx.current_label << std::endl; 
    //     }
    //     // std::cout << "nº boundary " << cnt_b << std::endl;
    // }

    MPI_Barrier(MPI_COMM_WORLD); 
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
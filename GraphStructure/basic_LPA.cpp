#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <vector>

#include "comm.h"
#include "graph_structure.h"
#include <unordered_set>

void run_LPA(DistributedGraph *graph, CommunicationHandler *cm){
    int nsteps = 0; 
    while(nsteps < 20){
        for(auto local_vtx : (*graph->local_nodes)){
            // init label counter & maximal labels 
            std::unordered_map<T,T> label_cnt;
            std::vector<T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;

            // for all neighbors of the local vtx 
            for(auto neigh : (*local_vtx.edges)){
                T n_idx = neigh.target; 
                LABEL_T n_label = -1; 
                
                if(n_idx < graph->get_local_vtx()){
                    n_label = (*graph->local_nodes)[n_idx].current_label; 
                }else{
                    int local_idx =  graph->from_local_ghost_to_index(n_idx);
                    n_label = (*graph->ghost_nodes)[local_idx].current_label; 
                }

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
                int rng_label = std::rand() % max_labels.size();   // pick randomly one label

                if(max_labels[rng_label] != local_vtx.current_label){
                    graph->set_next_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                    // append changed vtx to queue to be sent ? 
                    cm->addToSend(graph, local_vtx); 
                }
            }
        }
        cm->send_recv_data(graph); 
        graph->update_local_labels(); 

        nsteps++; //number of LPA steps, change later 
    }
}

void write_output(std::string output_filename, DistributedGraph * graph, int rank, int world_size){
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

            // write into the output file its portion of the graph + labels 
            for( auto local_vtx : (*graph->local_nodes) ){
                T global_id = graph->from_local_to_global(local_vtx.id); 
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
    myOutputFile.close(); 
}

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::string filename, output_filename = "../output_small_test.txt";

    if(argc < 2){
        std::cout << "Please pass the file path as an argument to the program." << std::endl;
        return 1;
    }else{ 
        filename = argv[1];
    }

    /* Create graph and communication handler */
    DistributedGraph graph; 
    CommunicationHandler cm; 

    /* Init graph from input file */
    graph.create_graph_from_METIS(filename);
    
    /* Init communication structure */
    cm.init_communications(graph.ghost_nodes); 

    // RUN LPA 
    run_LPA(&graph, &cm);

    /* Writes the output of LPA to a file */
    write_output(output_filename, &graph, rank, world_size);
    
    /////////////////////////////////////
    /////* Calculate modularity ?? */////
    /////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD); 
    MPI_Finalize();
}
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <vector>

#include "CommunicationHandler.h"
#include "DistributedGraph.h"
#include <unordered_set>



void run_LPA(DistributedGraph *graph, CommunicationHandler *cm){
    int nsteps = 0; 
    int end_condition = 0;

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    std::vector<ID_T>end_c((*graph->get_local_vertices()).size());
    while(end_condition == 0){
        end_c.clear(); 

        for(auto local_vtx : (*graph->get_local_vertices())){
            if(local_vtx.is_boundary) // for non boundary vtx 
                continue; 
            // init label counter & maximal labels 
            std::unordered_map<ID_T,ID_T> label_cnt;
            std::vector<ID_T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;
            
            // label_cnt[local_vtx.current_label] = 1; 

            // for all neighbors of the local vtx 
            for(auto neigh : (*local_vtx.edges)){
                ID_T n_idx = neigh.target; 
                LABEL_T n_label = -1; 
                
                n_label = (*graph->get_local_vertices())[n_idx].current_label; 

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
                }
            }

            // Evaluate if LPA it finished 
            // nº of thimes current max label = max_label count -> we can stop 
            if(label_cnt[local_vtx.current_label] == max_label_value){
                end_c.push_back(1);
            }else{
                end_c.push_back(0);
            }
        }

        if(nsteps > 0){
            cm->recv_data();
            graph->update_ghost_labels(cm->get_recv_buffer()); 
            cm->wait_requests(); 
        }

        for(auto local_vtx : (*graph->get_local_vertices())){
            if(!local_vtx.is_boundary) // for boundary vtx 
                continue; 

            std::unordered_map<ID_T,ID_T> label_cnt;
            std::vector<ID_T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;

            // if(rank == 1 && nsteps == 1)
            //     std::cout << "I'm vtx " << local_vtx.id << " with neighbors : "; 

            // for all neighbors of the local vtx 
            for(auto neigh : (*local_vtx.edges)){
                ID_T n_idx = neigh.target; 
                LABEL_T n_label = -1; 
                
                if(n_idx < graph->get_local_vtx()){
                    n_label = (*graph->get_local_vertices())[n_idx].current_label; 
                }else{
                    int local_idx =  graph->from_local_ghost_to_index(n_idx);
                    n_label = (*graph->get_ghost_vertices())[local_idx].current_label; 
                    int ghost_pe_id = graph->get_ghost_vertices()->at(local_idx).pe_id;

                    if(boundary_neighbor_PEs.find(ghost_pe_id) == boundary_neighbor_PEs.end())
                        boundary_neighbor_PEs.insert(ghost_pe_id); 
                }

                // if(rank == 1 && nsteps == 1)
                //     std::cout << "[" << n_idx << ", " << n_label << "]" << " ";

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

            // if(rank == 1 && nsteps == 1)
            //     std::cout <<  std::endl; 

            // if there is at least one max label -> change label 
            if(max_labels.size() >= 1){
                int rng_label = std::rand() % max_labels.size();   // pick randomly one label
                if(max_labels[rng_label] != local_vtx.current_label){
                    graph->set_next_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                    // if(rank == 2 && nsteps == 0)
                    //     std::cout << "I'm vtx " << local_vtx.id << " updating label to : " << max_labels[rng_label] << std::endl; 
                    cm->add_all_to_send(&boundary_neighbor_PEs, graph->from_local_to_global(local_vtx.id), max_labels[rng_label]); 
                }
            }

            // Evaluate if LPA it finished 
            // nº of thimes current max label = max_label count -> we can stop 
            if(label_cnt[local_vtx.current_label] == max_label_value){
                end_c.push_back(1);
            }else{
                end_c.push_back(0);
            }
        }

        
        /* TO DO */
        // Evaluate if LPA it finished 
        if(nsteps > 0){
            // evaluate end condition for the local vtx on this PE  
            end_condition = ( graph->get_local_vtx() == std::accumulate(end_c.begin(),end_c.end(),0)) ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &end_condition, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

            if(end_condition == world_size){
                end_condition = 1; 
            }else 
                end_condition = 0; 
        }

        if(end_condition == 0){
            // cm->send_recv_data(graph); 
            cm->send_data();
            graph->update_local_labels(); 
            nsteps++; //number of LPA steps, change later 
        }
    }

    if(rank == 0)
        std::cout << " LPA ended with " << nsteps << " steps" << std::endl;
    // Need something to handle degree 0 vtx <- 
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
            for( auto local_vtx : (*graph->get_local_vertices()) ){
                ID_T global_id = graph->from_local_to_global(local_vtx.id); 
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
    cm.init_communications(graph.get_ghost_vertices()); 

    // cm.order_ghosts(graph.get_ghost_vertices());

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
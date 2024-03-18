#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <vector>

#include "../../GraphStructure/CommunicationHandler.h"
#include "../../GraphStructure/DistributedGraph.h"
#include "../../GraphStructure/DistributedCompressedGraph.h"
#include <unordered_set>


void run_with_compression_LPA(DistributedCompressedGraph *graph, CommunicationHandler *cm, int ns){
    int nsteps = 0; 
    int end_condition = 0;

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* TO DO - Send / recv boundary that are not active */
    cm->send_rcv_inactive(graph); 
    graph->update_inactive_ghosts(cm->get_recv_buffer()); 

    // // std::cout << "Getting the following ghosts on rank " << rank << " :"; 
    // int inactive_cnt = 0, bndry_cnt = 0, b_inactive=0, inctv_gh=0; 
    // for( auto gh : (*graph->get_ghost_vertices()) ){
    //     if(!gh.active)
    //         inctv_gh++;
    //     // std::cout << " [state : " << gh.active << ", id :  " << gh.id << "] " << std::endl; 
    // }
    // std::cout << "inactive ghosts : " << inctv_gh << " from " << (*graph->get_ghost_vertices()).size() << std::endl;
    
    // for( auto gh : (*graph->get_local_vertices()) ){
    //     if(!gh.active)
    //         inactive_cnt++;
    //     if(gh.is_boundary && !gh.active)
    //         b_inactive++;
    //     if (gh.is_boundary)
    //         bndry_cnt++;
    //     // std::cout << " [state : " << gh.active << ", id :  " << gh.id << ", label: " << gh.current_label << ",] " << std::endl; 
    // }
    // std::cout << "total inactive vtx : " <<  inactive_cnt << " bndry inactive: " << b_inactive << " and bndry : " << bndry_cnt << " from : " << (*graph->get_local_vertices()).size() << std::endl;
    // std::cout << std::endl;

    // int cnt_deg_0 = 0, cnt_deg_1 = 0, cnt_deg_2 = 0; 
    // for (int i = 0; i < (*graph->get_local_vertices()).size(); i++){
    //     if( (*graph->get_degrees())[i] == 0 ){
    //         cnt_deg_0++;
    //     }else if ( (*graph->get_degrees())[i] == 1 ){
    //         cnt_deg_1++;
    //     }else if ( (*graph->get_degrees())[i] == 2 )
    //         cnt_deg_2++;
    // }

    // std::cout << "We have a total of vtx of degree 0 " << cnt_deg_0 << " degree 1 : " << cnt_deg_1 << " and degree 2 : " << cnt_deg_2 << std::endl; 

    
    std::vector<ID_T>end_c((*graph->get_local_vertices()).size());
    while(end_condition == 0 && nsteps < ns){
        end_c.clear(); 
        for(auto local_vtx : (*graph->get_local_vertices())){
            if(local_vtx.is_boundary){ // for non boundary vtx and active 
                continue; 
            }else if(!local_vtx.active){
                end_c.push_back(1); 
                continue; 
            }

            // init label counter & maximal labels 
            std::unordered_map<ID_T,ID_T> label_cnt;
            std::vector<ID_T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;
            
            label_cnt[local_vtx.current_label] = 1; 

            // for all neighbors of the local vtx 
            for(auto neigh : (*local_vtx.edges)){
                ID_T n_idx = neigh.target; 
                LABEL_T n_label = -1; 
               
                if(!graph->get_local_vertices()->at(n_idx).active)
                    continue; 
                
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
            if(label_cnt[local_vtx.current_label] >= max_label_value){
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
            if(!local_vtx.is_boundary){ // for boundary vtx  
                continue; 
            }else if(!local_vtx.active){
                end_c.push_back(1);
                continue; 
            }

            std::unordered_map<ID_T,ID_T> label_cnt;
            std::vector<ID_T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;

            label_cnt[local_vtx.current_label] = 1; 

            // if(rank == 1 && nsteps == 1)
            //     std::cout << "I'm vtx " << local_vtx.id << " with neighbors : "; 

            // for all neighbors of the local vtx 
            for(auto neigh : (*local_vtx.edges)){
                                ID_T n_idx = neigh.target; 
                LABEL_T n_label = -1; 
                
                if(n_idx < graph->get_local_vtx()){
                    n_label = (*graph->get_local_vertices())[n_idx].current_label; 

                    if(!graph->get_local_vertices()->at(n_idx).active)
                        continue; 
                }else{
                    int local_idx =  graph->from_local_ghost_to_index(n_idx);
                    n_label = (*graph->get_ghost_vertices())[local_idx].current_label; 
                    int ghost_pe_id = graph->get_ghost_vertices()->at(local_idx).pe_id;

                    if(!graph->get_ghost_vertices()->at(local_idx).active)
                        continue; 

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
                    // if(rank == 2)
                    //     std::cout << "I'm vtx " << local_vtx.id << " updating label to : " << max_labels[rng_label] << " at " << rank << std::endl; 
                    // cm->add_all_to_send(&boundary_neighbor_PEs, graph->from_local_to_global(local_vtx.id), max_labels[rng_label]); 
                }
            }

            // Evaluate if LPA it finished 
            // nº of thimes current max label = max_label count -> we can stop 
            if(label_cnt[local_vtx.current_label] >= max_label_value){
                end_c.push_back(1);
            }else{
                end_c.push_back(0);
                // if(nsteps > 10 ){
                //     if( rank == 0 ){
                //         std::cout << "[ " << local_vtx.current_label << " -- " << label_cnt[local_vtx.current_label] << "/" << max_label_value << "] VTX " << local_vtx.id << " has neighbors with labels : "; 
                //         for(auto e : (*local_vtx.edges)){
                //             ID_T id = e.target; 
                //             int lab; 
                //             if(id < graph->get_local_vtx()){
                //                 lab = (*graph->get_local_vertices())[id].current_label; 
                //             }else{
                //                 int local_idx =  graph->from_local_ghost_to_index(id);
                //                 lab = (*graph->get_ghost_vertices())[local_idx].current_label;
                //             }
                //             std::cout << "[" << id << ", " << lab << "] ";
                //         }
                //         std::cout << std::endl; 
                //     }
                // }
            }       
        }

        
        /* TO DO */
        // Evaluate if LPA it finished 
        if(nsteps > 0){
            // evaluate end condition for the local vtx on this PE  
            end_condition = ( graph->get_local_vtx() == std::accumulate(end_c.begin(),end_c.end(),0)) ? 1 : 0;
            MPI_Allreduce(MPI_IN_PLACE, &end_condition, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 

            // if(nsteps > 10 ){
            //     // if( rank == 0 )
            //     std::cout << "At rank " << rank << " I have " << end_condition << " for step " << nsteps << " : ";
            //     for(auto e : end_c){
            //         std::cout << e << " ";
            //     }
            //     std::cout << std::endl; 
            // }

            if(end_condition == world_size){
                end_condition = 1; 
            }else 
                end_condition = 0; 
        }

        if(end_condition == 0 && nsteps < ns-1){
            // cm->send_recv_data(graph); 
            cm->send_data();
            graph->update_local_labels(); 
        }
        nsteps++; //number of LPA steps, change later 
    }

    if(rank == 0)
        std::cout << " LPA ended with " << nsteps << " steps" << std::endl;
    // Need something to handle degree 0 vtx <- 
}


void run_LPA(DistributedGraph *graph, CommunicationHandler *cm, int ns){
    int nsteps = 0; 
    int end_condition = 0;

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int updated = graph->get_total_vtx();

    int updt_thresh = ceil(graph->get_total_vtx()*0.0001); 

    std::cout << "Thresh : " << updt_thresh << " nº vtx : " <<  updated << std::endl; 
    std::vector<ID_T>end_c((*graph->get_local_vertices()).size());
    while(updated > updt_thresh && nsteps < ns){
        // end_c.clear(); 
        updated = 0; 
        for(auto local_vtx : (*graph->get_local_vertices())){
            if(local_vtx.is_boundary) // for non boundary vtx 
                continue; 
            // init label counter & maximal labels 
            std::unordered_map<ID_T,ID_T> label_cnt;
            std::vector<ID_T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;
            
            label_cnt[local_vtx.current_label] = 1; 

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
                if(find(max_labels.begin(), max_labels.end(), local_vtx.current_label) != max_labels.end())
                    continue; 
                int rng_label = std::rand() % max_labels.size();   // pick randomly one label
                if(max_labels[rng_label] != local_vtx.current_label){
                    graph->set_next_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                    updated++;
                }
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

            label_cnt[local_vtx.current_label] = 1; 

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


            // if there is at least one max label -> change label 
            if(max_labels.size() >= 1){
                if(find(max_labels.begin(), max_labels.end(), local_vtx.current_label) != max_labels.end())
                    continue; 
                int rng_label = std::rand() % max_labels.size();   // pick randomly one label
                if(max_labels[rng_label] != local_vtx.current_label){
                    graph->set_next_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                    cm->add_all_to_send(&boundary_neighbor_PEs, graph->from_local_to_global(local_vtx.id), max_labels[rng_label]);
                    updated++; 
                }
            }

        }

        
        /* TO DO */
        // Evaluate if LPA it finished 
        // if(nsteps > 0){
            // evaluate end condition for the local vtx on this PE  
            // end_condition = ( graph->get_local_vtx() == std::accumulate(end_c.begin(),end_c.end(),0)) ? 1 : 0;
            // MPI_Allreduce(MPI_IN_PLACE, &end_condition, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
            MPI_Allreduce(MPI_IN_PLACE, &updated, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
            std::cout << "At step " << nsteps << " nº vtx updated : " <<  updated << std::endl; 

            // if(nsteps > 10 ){
            //     // if( rank == 0 )
            //     std::cout << "At rank " << rank << " I have " << end_condition << " for step " << nsteps << " : ";
            //     for(auto e : end_c){
            //         std::cout << e << " ";
            //     }
            //     std::cout << std::endl; 
            // }

            // if(end_condition == world_size){
            //     end_condition = 1; 
            // }else 
            //     end_condition = 0; 
        // }

        if(updated > updt_thresh && nsteps < ns-1){
            // cm->send_recv_data(graph); 
            cm->send_data();
            graph->update_local_labels(); 
        }
        nsteps++; //number of LPA steps, change later 
        if(rank == 0)
            std::cout << " Current step : " << nsteps << std::endl;
    }

    if(rank == 0)
        std::cout << " LPA ended with " << nsteps << " steps" << std::endl;
    // Need something to handle degree 0 vtx <- 
}

void run_asynchronous_LPA(DistributedGraph *graph, CommunicationHandler *cm, int ns){
    int nsteps = 0; 
    int end_condition = 1;

    int world_size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double updt_thresh = graph->get_total_vtx()*0.95; 

    // int running = 0; 
    int eval = false; 
    std::vector<ID_T>end_c((*graph->get_local_vertices()).size());
    while(end_condition > 0 && nsteps < ns){
        end_c.clear(); 
        if(eval == true)
            end_condition = 0;

        for(auto local_vtx : (*graph->get_local_vertices())){
            if(local_vtx.is_boundary) // for non boundary vtx 
                continue; 
            // init label counter & maximal labels 
            std::unordered_map<ID_T,ID_T> label_cnt;
            std::vector<ID_T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;
            
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
                if(eval == false){
                    int rng_label = std::rand() % max_labels.size();   // pick randomly one label
                    // if(max_labels[rng_label] != local_vtx.current_label){
                    graph->set_next_current_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                    // }
                }else{
                    if(label_cnt[(*graph->get_local_vertices())[local_vtx.id].current_label] < max_label_value){
                        end_condition = 1;
                    }
                }
            }

            // Evaluate if LPA it finished 
            // nº of thimes current max label = max_label count -> we can stop 
            // if(label_cnt[local_vtx.current_label] >= max_label_value){
                
            //     end_c.push_back(1);
            // }else{
            //     end_c.push_back(0);
            // }
        }

        if(nsteps > 0 && eval == false){
            cm->recv_data();
            graph->update_ghost_labels(cm->get_recv_buffer()); 
            cm->wait_requests(); 
            MPI_Barrier(MPI_COMM_WORLD);
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

                    if(eval == false)
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
                if(eval == false){
                    int rng_label = std::rand() % max_labels.size();   // pick randomly one label
                    if(max_labels[rng_label] != local_vtx.current_label){
                        graph->set_next_current_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                        // if(rank == 2 && nsteps == 0)
                        //     std::cout << "I'm vtx " << local_vtx.id << " updating label to : " << max_labels[rng_label] << std::endl; 
                        cm->add_all_to_send(&boundary_neighbor_PEs, graph->from_local_to_global(local_vtx.id), max_labels[rng_label]); 
                    }
                }else{
                    if(label_cnt[(*graph->get_local_vertices())[local_vtx.id].current_label] < max_label_value){
                        end_condition = 1;
                    }
                }
            }

            // Evaluate if LPA it finished 
            // nº of thimes current max label = max_label count -> we can stop 
            // if(label_cnt[local_vtx.current_label] >= max_label_value){
            // if(label_cnt[(*graph->get_local_vertices())[local_vtx.id].current_label] >= max_label_value){
            //     end_c.push_back(1);
            // }else{
            //     end_c.push_back(0);
            //     // if(nsteps > 10 ){
            //     //     if( rank == 0 ){
            //     //         std::cout << "[ " << local_vtx.current_label << " -- " << label_cnt[local_vtx.current_label] << "/" << max_label_value << "] VTX " << local_vtx.id << " has neighbors with labels : "; 
            //     //         for(auto e : (*local_vtx.edges)){
            //     //             ID_T id = e.target; 
            //     //             int lab; 
            //     //             if(id < graph->get_local_vtx()){
            //     //                 lab = (*graph->get_local_vertices())[id].current_label; 
            //     //             }else{
            //     //                 int local_idx =  graph->from_local_ghost_to_index(id);
            //     //                 lab = (*graph->get_ghost_vertices())[local_idx].current_label;
            //     //             }
            //     //             std::cout << "[" << id << ", " << lab << "] ";
            //     //         }
            //     //         std::cout << std::endl; 
            //     //     }
            //     // }
            // }       
        }

        
        /* TO DO */
        // Evaluate if LPA it finished 
        // if(nsteps > 0){
            // evaluate end condition for the local vtx on this PE  
            // end_condition = ( graph->get_local_vtx() == std::accumulate(end_c.begin(),end_c.end(),0)) ? 1 : 0;
        if(eval == true){
            std::cout << "EC [" << rank << "] : " << end_condition << std::endl; 
            MPI_Allreduce(MPI_IN_PLACE, &end_condition, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
        }

            // if(nsteps > 10 ){
            //     // if( rank == 0 )
            //     std::cout << "At rank " << rank << " I have " << end_condition << " for step " << nsteps << " : ";
            //     for(auto e : end_c){
            //         std::cout << e << " ";
            //     }
            //     std::cout << std::endl; 
            // }

            // if(end_condition == world_size){
            //     end_condition = 1; 
            // }else 
            //     end_condition = 0; 
        // }

        if(end_condition > 0 && nsteps < ns-1 && eval == false){
            cm->send_recv_data(); 
            // cm->send_data();
            // graph->update_local_labels(); 
            nsteps++; //number of LPA steps, change later 
            if(rank == 0)
                std::cout << " Current step : " << nsteps << std::endl;
        }
        eval = !eval; 
    }

    if(rank == 0)
        std::cout << " LPA ended with " << nsteps << " steps" << std::endl;
    // Need something to handle degree 0 vtx <- 
}



void run_basic_asynchronous_LPA(DistributedGraph *graph, CommunicationHandler *cm, int ns){
    int nsteps = 0; 
    bool eval = false; 
    int end_condition = 1; 
    int updated = graph->get_total_vtx(); 

    int updt_thresh = graph->get_total_vtx()*0.00001;

    std::cout << "Threshold : " << updt_thresh << " total : " << graph->get_total_vtx() << "["<< (*graph->get_local_vertices()).size() <<  "]" << std::endl;
    while(updated > updt_thresh && nsteps < ns){
        // updated=110;
        // if(eval)
        //     end_condition = 0;

        for(auto local_vtx : (*graph->get_local_vertices())){
            // init label counter & maximal labels 
            std::unordered_map<ID_T,ID_T> label_cnt;
            std::vector<ID_T> max_labels; 
            LABEL_T max_label_value = 0;
            std::unordered_set<int> boundary_neighbor_PEs;

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

                    // if(!eval)
                    if(boundary_neighbor_PEs.find(ghost_pe_id) == boundary_neighbor_PEs.end())
                        boundary_neighbor_PEs.insert(ghost_pe_id); 
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
                }
                // else if (label_cnt[n_label] == max_label_value){
                //     max_labels.push_back(n_label);
                // }
            }

            // if there is at least one max label -> change label 
            if(max_labels.size() >= 1){
                // if(!eval){
                    int rng_label = std::rand() % max_labels.size();   // pick randomly one label

                    if(max_labels[rng_label] != local_vtx.current_label){
                        graph->set_next_current_label(max_labels[rng_label], local_vtx.id);  // assign new label to local vtx 
                        // updated++; 
                        // append changed vtx to queue to be sent ? 
                        if(local_vtx.is_boundary)
                            cm->add_all_to_send(&boundary_neighbor_PEs, local_vtx.id, max_labels[rng_label]);  
                    } 

            }
        }

        // if(eval == true){
            MPI_Allreduce(MPI_IN_PLACE, &updated, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); 
            
            std::cout << " Updated on step " << nsteps << "   : " << updated << "" << std::endl; 
        // }else{
            // printf("Sending,,,\n");
            cm->send_recv_data(); 
            graph->update_ghost_labels(cm->get_recv_buffer()); 
            cm->wait_requests(); 
            std::cout << "Step : " << nsteps << std::endl;
            nsteps++;   //number of LPA steps, change later 
        // }
    }
    std::cout << "Required steps : " << nsteps << std::endl;
}

void write_output(std::string output_filename, DistributedGraph * graph, int rank, int world_size){
    std::ofstream myOutputFile;

    /*Rewrite labels for lpa*/

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
    int ns; // debug purposes 

    // srand(0); // Fix seed ?

    if(argc < 2){
        std::cout << "Please pass the file path as an argument to the program." << std::endl;
        return 1;
    }else{ 
        filename = argv[1];
        ns = atoi(argv[2]); 
        // if(rank == 0)
        //     std::cout << "Working with the file : " << filename << std::endl; 
    }

    /* Create graph and communication handler */
    // DistributedGraph graph; 
    DistributedCompressedGraph graph; 
    CommunicationHandler cm; 

    /* Init graph from input file */
    graph.create_graph_from_METIS(filename);
    
    // /* Init communication structure */
    cm.init_communications(graph.get_ghost_vertices()); 

    // cm.order_ghosts(graph.get_ghost_vertices());

    // RUN LPA 
    // run_LPA(&graph, &cm, ns);
    // run_asynchronous_LPA(&graph, &cm, ns);
    run_with_compression_LPA(&graph, &cm, ns);
    // run_basic_asynchronous_LPA(&graph, &cm, ns);
    

    /* Writes the output of LPA to a file */
    write_output(output_filename, &graph, rank, world_size);
    
    /////////////////////////////////////
    /////* Calculate modularity ?? */////
    /////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD); 
    MPI_Finalize();
}
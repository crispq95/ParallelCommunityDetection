#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <vector>

#include "../../GraphStructure/CommunicationHandler.h"
#include "../../GraphStructure/DistributedGraph.h"
#include "../../GraphStructure/DistributedGraphDPC.h"
#include <unordered_set>


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
    int ns; // debug purposes 

    srand(0); // Fix seed ?

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
    DistributedGraphDPC graph; 
    CommunicationHandler cm; 

    /* Init graph from input file */
    graph.create_graph_for_DCP_from_METIS(filename);
    
    // /* Init communication structure */
    cm.init_communications(graph.get_ghost_vertices()); 

    graph.get_seeds(ns, &cm); 

    // RUN LPA 
    // run_DPC_LPA(&graph, &cm, ns);

    /* Writes the output of LPA to a file */
    // write_output(output_filename, &graph, rank, world_size);
    
    /////////////////////////////////////
    /////* Calculate modularity ?? */////
    /////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD); 
    MPI_Finalize();
}
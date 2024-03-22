#ifndef DistributedGraphDPC_H 
#define DistributedGraphDPC_H

#include "define.h"

// define a class ID_To handle ID_The graph
class DistributedGraphDPC : public DistributedGraph {
    private: 
        ID_T max_degree, max_weight, num_seeds; 

        std::vector<ID_T> degrees, community_seeds, seed_candidates, decreasing_sequence; 
        std::vector<double> quality, density, centrality_index; 
        std::vector<double> second_order_difference;

        std::vector<double> ck_distance_weights, avg_deg; 
        // std::vector<bool> is_locked_by_degree; 

        std::queue<ID_T> candidate_queue; 


        void calculate_quality();
        void calculate_density();
        void calculate_centrality_index();
        void construct_second_order_diff_decreasing_sequence(); 
        void get_number_seeds();
    public: 
        DistributedGraphDPC();
        ~DistributedGraphDPC();

        void set_max_degree( ID_T md ) { max_degree = md; }; 
        ID_T get_max_degree() { return max_degree; };
        // std::vector<ID_T> *get_degrees() { return &degrees; };
        // std::vector<double> *get_quality() { return &quality; };
        // std::vector<double> *get_density() { return &density; };

        void create_graph_for_DCP_from_METIS(std::string filename);

        /* TO-DO : DCP + LPA */
        void find_cores(); 
        void get_seeds( int num_seeds, CommunicationHandler* cm ); 
        void update_ghost_seeds(const std::vector<std::vector<ID_T>>& recv_buffer_ids, const std::vector<std::vector<double>>& recv_buffer_weight); 
        // void run_DPC_LPA_step(); 
};

#endif 
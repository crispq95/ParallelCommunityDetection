#ifndef DEFINE_H
#define DEFINE_H

#include <mpi.h> 
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <math.h> 
#include <iostream>
#include <iomanip>  
#include <ostream> 
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>


typedef unsigned long int T; // unsigned long int ? 
typedef long int LABEL_T;

struct Edge {
    T target; 
    T edge_weight;
};

struct Node{
    T id;
    T node_weight;
    LABEL_T current_label; 
    // bool active; // maybe we can mark unactive nodes so they dont update anymore 
};

struct LocalNode : Node {
    LABEL_T next_label; 
    std::vector<Edge> *edges; 
    bool is_boundary;
};

struct GhostNode : Node {
    int pe_id; // id of the PE it belongs to 
}; 


#endif 
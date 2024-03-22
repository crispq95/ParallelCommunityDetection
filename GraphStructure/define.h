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
#include <numeric>
#include <queue>


typedef unsigned long int ID_T; // unsigned long int ? 
typedef long int LABEL_T;

#define DEFAULT_WEIGTH 1
#define SEED_WEIGTH 2

struct Edge {
    ID_T target; 
    ID_T edge_weight;
};

struct Node{
    ID_T id;
    ID_T node_weight;
    LABEL_T current_label; 
    bool active; // maybe we can mark unactive nodes so they dont update anymore 
};

struct LocalNode : Node {
    LABEL_T next_label; 
    std::vector<Edge> *edges; 
    bool is_boundary;
};

struct GhostNode : Node {
    int pe_id; // id of the PE it belongs to 
}; 

template <typename T>
void sort_indexes_by_label(const std::vector<T> &v, std::vector<int> &ordered_ghosts) {
  ordered_ghosts.resize(v.size());
  iota(ordered_ghosts.begin(), ordered_ghosts.end(), 0);

  stable_sort(ordered_ghosts.begin(), ordered_ghosts.end(),
       [&v](size_t i1, size_t i2) {return v[i1].current_label < v[i2].current_label;});
}

template <typename T>
void sort_decreasing_indexes(const std::vector<T> &v, std::vector<ID_T> &ordered_ghosts) {
  ordered_ghosts.resize(v.size());
  iota(ordered_ghosts.begin(), ordered_ghosts.end(), 0);

  stable_sort(ordered_ghosts.begin(), ordered_ghosts.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
}

#endif 
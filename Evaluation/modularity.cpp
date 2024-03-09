#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;


void addEdge(vector<vector<int>>& adjacencyList, int vertex, int neighbor) {
    adjacencyList[vertex].push_back(neighbor);
    adjacencyList[neighbor].push_back(vertex);  
}


void printGraph(const vector<vector<int>>& adjacencyList) {
    for (int i = 0; i < adjacencyList.size(); ++i) {
        cout << i << ": ";
        for (int j : adjacencyList[i]) {
            cout << j << " ";
        }
        cout << endl;
    }
}

vector<int> readCommunityAssignments(const string& filename) {
    vector<int> communityAssignments;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        return communityAssignments; 
    }

    int community;
    while (file >> community) {
        communityAssignments.push_back(community);
    }

    file.close();
    return communityAssignments;
}


double calculateModularity(const vector<vector<int>>& adjacencyList, const vector<int>& communityAssignments) {
    int numNodes = adjacencyList.size();
    double modularity = 0.0;
    int numEdges = 0;

    for (int i = 0; i < numNodes; ++i) {
        numEdges += adjacencyList[i].size();
    }
    numEdges /= 2;  // Since each edge is counted twice in an undirected graph

    for (int i = 0; i < numNodes; ++i) {
        int community_i = communityAssignments[i];
        for (int j : adjacencyList[i]) {
            int community_j = communityAssignments[j];
            modularity += (community_i == community_j) ? 1.0 : 0.0 - (static_cast<double>(adjacencyList[i].size()) * adjacencyList[j].size()) / (2.0 * numEdges);
        }
    }

    modularity /= (2.0 * numEdges);
    return modularity;
}

int main() {
   
    string filename = "/mnt/c/Users/magda/OneDrive/Desktop/MT/lpa/ParallelCommunityDetection/GraphExamples/0_karate_club_metis.txt";
    string communityFilename = "/mnt/c/Users/magda/OneDrive/Desktop/MT/lpa/ParallelCommunityDetection/output_small_test.txt";

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        return 1;
    }

    int numVertices, numEdges;
    file >> numVertices >> numEdges;

    vector<vector<int>> adjacencyList(numVertices);

    string line;
    getline(file, line);  // To skip to next line

    for (int i = 0; i < numVertices; ++i) {
        getline(file, line);
        istringstream iss(line);
        int neighbor;
        while (iss >> neighbor) {
            addEdge(adjacencyList, i, neighbor);
        }
    }

    file.close();


    vector<int> communityAssignments = readCommunityAssignments(communityFilename);

 
    printGraph(adjacencyList);

    cout << "Community Assignments: ";
    for (int i = 0; i < communityAssignments.size(); ++i) {
        cout << "Node " << i << ": " << communityAssignments[i] << " ";
    }
    cout << endl;

    double modularity = calculateModularity(adjacencyList, communityAssignments);
    cout << "Modularity: " << modularity << endl;

    return 0;
}

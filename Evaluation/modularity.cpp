#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

using namespace std;


void printGraph(const vector<vector<int>>& adjacencyList) {
    for (int i = 0; i < adjacencyList.size(); i++) {
        cout << i << ": ";
        for (int neighbor : adjacencyList[i]) {
            cout << neighbor << " ";
        }
        cout << endl;
    }
}

void printCommunities(const vector<int>& communities) {
    cout << "Communities: ";
    for (size_t i = 0; i < communities.size(); i++) {
        cout << i << ": " << communities[i]  << endl;
    }
}


double calculateModularity(const vector<vector<int>>& adjacencyList, const vector<int>& communityAssignment, int numVertex, int numEdges) {

    double modularity = 0.0;
    vector<int> degrees(numVertex);

    for (int i = 0; i < adjacencyList.size(); i++) {
        int neighbors = adjacencyList[i].size();
        degrees[i] = neighbors;
    }

    /*
    cout << "Degrreess: ";
    for (size_t i = 0; i < degrees.size(); i++) {
        cout << i << ": " << degrees[i]  << endl;
    }
    */

     //Newman-Girvan modularity formula? Only for partitioning to two communities only?

     for (int i = 0; i < adjacencyList.size(); ++i) {
        int ki = adjacencyList[i].size();  // Degree of vertex i
        for (const int& j : adjacencyList[i]) {
            int kj = adjacencyList[j].size();  // Degree of vertex j
            int a = (communityAssignment[i] == communityAssignment[j]) ? 1 : 0;
            modularity += (a - static_cast<double>(ki * kj) / (2 * numEdges));
        }
    }

    modularity /= (2 * numEdges);

    cout<< modularity << endl;

    
  

}

int main() {
    const string filename = "/mnt/c/Users/magda/OneDrive/Desktop/MT/lpa/ParallelCommunityDetection/GraphExamples/0_karate_club_metis.txt";
    const string communityFilename = "/mnt/c/Users/magda/OneDrive/Desktop/MT/lpa/ParallelCommunityDetection/output_small_test.txt";
    const string communityFilename2 = "/mnt/c/Users/magda/OneDrive/Desktop/MT/lpa/ParallelCommunityDetection/out_papers.txt";

    ifstream inputFile(filename);
    ifstream communityFile(communityFilename2);

     if (!inputFile.is_open() || !communityFile.is_open()) {
        cerr << "Error opening the file!" << endl;
        return 1;
    }

    int vertices, edges;
    inputFile >> vertices >> edges;

    string line;
    getline(inputFile, line);  // fix it later 

    vector<vector<int>> adjacencyList(vertices);

    for (int i = 0; i < vertices; i++) {
        string line;
        getline(inputFile, line);

        istringstream iss(line);
        int neighbor;

        while (iss >> neighbor) {
            adjacencyList[i].push_back(neighbor);
        }
    }

    //printGraph(adjacencyList);

    vector<int> communities(vertices, -1); 

    int vertex, community;
    while (communityFile >> vertex >> community) {
        communities[vertex] = community;
    }

    //printCommunities(communities);

    inputFile.close();
    communityFile.close();


    calculateModularity(adjacencyList, communities, vertices, edges);



    return 0;
}

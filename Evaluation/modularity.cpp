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
    for (size_t i = 0; i < adjacencyList.size(); i++) {
        cout << i << ": ";
        for (int neighbor : adjacencyList[i]) {
            cout << neighbor << " ";
        }
        cout << endl;
    }
}

int main() {
    const string filename = "/mnt/c/Users/magda/OneDrive/Desktop/MT/lpa/ParallelCommunityDetection/GraphExamples/0_karate_club_metis.txt";
    const string communityFilename = "/mnt/c/Users/magda/OneDrive/Desktop/MT/lpa/ParallelCommunityDetection/output_small_test.txt";

    ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        cerr << "Error opening the file: " << filename << endl;
        return 1;
    }

    int vertices, edges;
    inputFile >> vertices >> edges;

     string line;
     getline(inputFile, line);

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

    printGraph(adjacencyList);

    inputFile.close();

    return 0;
}

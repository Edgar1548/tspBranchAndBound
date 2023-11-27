#include <iostream>
#include <string.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <climits>
#include <list>

using namespace std;
double calculateDistance(int x1, int y1, int x2, int y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

vector<vector<double>> convertToMatrix(const vector<pair<int, int>>& coordinates) {
    int n = coordinates.size();
    vector<vector<double>> distanceMatrix(n, vector<double>(n, 0.0));
    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
            if (a == b) distanceMatrix[a][b] = -1;
            else distanceMatrix[a][b] = calculateDistance(coordinates[a].first, coordinates[a].second, coordinates[b].first, coordinates[b].second);
        }
    }
    return distanceMatrix;
}

bool compareLists(const std::list<int>& list1, const std::list<int>& list2) {
    return list1 == list2;
}

void includeEdge(vector<vector<double>> &matrix, pair<int, int> edge, vector<list<int>>& linkedLists){
    int n = matrix.size(), a = edge.first, b = edge.second;
    if (matrix[a][b] == -1) return;
    for (int i = 0; i < n; ++i) matrix[a][i] = -1;
    for (int i = 0; i < n; ++i)  matrix[i][b] = -1;
    matrix[b][a] = -1;
    int aFound = -1, bFound = -1;
    for (int i=0; i<linkedLists.size(); i++) {
        if (!linkedLists[i].empty() && linkedLists[i].back() == a) {
            matrix[b][linkedLists[i].front()] = -1;
            aFound = i;
            break;
        }
    }
    for (int i=0; i<linkedLists.size(); i++) {
        if (!linkedLists[i].empty() && linkedLists[i].front() == b) {
            matrix[linkedLists[i].back()][a] = -1;
            bFound = i;
            break;
        }
    }
    if (aFound != -1 && bFound != -1) {
        linkedLists[aFound].insert(linkedLists[aFound].end(), linkedLists[bFound].begin(), linkedLists[bFound].end());
        matrix[linkedLists[aFound].back()][linkedLists[aFound].front()] = -1;
        if (aFound != bFound){
            auto it = linkedLists.begin();
            while (it != linkedLists.end()) {
                if (compareLists(*it, linkedLists[bFound])) {
                    linkedLists.erase(it);
                    break;
                } else {
                    ++it;
                }
            }
        }
    } else if (aFound != -1) {
        linkedLists[aFound].push_back(b);
    } else if (bFound != -1) {
        linkedLists[bFound].push_front(a);
    } else {
        linkedLists.push_back({a, b});
    }
}

void notIncludeEdge(vector<vector<double>> &matrix, pair<int, int> edge){
    matrix[edge.first][edge.second] = -1;
}

int reduceMatrix(vector<vector<double>> &matrix, vector<pair<int, int>> &edges, bool makeOtherThing = false){
    int cost = 0;
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        double minDist = INT_MAX;
        for (int j = 0; j < n; ++j) {
            if (matrix[i][j] != -1 && matrix[i][j] < minDist) minDist = matrix[i][j];
        }
        if (minDist != INT_MAX){
            for (int j = 0; j < n; ++j) {
                if (matrix[i][j] != -1) matrix[i][j] -= minDist;
            }
            cost += minDist;
        }
    }
    for (int j = 0; j < n; ++j) {
        double minDist = INT_MAX;
        for (int i = 0; i < n; ++i) {
            if (matrix[i][j] != -1 && matrix[i][j] < minDist) minDist = matrix[i][j];
        }
        if (minDist != INT_MAX){
            for (int i = 0; i < n; ++i) {
                if (matrix[i][j] != -1){
                    matrix[i][j] -= minDist;
                    if (makeOtherThing && matrix[i][j] == 0) edges.push_back({i, j});
                }
            }
            cost += minDist;
        }
    }
    return cost;
}
int findCost(vector<vector<double>> &matrix, pair<int, int> &edge){
    int cost = 0;
    int n = matrix.size();
    int i = edge.first;
    double minDist = INT_MAX;
    for (int j = 0; j < n; ++j) {
        if (matrix[i][j] != -1 && matrix[i][j] < minDist) minDist = matrix[i][j];
    }
    if (minDist != INT_MAX){
        cost += minDist;
    }
    int j = edge.second;
    minDist = INT_MAX;
    for (int i = 0; i < n; ++i) {
        if (matrix[i][j] != -1 && matrix[i][j] < minDist) minDist = matrix[i][j];
    }
    if (minDist != INT_MAX){
        cost += minDist;
    }
    
    return cost;
}

double findMinPath(vector<vector<double>> distanceMatrix){
    auto originalMatrix = distanceMatrix;
    vector<pair<int, int>> edges;
    vector<list<int>> pathEdges;
    double minCostPath = 0;
    double matrixCost = reduceMatrix(distanceMatrix, edges, true);
    while(edges.size() > 0){
        double maxCost = 0;
        pair<int, int> maxEdge;
        for (auto it: edges){
            auto distanceMatrixModified = distanceMatrix;
            distanceMatrixModified[it.first][it.second] = -1;
            double tempCost = findCost(distanceMatrixModified, it);
            if (maxCost < tempCost) {
                maxCost = tempCost;
                maxEdge = it;
            }
        }   
        if (maxCost == 0) {
            cout << edges.size() << "\n";
            for (auto it: edges){
                cout << it.first << " " << it.second << " cost: ";
                cout << originalMatrix[it.first][it.second] << "\n";
            }
            break;  
        }
        auto distanceMatrixInclude = distanceMatrix;
        auto distanceMatrixNotInclude = distanceMatrix;
        vector<pair<int, int>> edgesInclude;
        vector< pair<int, int>> edgesNotInclude;
        auto includePathEdge = pathEdges;
        auto notIncludePathEdge = pathEdges;
        includeEdge(distanceMatrixInclude, maxEdge, includePathEdge);
        double includeCost = reduceMatrix(distanceMatrixInclude, edgesInclude, true);
        includeCost += matrixCost;

        notIncludeEdge(distanceMatrixNotInclude, maxEdge);
        double notIncludeCost = reduceMatrix(distanceMatrixNotInclude, edgesNotInclude, true);
        notIncludeCost += matrixCost;

        distanceMatrix = includeCost <= notIncludeCost ? distanceMatrixInclude : distanceMatrixNotInclude;
        edges = includeCost <= notIncludeCost ?  edgesInclude : edgesNotInclude; 
        matrixCost = includeCost <= notIncludeCost ? includeCost : notIncludeCost;
        pathEdges = includeCost <= notIncludeCost ? includePathEdge : notIncludePathEdge;
        cout << includeCost << " - " << notIncludeCost << "\n";
    }
    int i;
    for (auto it: pathEdges) {
        list<int>::const_iterator iterator;
        list<int>::const_iterator iteratorAfter = it.begin();
        for (iterator = it.begin(), i = 0; iterator != it.end() && i < it.size()-1; ++iterator, ++i) {
            ++iteratorAfter;
            minCostPath += originalMatrix[*iterator][*iteratorAfter];
        }
    }
    cout << matrixCost << "\n";
    return minCostPath;
}

int main(int argc, char** argv) {    
    vector<pair<int, int>> coords;
	string file_name;
    if (argc != 2) {
		cout << "No argument provided." << endl;
		return -1;
	}
	else {
		string line;
		ifstream myfile(argv[1]);
		if(!myfile) {
			cout << "Cannot open the file." << endl;
			return -1;
		}
		int line_count = 0;
		int a, b, c;
		string temp, inp;
		while(getline(myfile,line)) {
			if (line_count < 6) {
				istringstream iss(line);
				iss >> temp >> temp >> inp;
				if (line_count == 0) {
					file_name = inp;
				}
				line_count +=1;
				continue;
			}
			istringstream iss(line);
			iss >> a >> b >> c;
			if(a!=0) {
				pair<int, int> coord = {b, c};
				coords.push_back(coord);
			}
			line_count += 1;
		}
		myfile.close(); 
	}
    auto distanceMatrix = convertToMatrix(coords);
    auto minPath = findMinPath(distanceMatrix);
    cout << minPath;
}
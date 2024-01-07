#include "Graph.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>

/**
* Constructor implementations
*/
Graph::Graph(vector<vector<int>>& adjlist)
{
    this->adjlist = adjlist;
}

Graph::Graph(const string& graphString)
{
    string graphStringCopy = graphString;

    graphStringCopy.erase(remove(graphStringCopy.begin(), graphStringCopy.end(), ' '), graphStringCopy.end());
    graphStringCopy.erase(remove(graphStringCopy.begin(), graphStringCopy.end(), ','), graphStringCopy.end());

    istringstream iss(graphStringCopy);
    vector<string> lines;
    string line;

    while (getline(iss, line)) {
        lines.push_back(line);
    }

    sort(lines.begin(), lines.end());

    for (int i = 0; i < lines.size(); i++) {
        vector<int> currentAdjacencies;
        for (int j = 2; j < lines[i].size(); j++) {
            char current = lines[i][j];
            currentAdjacencies.push_back(static_cast<int>(current - 'a'));
        }
        adjlist.push_back(currentAdjacencies);
    }
}

/**
* Getter
*/
const vector<vector<int>>& Graph::getAdjlist() const
{
    return adjlist;
}

/**
* Returns a graph, that only contains the vertices in a given set X. All edges
* are calculated accordingly.
*
* @param X The set
* @param nodeMapping Mapping from old indices to new indices (updated in the method)
* @return A subgraph that only contains the vertices given in X.
*/
Graph Graph::getSubGraph(const unordered_set<int>& X, vector<int>& nodeMapping) const
{
    vector<vector<int>> newAdjList;

    // Create a mapping from original indices to subgraph indices
    int subgraphIndex = 0;

    for (int i : X) {
        nodeMapping.push_back(i);
    }

    for (int i : X) {
        vector<int> neighbors;
        for (int neighbor : adjlist[i]) {
            auto it = find(X.begin(), X.end(), neighbor);
            if (it != X.end()) {
                neighbors.push_back(distance(X.begin(), it));
            }
        }
        newAdjList.push_back(neighbors);
    }

    return Graph(newAdjList);
}

/**
* Uses BFS to check if the graph is connected
*/
bool Graph::isConnected() const {
    vector<bool> visited(adjlist.size(), false);
    queue<int> q;

    q.push(0);
    visited[0] = true;

    while (!q.empty()) {
        int current = q.front();
        q.pop();

        for (int neighbor : adjlist[current]) {
            if (!visited[neighbor]) {
                q.push(neighbor);
                visited[neighbor] = true;
            }
        }
    }

    return all_of(visited.begin(), visited.end(), [](bool v) { return v; });
}

/**
* Uses DFS to retrieve all components of the graph.
*/
vector<unordered_set<int>> Graph::getConnectedComponents() const {
    vector<unordered_set<int>> components;
    vector<bool> visited(adjlist.size(), false);

    for (int i = 0; i < adjlist.size(); ++i) {
        if (!visited[i]) {
            unordered_set<int> componentNodes;
            dfs(i, visited, componentNodes);
            components.push_back(componentNodes);
        }
    }

    return components;
}

/**
* A simple recursive DFS algorithm.
*
* @param node The current node.
* @param visited Contains information about which nodes have already been visited.
* @param componentNodes This set is used to store the results.
*/
void Graph::dfs(int node, vector<bool>& visited, unordered_set<int>& componentNodes) const {
    visited[node] = true;
    componentNodes.insert(node);

    for (int neighbor : adjlist[node]) {
        if (!visited[neighbor]) {
            dfs(neighbor, visited, componentNodes);
        }
    }
}


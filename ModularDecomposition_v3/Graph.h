#pragma once
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

using namespace std;

class Graph
{
private:
    vector<vector<int>> adjlist;

public:
    Graph(vector<vector<int>>& adjlist);
    Graph(const string& graphString);
    const vector<vector<int>>& getAdjlist() const;
    Graph getSubGraph(const unordered_set<int>& X, vector<int>& nodeMapping) const;

    bool isConnected() const;
    vector<unordered_set<int>> getConnectedComponents() const;

private:
    void dfs(int node, vector<bool>& visited, unordered_set<int>& componentNodes) const;
};

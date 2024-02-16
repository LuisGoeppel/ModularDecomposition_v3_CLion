#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <vector>
#include <random>

#include "MD_Tree.h"
#include "Graph.h"

using namespace std;

class Util {
public:
    static void sortTree(MD_Tree& tree);
    static vector<int> rewriteAdjacencyList(string& adjList);
    static bool isConsecutivelyOrdered(const string& input);
    static bool testModularDecompositionTree(const Graph& graph, const MD_Tree& tree);
    static void checkChildNodeValues(MD_Tree& tree);
    static int getNumberVertices(const Graph& graph);
    static int getNumberEdges(const Graph& graph);
    static MD_Tree createRandomModularDecompositionTree(int nVertices, bool isCoGraph);
    static Graph createGraphFromTree(const MD_Tree& tree);

private:
    static void sortTreeHelper(TreeNode* node);
    static vector<int> getInverse(const vector<int>& vec);
    static TreeNode* getCommonAncestorChildLhs(const MD_Tree& tree, int lhs, int rhs);
    static TreeNode* getCorrspondingTreeNode(TreeNode* currentNode, int node);
    static void checkChildNodeValuesHelper(TreeNode* node);
    static void removeFalseInnerNodes(TreeNode* currentNode, TreeNode* callingNode, bool isParent, bool& removedNode);
    static void getHighestVertexValue(TreeNode* node, int& highestValue);
    static Label getLabel(int n);
};


#pragma once
#include <iostream>
#include <vector>
#include <unordered_set>
#include <set>

using namespace std;

enum Label {
    PARALLEL,
    SERIES,
    PRIME,
    LEAF
};


struct TreeNode {
    Label label;
    int value;

    TreeNode* child;
    TreeNode* sibling;
    TreeNode* parent;

    bool markedLeft;
    bool markedRight;

    TreeNode(int val);
    TreeNode(Label l);

};

struct MD_Tree {
    TreeNode* root;
    
    MD_Tree* left;
    MD_Tree* right;

    int leftIndex;
    int rightIndex;

    MD_Tree(TreeNode* r);
};

bool compareTreeNodePointers(const TreeNode* lhs, const TreeNode* rhs);

bool operator==(const MD_Tree& lhs, const TreeNode* rhs);
bool operator!=(const MD_Tree& lhs, const TreeNode* rhs);
bool operator==(const TreeNode& lhs, const TreeNode& rhs);

void setChild(TreeNode* lhs, TreeNode* rhs);
void setSibling(TreeNode* lhs, TreeNode* rhs);

void setNeighbor(MD_Tree* lhs, MD_Tree* rhs);

void printTree(const TreeNode* node, int depth = 0);
vector<int> getPreOrderLeafs(const TreeNode* root);
void getMaxContSubTrees(TreeNode* node, unordered_set<int>& X,
    unordered_set<TreeNode*>& subTrees);

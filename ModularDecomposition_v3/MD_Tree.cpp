#include "MD_Tree.h"

/**
* Constructor implementations
*/
TreeNode::TreeNode(int val) : value(val), label(LEAF), markedLeft(false),
markedRight(false), child(nullptr), sibling(nullptr), parent(nullptr) {}

TreeNode::TreeNode(Label l) : value(-1), label(l), markedLeft(false),
markedRight(false), child(nullptr), sibling(nullptr), parent(nullptr) {}

MD_Tree::MD_Tree(TreeNode* r) : root(r), leftIndex(-1), rightIndex(-1),
left(nullptr), right(nullptr) {}

/**
* Operator implementations
*/

bool compareTreeNodePointers(const TreeNode* lhs, const TreeNode* rhs)
{
    if (lhs->value == rhs->value) {
        return lhs->label < rhs->label;
    }
    return lhs->value < rhs->value;
}
bool operator==(const MD_Tree& lhs, const TreeNode* rhs) {
    return lhs.root == rhs;
}

bool operator!=(const MD_Tree& lhs, const TreeNode* rhs) {
    return !(lhs == rhs);
}

bool operator==(const TreeNode& lhs, const TreeNode& rhs) {
    return lhs.value == rhs.value && lhs.label == rhs.label && lhs.child == rhs.child
        && lhs.sibling == rhs.sibling && lhs.parent == rhs.parent;
}

/**
* Sets the child of a given TreeNode.
*
* @param lhs The parent
* @param rhs The child
*/
void setChild(TreeNode* lhs, TreeNode* rhs)
{
    lhs->child = rhs;
    rhs->parent = lhs;
}

/**
* Sets the sibling of a given TreeNode.
*
* @param lhs The "left" sibling.
* @param rhs The "right" sibling.
*/
void setSibling(TreeNode* lhs, TreeNode* rhs)
{
    lhs->sibling = rhs;
    rhs->parent = lhs->parent;
}

void setNeighbor(MD_Tree* lhs, MD_Tree* rhs)
{
    if (lhs != nullptr) {
        lhs->right = rhs;
    }
    if (rhs != nullptr) {
        rhs->left = lhs;
    }
}

/**
* Recursively prints a MD_Tree on the command line
*
* @param node The starting node for the recursive process.
* @param depth The current depth of the tree (default: 0)
*/
void printTree(const TreeNode* node, int depth) {
    for (int i = 0; i < depth; ++i) {
        cout << "  ";
    }

    string mark = "";
    if (node->markedLeft && node->markedRight) {
        mark = ", Mark: Left & Right";
    }
    else if (node->markedLeft) {
        mark = ", Mark: Left";
    }
    else if (node->markedRight) {
        mark = ", Mark: Right";
    }

    if (node->label == LEAF) {
        cout << "Leaf: " << static_cast<char>(node->value + 'a') << mark << endl;
    }
    else {
        cout << "Label: " << node->label << mark << endl;
    }

    if (node->child) {
        printTree(node->child, depth + 1);
    }
    if (depth > 0) {
        if (node->sibling) {
            printTree(node->sibling, depth);
        }
    }
}

/**
* A helper function to get a pre-Order of all leafs in the tree.
*
* @param node The current node.
* @param checkSiblings If the pre-Order of the nodes siblings should be returned.
* @return A pre-Order of the tree's leafs.
*/
vector<int> getPreOrderLeafsHelper(const TreeNode* node, bool checkSiblings) {
    vector<int> result;

    if (node == nullptr) {
        return result;
    }

    if (node->label == LEAF) {
        result.push_back(node->value);
    }

    auto childResult = getPreOrderLeafsHelper(node->child, true);
    result.insert(result.end(), childResult.begin(), childResult.end());

    if (checkSiblings) {
        auto siblingResult = getPreOrderLeafsHelper(node->sibling, true);
        result.insert(result.end(), siblingResult.begin(), siblingResult.end());
    }

    return result;
}

/**
* Return a pre-Order of the given tree's leafs.
*
* @param root The root node of the tree.
* @return A pre-Order of the tree's leafs.
*/
vector<int> getPreOrderLeafs(const TreeNode* root) {
    return getPreOrderLeafsHelper(root, false);
}

/**
* Returns the number of matching arguments of a vector and an unordered_set.
*
* @tparam T The type of elements in the vector and set.
* @param lhs The vector
* @param rhs The set
* @return The number of matching parameters.
*/
template <typename t>
int getNMatchingArguments(const vector<t>& lhs, const unordered_set<t>& rhs) {
    int count = 0;
    for (int element : lhs) {
        if (rhs.find(element) != rhs.end()) {
            count++;
        }
    }
    return count;
}

/**
* Calculates the maximal containing subtrees, based on a set X. This means, that all nodes are returned,
* whose children can all be found in a set X. This cannot be true for the node's parent.
*
* @param node The root node of the tree to check
* @param X The set
* @param subTrees This set will be filled with the maximal containing subtrees.
*/
void getMaxContSubTrees(TreeNode* node, unordered_set<int>& X, unordered_set<TreeNode*>& subTrees)
{
    vector<int> leafs = getPreOrderLeafs(node);
    int nMatching = getNMatchingArguments(leafs, X);
    if (nMatching == leafs.size()) {
        subTrees.insert(node);
        for (int leaf : leafs) {
            X.erase(leaf);
        }
    }
    else if (nMatching > 0 && node->child != nullptr) {
        getMaxContSubTrees(node->child, X, subTrees);
    }
    if (node->sibling != nullptr) {
        getMaxContSubTrees(node->sibling, X, subTrees);
    }
}

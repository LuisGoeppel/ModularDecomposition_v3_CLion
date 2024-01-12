#include "Util.h"


/**
* Sorts the given tree alphabetically.
*
* @param tree The tree to sort.
*/
void Util::sortTree(MD_Tree& tree) {
    sortTreeHelper(tree.root);
}

/**
 * Rewrites the provided adjacency list, relabeling nodes consecutively and returning
 * a vector of differences between ASCII values of old and new node names.
 * The function modifies the input adjacency list to maintain the graph structure while
 * assigning consecutive labels to nodes.
 *
 * @param adjList A reference to the input adjacency list (modified in the method).
 * @return A vector of integers representing differences between ASCII values of old and new node names.
 */
vector<int> Util::rewriteAdjacencyList(string& adjList)
{
    vector<int> indexMapping;
    string updatedAdjlist = "";

    istringstream iss(adjList);
    vector<string> lines;
    string line;

    while (getline(iss, line)) {
        lines.push_back(line);
    }

    sort(lines.begin(), lines.end());

    for (int i = 0; i < lines.size(); i++) {
        char lineChar = lines[i][0];
        int lineLetter = lineChar - 'a';
        indexMapping.push_back(lineLetter);
        updatedAdjlist += lines[i] + "\n";
    }
    vector<int> inverseIndices = getInverse(indexMapping);

    for (int i = 0; i < updatedAdjlist.size(); i++) {
        int current = updatedAdjlist[i] - 'a';
        if (current >= 0 && current < inverseIndices.size()) {
            updatedAdjlist[i] = 'a' + inverseIndices[current];
        }
    }
    adjList = updatedAdjlist;

    return indexMapping;
}

bool Util::isConsecutivelyOrdered(const string& input)
{
    stringstream ss(input);
    string line;
    char expectedChar = 'a';

    while (std::getline(ss, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] != expectedChar) {
            return false;
        }
        ++expectedChar;
    }

    return true;
}

bool Util::testModularDecompositionTree(const Graph& graph, const MD_Tree& tree)
{
    vector<vector<int>> adjlist = graph.getAdjlist();
    vector<vector<int>> isConnectedByPrimeNode(adjlist.size(), vector<int>(adjlist.size(), 0));

    for (int lhs = 0; lhs < adjlist.size(); lhs++) {
        for (int rhs = 0; rhs < adjlist.size(); rhs++) {
            if (lhs != rhs) {

                bool hasConnectionGraph = false;
                for (int adj : adjlist[lhs]) {
                    if (adj == rhs) {
                        hasConnectionGraph = true;
                        break;
                    }
                }

                bool hasConnectionTree = hasConnectionGraph;
                TreeNode* ancestorChildLhs = getCommonAncestorChildLhs(tree, lhs, rhs);
                TreeNode* ancestorChildRhs = getCommonAncestorChildLhs(tree, rhs, lhs);
                Label commonAncestorLabel = ancestorChildLhs->parent->label;

                if (commonAncestorLabel == PARALLEL) {
                    hasConnectionTree = false;
                }
                else if (commonAncestorLabel == SERIES) {
                    hasConnectionTree = true;
                }
                else {
                    if (isConnectedByPrimeNode[lhs][rhs] != 0) {
                        hasConnectionTree = (isConnectedByPrimeNode[lhs][rhs] == 1);
                    }
                    else {
                        vector<int> valuesLhs = getPreOrderLeafs(ancestorChildLhs);
                        vector<int> valuesRhs = getPreOrderLeafs(ancestorChildRhs);
                        int matrixValue = hasConnectionGraph ? 1 : -1;
                        for (int i : valuesLhs) {
                            for (int j : valuesRhs) {
                                isConnectedByPrimeNode[i][j] = matrixValue;
                                isConnectedByPrimeNode[j][i] = matrixValue;
                            }
                        }
                    }
                }

                if (hasConnectionGraph != hasConnectionTree) {
                    cout << lhs << ", " << rhs << endl;
                    return false;
                }
            }
        }
    }
    return true;
}


/**
* A helper function for the sortTree algorithm. It takes a node, sorts its children
* alphabetically, and then recurisvely calls itself for the node' children and siblings.
*
* @param node The current node, of which the children should get sorted.
*/
void Util::sortTreeHelper(TreeNode* node) {
    vector<TreeNode*> children;
    TreeNode* next = node->child;
    while (next != nullptr) {
        children.push_back(next);
        next = next->sibling;
    }
    if (children.size() > 0) {

        sort(children.begin(), children.end(), compareTreeNodePointers);
        for (int i = 0; i < children.size(); i++) {
            if (i == 0) {
                setChild(node, children[i]);
            }
            else {
                setSibling(children[i - 1], children[i]);
            }
        }
        children[children.size() - 1]->sibling = nullptr;
    }

    if (node->sibling != nullptr) {
        sortTreeHelper(node->sibling);
    }
    if (node->child != nullptr) {
        sortTreeHelper(node->child);
    }
}

vector<int> Util::getInverse(const vector<int>& vec)
{
    int maxElement = 0;
    for (int element : vec) {
        if (element > maxElement) {
            maxElement = element;
        }
    }
    vector<int> inverseVector(maxElement + 1, -1);
    for (int i = 0; i < vec.size(); i++) {
        inverseVector[vec[i]] = i;
    }
    return inverseVector;
}

TreeNode* Util::getCommonAncestorChildLhs(const MD_Tree& tree, int lhs, int rhs)
{
    TreeNode* lhsNode = getCorrspondingTreeNode(tree.root, lhs);
    TreeNode* rhsNode = getCorrspondingTreeNode(tree.root, rhs);
    TreeNode* currentAncestor = lhsNode->parent;
    TreeNode* ancestorChildLhs = lhsNode;

    while (currentAncestor != nullptr) {
        vector<int> currentLeafNodes = getPreOrderLeafs(currentAncestor);
        if (find(currentLeafNodes.begin(), currentLeafNodes.end(), rhs) != currentLeafNodes.end()) {
            return ancestorChildLhs;
        }
        ancestorChildLhs = currentAncestor;
        currentAncestor = currentAncestor->parent;
    }
}

TreeNode* Util::getCorrspondingTreeNode(TreeNode* currentNode, int node)
{
    if (currentNode->value == node) {
        return currentNode;
    }
    TreeNode* out = nullptr;
    if (currentNode->child != nullptr) {
        out = getCorrspondingTreeNode(currentNode->child, node);
    }
    if (currentNode->sibling != nullptr && out == nullptr) {
        out = getCorrspondingTreeNode(currentNode->sibling, node);
    }
    return out;
}





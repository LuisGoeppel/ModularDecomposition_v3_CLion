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
vector<int> Util::rewriteAdjacencyList(string& adjList) {
    vector<int> indexMapping;
    string updatedAdjlist = "";

    istringstream iss(adjList);
    vector<string> lines;
    string line;

    while (getline(iss, line)) {
        lines.push_back(line);
    }

    sort(lines.begin(), lines.end(), [](const string& a, const string& b) {
        return stoi(a) < stoi(b);
    });

    for (int i = 0; i < lines.size(); i++) {
        int lineNum = stoi(lines[i]);
        indexMapping.push_back(lineNum);
        updatedAdjlist += to_string(i) + ": " + lines[i].substr(lines[i].find(':') + 2) + "\n";
    }

    adjList = updatedAdjlist;

    return indexMapping;
}

/**
 * Checks if the provided string of node names is consecutively ordered from 'a' onwards.
 *
 * @param input The string of node names.
 * @return True if consecutively ordered, false otherwise.
 */
bool Util::isConsecutivelyOrdered(const string& input) {
    istringstream iss(input);
    string line;
    int prevNumber = -1;

    while (getline(iss, line)) {
        std::istringstream lineStream(line);
        int currentNumber;
        char colon;
        if (lineStream >> currentNumber >> colon) {
            if (currentNumber != prevNumber + 1) {
                return false;
            }
            prevNumber = currentNumber;
        } else {
            return false;
        }
    }
    return true;
}

/**
 * Tests if the provided Modular Decomposition tree is valid for the given graph.
 *
 * @param graph The graph.
 * @param tree The Modular Decomposition tree.
 * @return True if the tree is valid, false otherwise.
 */
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
                    cout << "Error! Nodes " << lhs << " and " << rhs  << " have different connections in the graph and the MD_Tree!" << endl;
                    return false;
                }
            }
        }
    }
    return true;
}

/**
 * Checks the child node values of the provided Modular Decomposition tree.
 *
 * @param tree The Modular Decomposition tree.
 */
void Util::checkChildNodeValues(MD_Tree& tree)
{
    checkChildNodeValuesHelper(tree.root);
}

/**
 * @brief A helper function for the sortTree algorithm. It takes a node, sorts its children
 * alphabetically, and then recursively calls itself for the node's children and siblings.
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

/**
 * Gets the inverse of the provided vector.
 *
 * @param vec The input vector.
 * @return The inverse vector.
 */
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

/**
* Gets the common ancestor child node for the given nodes in the Modular Decomposition tree.
*
* @param tree The Modular Decomposition tree.
* @param lhs The left - hand side node.
* @param rhs The right - hand side node.
* @return The common ancestor child node.
*/
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
    return nullptr;
}

/**
 * Gets the corresponding tree node for the given node in the Modular Decomposition tree.
 *
 * @param currentNode The current node to search from.
 * @param node The target node.
 * @return The corresponding tree node.
 */
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

/**
 * A helper function to check child node values of the Modular Decomposition tree.
 *
 * @param node The current node.
 */
void Util::checkChildNodeValuesHelper(TreeNode* node)
{
    TreeNode* currentChild = node->child;
    int childCount = 0;

    while (currentChild != nullptr) {
        childCount++;
        currentChild = currentChild->sibling;
    }

    if (childCount != node->nChildNodes) {
        cout << "Error at child count: Expected " << node->nChildNodes
            << " but found " << childCount << " nodes!";
    }

    if (node->sibling != nullptr) {
        checkChildNodeValuesHelper(node->sibling);
    }
    if (node->child != nullptr) {
        checkChildNodeValuesHelper(node->child);
    }
}

int Util::getNumberVertices(const Graph &graph) {
    return graph.getAdjlist().size();
}

int Util::getNumberEdges(const Graph &graph) {
    int count = 0;
    for (int i = 0; i < graph.getAdjlist().size(); i++) {
        count += graph.getAdjlist()[i].size();
    }
    return count/2;
}

MD_Tree Util::createRandomModularDecompositionTree(int nVertices, bool isCoGraph) {

    random_device rd;
    mt19937 gen(rd());
    int upperBoundNodeTypes = isCoGraph ? 1 : 2;
    uniform_int_distribution<int> randomNodeType(0, upperBoundNodeTypes);
    uniform_int_distribution<int> randomRootNodeType(1, upperBoundNodeTypes);
    uniform_int_distribution<int> randomAmountChildNodes(2, 5);
    uniform_int_distribution<int> randomBool(0, 1);
    int vertexCount;
    TreeNode* rootNode;

    do {
        vertexCount = 0;
        vector<int> remainingVertices(nVertices);
        iota(remainingVertices.begin(), remainingVertices.end(), 0);
        queue<TreeNode*> remainingInnerNodes;
        rootNode = new TreeNode(getLabel(randomRootNodeType(gen)));
        remainingInnerNodes.push(rootNode);

        while (!remainingInnerNodes.empty()) {
            TreeNode* currentNode = remainingInnerNodes.front();
            remainingInnerNodes.pop();

            int nChildren = min(randomAmountChildNodes(gen), static_cast<int>(remainingVertices.size()));
            if (nChildren <= 2 && currentNode->label == PRIME) {
                currentNode->label = SERIES;
            }
            currentNode->nChildNodes = nChildren;
            TreeNode* currentChild;
            TreeNode* lastChild;
            for (int i = 0; i < nChildren; i++) {
                int isLeafNode = remainingVertices.size() >= 2 ? randomBool(gen) : 1;
                if (isLeafNode == 1) {
                    uniform_int_distribution<int> randomVertex(0, remainingVertices.size() - 1);
                    int vertexIndex = randomVertex(gen);
                    int vertex = remainingVertices[vertexIndex];
                    remainingVertices.erase(remainingVertices.begin() + vertexIndex);
                    vertexCount++;
                    currentChild = new TreeNode(vertex);
                } else {
                    Label currentLabel;
                    do {
                        currentLabel = getLabel(randomNodeType(gen));
                    } while (currentLabel == currentNode->label);
                    currentChild = new TreeNode(currentLabel);
                    remainingInnerNodes.push(currentChild);
                }
                if (i == 0) {
                    setChild(currentNode, currentChild);
                } else {
                    setSibling(lastChild, currentChild);
                }
                lastChild = currentChild;
            }
        }
    } while (vertexCount < nVertices);

    bool removedNodes = false;
    do {
        removedNodes = false;
        removeFalseInnerNodes(rootNode, rootNode, false, removedNodes);
    } while (removedNodes);

    MD_Tree* result = new MD_Tree(rootNode);
    return *result;
}

Graph Util::createGraphFromTree(const MD_Tree& tree) {
    int nVertices = 0;
    getHighestVertexValue(tree.root, nVertices);
    nVertices++;
    vector<vector<int>> adjList(nVertices);
    for (int lhs = 0; lhs < nVertices; lhs++) {
        for (int rhs = 0; rhs < nVertices; rhs++) {
            if (lhs != rhs) {

                TreeNode* lhsNode = getCorrspondingTreeNode(tree.root, lhs);
                TreeNode* currentAncestor = lhsNode->parent;
                while (currentAncestor != nullptr) {
                    vector<int> currentLeafNodes = getPreOrderLeafs(currentAncestor);
                    if (find(currentLeafNodes.begin(), currentLeafNodes.end(), rhs) != currentLeafNodes.end()) {
                        break;
                    }
                    currentAncestor = currentAncestor->parent;
                }

                switch(currentAncestor->label) {
                    case PARALLEL:
                        // do nothing
                        break;
                    case SERIES:
                        adjList[lhs].push_back(rhs);
                        break;
                    case PRIME:
                        adjList[lhs].push_back(rhs);
                        break;
                    case LEAF:
                        // do nothing;
                        break;
                }
            }
        }
    }
    return Graph(adjList);
}

Label Util::getLabel(int n) {
    switch (n) {
        case 0:
            return PARALLEL;
        case 1:
            return SERIES;
        case 2:
            return PRIME;
        default:
            return static_cast<Label>(NULL);
    }
}

void Util::removeFalseInnerNodes(TreeNode* currentNode, TreeNode* callingNode, bool isParent, bool& removedNode) {
    if (currentNode->label != LEAF) {
        int nChildren = 0;
        TreeNode* nextChild = currentNode->child;
        while(nextChild != nullptr) {
            nChildren++;
            nextChild = nextChild->sibling;
        }
        if (nChildren == 0) {
            if (isParent) {
                setChild(callingNode, currentNode->sibling);
            } else {
                setSibling(callingNode, currentNode->sibling);
            }
            removedNode = true;
        } else if (nChildren == 1) {
            if (isParent) {
                setChild(callingNode, currentNode->child);
                setSibling(currentNode->child, currentNode->sibling);
            } else {
                setSibling(callingNode, currentNode->child);
                setSibling(currentNode->child, currentNode->sibling);
            }
            removedNode = true;
        }
    }

    if (currentNode->sibling != nullptr) {
        removeFalseInnerNodes(currentNode->sibling, currentNode, false, removedNode);
    }
    if (currentNode->child != nullptr) {
        removeFalseInnerNodes(currentNode->child, currentNode, true, removedNode);
    }
}

void Util::getHighestVertexValue(TreeNode* node, int& highestValue) {
    if (node->label == LEAF) {
        if (node->value > highestValue) {
            highestValue = node->value;
        }
    }

    if (node->sibling != nullptr) {
        getHighestVertexValue(node->sibling, highestValue);
    }
    if (node->child != nullptr) {
        getHighestVertexValue(node->child, highestValue);
    }
}

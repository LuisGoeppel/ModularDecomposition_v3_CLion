#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <memory>
#include <queue>
#include <chrono>

#include "Graph.h"
#include "MD_Tree.h"
#include "Util.h"
#include "TreeList.h"


using namespace std;

MD_Tree getModularDecomposition(const Graph& graph);

/**
 * Returns the contents of a file as a string.
 *
 * @param filename The name of the file to be read.
 * @return A string containing the contents of the file, or an empty string if an error occurred.
 */
string readFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return "";
    }
    ostringstream oss;
    oss << file.rdbuf();
    file.close();

    return oss.str();
}


/**
* Normally, the modular decomposition for the subgraphs would be performed recursively.
* However, as long as the algorithm does not work completely, this helper function returns
* the hard-coded "recursively" computed MD_Trees for a specific example. This allows the
* algorithm to be tested, even though it is not completely finished.
*
* @param subgraph The subgraph of which the MD_Tree should be recursively computed
* @return the hard-coded MD_Tree for specific subgraphs
*/
MD_Tree getRecursiveComputation(const Graph& subgraph)
{
    if (subgraph.getAdjlist().size() == 4) {

        TreeNode* nodeA = new TreeNode(0);
        TreeNode* nodeC = new TreeNode(1);
        TreeNode* nodeD = new TreeNode(2);
        TreeNode* nodeE = new TreeNode(3);

        TreeNode* node1 = new TreeNode(SERIES);
        TreeNode* node0 = new TreeNode(PARALLEL);
        TreeNode* nodeRoot = new TreeNode(SERIES);

        setChild(nodeRoot, nodeA);
        setChild(node0, nodeC);
        setChild(node1, nodeD);

        setSibling(nodeC, node1);
        setSibling(nodeA, node0);
        setSibling(nodeD, nodeE);

        return MD_Tree(nodeRoot);
    }

    if (subgraph.getAdjlist().size() == 12) {
        TreeNode* nodeB = new TreeNode(0);
        TreeNode* nodeJ = new TreeNode(1);
        TreeNode* nodeF = new TreeNode(2);
        TreeNode* nodeG = new TreeNode(3);
        TreeNode* nodeH = new TreeNode(4);
        TreeNode* nodeI = new TreeNode(5);
        TreeNode* nodeK = new TreeNode(6);
        TreeNode* nodeL = new TreeNode(7);
        TreeNode* nodeM = new TreeNode(8);
        TreeNode* nodeN = new TreeNode(9);
        TreeNode* nodeP = new TreeNode(10);
        TreeNode* nodeQ = new TreeNode(11);

        TreeNode* nodeRoot = new TreeNode(PARALLEL);
        TreeNode* node1Left = new TreeNode(SERIES);
        TreeNode* node1Mid = new TreeNode(SERIES);
        TreeNode* node1Right = new TreeNode(SERIES);
        TreeNode* node0 = new TreeNode(PARALLEL);

        setChild(nodeRoot, node1Left);
        setChild(node1Mid, nodeK);
        setChild(node1Right, nodeQ);
        setChild(node0, nodeN);

        setChild(node1Left, nodeI);
        setSibling(nodeI, nodeG);
        setSibling(nodeG, nodeH);
        setSibling(nodeH, nodeF);

        setSibling(nodeK, nodeL);
        setSibling(nodeQ, nodeM);
        setSibling(nodeN, nodeP);
        setSibling(node1Left, nodeB);
        setSibling(nodeB, node1Mid);
        setSibling(node1Mid, nodeJ);
        setSibling(nodeJ, node1Right);
        setSibling(nodeM, node0);

        return MD_Tree(nodeRoot);
    }

    TreeNode* root = new TreeNode(0);
    return MD_Tree(root);
}

/**
* Prints a given modular - decomposition forest on the command line
*
* @param current The forest to be printed, passed by a vector of its trees
*/
void printForest(const vector<MD_Tree>& forest) {
    for (int i = 0; i < forest.size(); i++) {
        cout << "MD Tree " << (i + 1) << ": " << endl;
        printTree(forest[i].root);
        cout << endl;
    }
}


/**
* Marks a node and all its ancestor with either "left" or "right", based on the given parameter
*
* @param node The node that should be marked (along with its ancestors)
* @param markLeft If the node should be marked with "left" ("right" otherwise)
*/
void markNodeAndAncestors(TreeNode* node, bool markLeft) {
    TreeNode* ancestor = node;
    if (markLeft) {
        while (ancestor != nullptr && !ancestor->markedLeft) {
            ancestor->markedLeft = true;
            ancestor = ancestor->parent;
        }
    }
    else {
        while (ancestor != nullptr && !ancestor->markedRight) {
            ancestor->markedRight = true;
            ancestor = ancestor->parent;
        }
    }
}

/**
* Marks all children of a given node with either "left" or "right", based on the given parameter
*
* @param node The node which's children should be marked (along with its ancestors)
* @param markLeft If the node's children should be marked with "left" ("right" otherwise)
*/
void markChildren(TreeNode* node, bool markLeft) {
    TreeNode* next = node->child;
    while (next != nullptr) {
        if (markLeft) {
            next->markedLeft = true;
        }
        else {
            next->markedRight = true;
        }
        next = next->sibling;
    }
}

/**
* Constructs an MD_Tree, where the root has a given label and the given children
*
* @param X A set containing all children that the root of the new tree should have
* @param label The label of the new tree's root
*/
MD_Tree constructTree(vector<TreeNode*> X, Label label) {
    TreeNode* T;
    if (X.size() == 1) {
        T = X[0];
        T->sibling = nullptr;
        T->parent = nullptr;
    }
    else {
        T = new TreeNode(label);
        for (int i = 0; i < X.size(); i++) {
            if (i == 0) {
                setChild(T, X[i]);
            }
            else {
                setSibling(X[i - 1], X[i]);
            }
        }
        X[X.size() - 1]->sibling = nullptr;
    }
    return MD_Tree(T);
}

/**
* Refines a MD_Forest with a node that is not prime. For details on the refinement process
* see the algorithm description.
*
* @param forest The forest that should be refined, passed by its first element.
* @param p The non-Prime tree node on base of which the forest should be refined.
* @param maxSubTrees A set containing all maximal Subtrees.
* @param isLeftSplit If a left-split should be used (right-split otherwise).
*/
void refineByNonPrimeNode(TreeList& forest, TreeNode* p,
    const unordered_set<TreeNode*>& maxSubTrees, bool isLeftSplit) {

    vector<TreeNode*> A;
    vector<TreeNode*> B;

    TreeNode* next = p->child;
    while (next != nullptr) {
        if (maxSubTrees.find(next) != maxSubTrees.end()) {
            A.push_back(next);
        }
        else {
            B.push_back(next);
        }
        next = next->sibling;
    }

    if (A.size() > 0 && B.size() > 0) {
        MD_Tree* Ta = new MD_Tree(constructTree(A, p->label));
        MD_Tree* Tb = new MD_Tree(constructTree(B, p->label));

        if (p->parent == nullptr) {
            MD_Tree* currentTree = forest.getCorrespondingTree(p);

            if (isLeftSplit) {
                forest.replaceElement(currentTree, Ta, Tb);
            }
            else {
                forest.replaceElement(currentTree, Tb, Ta);
            }
        }
        else {
            setChild(p, Ta->root);
            setSibling(Ta->root, Tb->root);
        }
        markNodeAndAncestors(Ta->root, isLeftSplit);
        markNodeAndAncestors(Tb->root, isLeftSplit);
    }
}


/**
* Refines a MD_Forest with a node that is prime. For details on the refinement process
* see the algorithm description.
*
* @param forest The forest that should be refined.
* @param isLeftSplit If a left-split should be used (right-split otherwise).
*/
void refineByPrimeNode(TreeNode* node, bool isLeftSplit) {
    markNodeAndAncestors(node, isLeftSplit);
    markChildren(node, isLeftSplit);
}

/**
* Refines a MD_Forest by a given set of active edges. For details on the refinement process see
* the algorithm description.
*
* @param forest The MD_Forest to be refined, passed by its first element
* @param X The set of active edges that should be used for refinement.
* @param nodeIsLeft If the node that was used for calculating the active edges can be found
*   on the left side of the pivot element.
*/
void refineBySet(TreeList& forest, unordered_set<int>& X, int timestemp, bool nodeIsLeft) {
    unordered_set<TreeNode*> maxSubtrees;
    unordered_set<TreeNode*> maxSubtreeParents;

    MD_Tree* currentTree = forest.getStart();
    while (currentTree != nullptr) {
        unordered_set<TreeNode*> subtrees;
        getMaxContSubTrees(currentTree->root, currentTree->root, X, subtrees, timestemp);
        maxSubtrees.insert(subtrees.begin(), subtrees.end());

        currentTree = currentTree->right;
    }

    for (TreeNode* node : maxSubtrees) {
        if (node->parent == nullptr) {
            maxSubtreeParents.insert(node);
        }
        else {
            maxSubtreeParents.insert(node->parent);
        }
    }

    for (TreeNode* parentNode : maxSubtreeParents) {
        MD_Tree* correspondingTree = forest.getCorrespondingTree(parentNode);
        bool leftSplit = nodeIsLeft || (correspondingTree->left == nullptr);
        if (parentNode->label == PRIME) {
            refineByPrimeNode(parentNode, leftSplit);
        }
        else {
            refineByNonPrimeNode(forest, parentNode, maxSubtrees, leftSplit);
        }
    }
}

/**
* Executes the promotion algorithm for a specific tree, given by its root. For more details on promotion,
* see the algorithm description.
*
* @param root The root of the MD_Tree that should be promoted
* @return A list of all trees that are calculated in the promotion
*/
vector<MD_Tree> getPromotedTree(TreeNode* root) {
    vector<MD_Tree> forest;
    if (root->markedLeft) {
        TreeNode* markedChild = root->child;
        TreeNode* previous = root;

        while (markedChild != nullptr) {
            if (markedChild->markedLeft) {
                if (previous == root) {
                    root->child = markedChild->sibling;
                }
                else {
                    previous->sibling = markedChild->sibling;
                }
                markedChild->parent = nullptr;
                markedChild->sibling = nullptr;
                vector<MD_Tree> left = getPromotedTree(markedChild);
                forest.insert(forest.end(), left.begin(), left.end());
            }
            previous = markedChild;
            markedChild = markedChild->sibling;
        }
    }
    forest.push_back(MD_Tree(root));
    if (root->markedRight) {
        TreeNode* markedChild = root->child;
        TreeNode* previous = root;

        while (markedChild != nullptr) {
            if (markedChild->markedRight) {
                if (previous == root) {
                    root->child = markedChild->sibling;
                }
                else {
                    previous->sibling = markedChild->sibling;
                }
                markedChild->parent = nullptr;
                markedChild->sibling = nullptr;
                vector<MD_Tree> right = getPromotedTree(markedChild);
                forest.insert(forest.end(), right.begin(), right.end());
            }
            previous = markedChild;
            markedChild = markedChild->sibling;
        }
    }
    return forest;
}

/**
* Deletes all marks ("left" or "right"), starting at a given node.
*
* @param node The starting node for the mark-deletion process.
*/
void deleteMarks(TreeNode* node) {
    node->markedLeft = false;
    node->markedRight = false;

    if (node->sibling != nullptr) {
        deleteMarks(node->sibling);
    }
    if (node->child != nullptr) {
        deleteMarks(node->child);
    }
}

/**
* Cleans the forest after promotion is done. This includes:
*   - Every root with no child will be removed (Except the root is a leaf).
*   - Whenever a root as only one child, that child will take the place of the root.
*   - All markings will be delted.
*
* @param forest The forest to be cleaned up.
*/
void cleanUp(vector<MD_Tree>& forest) {
    vector<MD_Tree> newTreeList;
    for (MD_Tree& tree : forest) {
        if (tree.root->markedLeft || tree.root->markedRight) {
            if (tree.root->child != nullptr) {
                if (tree.root->child->sibling == nullptr) {
                    newTreeList.push_back(MD_Tree(tree.root->child));
                }
                else {
                    newTreeList.push_back(MD_Tree(tree.root));
                }
            }
            else if (tree.root->label == LEAF) {
                newTreeList.push_back(MD_Tree(tree.root));
            }
        }
        else {
            newTreeList.push_back(MD_Tree(tree.root));
        }
    }
    for (MD_Tree& tree : newTreeList) {
        deleteMarks(tree.root);
    }
    forest = newTreeList;
}

/**
* Updates the value of a given node based on a given list. Calls itself recursively.
*
* @param node The current node.
* @param updatedValues A vector containing the new values of the nodes.
*/
void updateTreeValues(TreeNode* node, const vector<int>& updatedValues) {
    if (node->label == LEAF && node->value < updatedValues.size()) {
        node->value = updatedValues[node->value];
    }
    if (node->sibling != nullptr) {
        updateTreeValues(node->sibling, updatedValues);
    }
    if (node->child != nullptr) {
        updateTreeValues(node->child, updatedValues);
    }
}


/**
* Inserts the given index in the list of indices, if a given node has a value.
* Afterwards, this function recursively calls the nodes children and siblings.
* This is a helper function for getMaxModuleIndices.
*
* @param node The current node.
* @param indices The list to insert into.
* @param index The value to insert.
*/
void insertIndex(TreeNode* node, vector<int>& indices, int index) {
    if (node->label == LEAF) {
        indices[node->value] = index;
    }
    if (node->sibling != nullptr) {
        insertIndex(node->sibling, indices, index);
    }
    if (node->child != nullptr) {
        insertIndex(node->child, indices, index);
    }
}

/**
* Computes a list that provides information about every vertex in the graph:
* The index of the module it belongs to at the start of the assembly-process.
*
* @param forest The MD-forest
* @param graphSize The number of elements in the graph.
* @return A list that contains the searched information.
*/
vector<int> getMaxModuleIndices(const vector<MD_Tree>& forest, int graphSize) {
    vector<int> indices(graphSize, -1);
    for (int i = 0; i < forest.size(); i++) {
        insertIndex(forest[i].root, indices, i);
    }
    return indices;
}


/**
* Creates a large list of adajencies that contains all adjacencies of the elements in a given list.
*
* @param vertices The vertices to sum up the adjacencies of.
* @param graph The graph.
* @return The summed-up list of adjacencies.
*/
vector<int> getTotalAdjacencyList(const vector<int>& vertices, const Graph& graph) {
    vector<vector<int>> adjlist = graph.getAdjlist();
    vector<int> totalAdjacencies;
    for (int vertex : vertices) {
        totalAdjacencies.insert(totalAdjacencies.end(), adjlist[vertex].begin(), adjlist[vertex].end());
    }
    return totalAdjacencies;
}
// Gleich die Benachabarten maxModules eintragen

/**
* Uses the connections given in the graph to insert left- and right pointers for every element in
* the forest. The left pointer of an element X is on the lowest index, so that all elements with lower
* indices are connected to X. The right pointer of this element is on the highest index, so that all
* elements with higher indices are disconnected to X.
*
* @param forest Contains the list of MD-trees.
* @param graph The graph.
*/
void insertLeftRightPointers(vector<MD_Tree>& forest, const Graph& graph) {
    int n = forest.size();
    vector<int> maxModuleIndices = getMaxModuleIndices(forest, graph.getAdjlist().size());
    vector<bool> connections(n, false);

    for (int i = 0; i < n; i++) {
        vector<int> vertices = getPreOrderLeafs(forest[i].root);
        vector<int> totalAdjacencies = getTotalAdjacencyList(vertices, graph);
        int maxConnectionIndex = numeric_limits<int>::min();
        vector<int> connectionIndices;

        for (int j = 0; j < totalAdjacencies.size(); j++) {
            int connectedModule = maxModuleIndices[totalAdjacencies[j]];
            if (connectedModule != -1) {
                connections[connectedModule] = true;
                connectionIndices.push_back(connectedModule);
                if (connectedModule > maxConnectionIndex) {
                    maxConnectionIndex = connectedModule;
                }
            }
        }
        int leftPointer = 0;
        while (leftPointer < i && connections[leftPointer]) {
            leftPointer++;
        }

        // Left and Right pointers are always to the left of the element with the given index
        forest[i].leftIndex = leftPointer;
        forest[i].rightIndex = maxConnectionIndex + 1;

        // Reset the connections vector
        for (int connection : connectionIndices) {
            connections[connection] = false;
        }
    }
}

/**
* Checks, if the subtree induced by a given node is connected to the pivot element.
*
* @param node The node to check.
* @param pivotAdj The adjacency list of the pivot element.
* @return If there is a connection to the pivot.
*/
bool isConnectedToPivot(TreeNode* node, const vector<int>& pivotAdj) {
    if (node->label == LEAF) {
        return find(pivotAdj.begin(), pivotAdj.end(), node->value) != pivotAdj.end();
    }
    if (node->child != nullptr) {
        return isConnectedToPivot(node->child, pivotAdj);
    }
    return isConnectedToPivot(node->sibling, pivotAdj);
}

/**
* Checks, if the next module in a list of trees only extends to the right, meaning that it is parallel.
*
* @param forest Contains the list of trees.
* @param currentLeft The left index of the current module.
* @param currentRight The right index of the current module.
* @return if the next module is parallel.
*/
bool checkForParallel(vector<MD_Tree>& forest, int currentLeft, int currentRight) {
    if (currentRight >= forest.size()) {
        return false;
    }
    int i = currentRight;
    currentRight++;
    while (i < currentRight) {
        int leftPointer = forest[i].leftIndex;
        int rightPointer = forest[i].rightIndex;

        if (leftPointer < currentLeft) {
            return false;
        }
        if (rightPointer > currentRight) {
            currentRight = rightPointer;
        }
        i++;
    }
    return true;
}

/**
* A helper method for the construction of the modular decompositon
* tree in the assembly step of the algorithm.
*
* @param currentNode The node to be added to the MD tree.
* @param currentModuleType The type of the current module
* @param lastNode The last node that was added to the tree.
*/
void addToMDTree(TreeNode*& currentNode, Label& currentModuleType, TreeNode*& lastNode) {
    if (currentNode->label == currentModuleType && currentNode->child != nullptr) {
        setSibling(lastNode, currentNode->child);
        lastNode = currentNode->child;
        while (lastNode->sibling != nullptr) {
            setSibling(lastNode, lastNode->sibling);
            lastNode = lastNode->sibling;
        }
    }
    else {
        setSibling(lastNode, currentNode);
        lastNode = currentNode;
    }
}


/**
* Performs the first step of the algorithm, the recursion. In this step, the vertices of
* the given graph are sorted by their distance to a given pivot element. The MD_Trees of each
* set are computed recursively and stored in a MD_Forest, that contains a list of theses MD_Trees.
* Furthermore, the active edges are computed as well.
*
* @param graph The graph on which the modular decomposition should be executed
* @param pivot The arbitrarly choosen vertex that works as pivot-element for the algorithm
* @param activeEdges This vector will contain the activeEdges of the given graph with the given pivot,
*   after the method has terminated
* @return The resulting MD_Forest (list of MD_Trees)
*/
TreeList recursion(const Graph& graph, int pivot, vector<unordered_set<int>>& activeEdges, vector<bool>& leftNodes) {
    int currentIndex = 0;
    vector<vector<int>> adjlist = graph.getAdjlist();
    activeEdges.assign(adjlist.size(), unordered_set<int>());
    vector<unordered_set<int>> N;
    unordered_set<int> remainingNodes;
    TreeList* output = new TreeList();

    do {
        unordered_set<int> currentSet;
        if (currentIndex == 0) {
            for (int i = 0; i < adjlist.size(); i++) {
                if (i != pivot) {
                    if (find(adjlist[pivot].begin(), adjlist[pivot].end(), i) != adjlist[pivot].end()) {
                        currentSet.insert(i);
                        activeEdges[i].insert(pivot);
                        activeEdges[pivot].insert(i);
                        leftNodes.push_back(true);
                    }
                    else {
                        remainingNodes.insert(i);
                        leftNodes.push_back(false);
                    }
                }
                else {
                    leftNodes.push_back(false);
                }
            }
        }
        else {
            for (int node : remainingNodes) {
                unordered_set<int> lastN = N[currentIndex - 1];
                for (int adj : adjlist[node]) {
                    if (lastN.find(adj) != lastN.end()) {
                        currentSet.insert(node);
                        activeEdges[adj].insert(node);
                        activeEdges[node].insert(adj);
                    }
                }
            }
            for (int node : currentSet) {
                remainingNodes.erase(node);
            }

            vector<int> subgraphIndexMapping;
            Graph subgraph = graph.getSubGraph(N[currentIndex - 1], subgraphIndexMapping);
            MD_Tree* tree = new MD_Tree(getModularDecomposition(subgraph));
            updateTreeValues(tree->root, subgraphIndexMapping);
            output->insert(tree);
        }
        N.push_back(currentSet);
        currentIndex++;
    } while (remainingNodes.size() > 0);

    vector<int> subgraphIndexMapping;
    Graph subgraph = graph.getSubGraph(N[currentIndex - 1], subgraphIndexMapping);
    MD_Tree* tree = new MD_Tree(getModularDecomposition(subgraph));
    updateTreeValues(tree->root, subgraphIndexMapping);
    output->insert(tree);

    return *output;
}

/**
* Performs the second step of the algorithm, the refinement. As step 1 (recursion) only used one node (the pivot) to
* separate all vertices into modules, refinement takes care of using all other nodes to separate the vertices even further.
* For more details on the refinement process, see the algorithm description.
*
* @param graph The graph that should be modular decomposed.
* @param previousPivot The pivot element that was selected in step 1 (recursion).
* @param forest The MD_Forest that was calculated in step 1 (recursion).
* @param activeEdges A set of all activeEdges.
*/
void refinement(const Graph& graph, int previousPivot, TreeList& forest,
    const vector<unordered_set<int>>& activeEdges, const vector<bool>& leftNodes) {

    for (int i = 0; i < graph.getAdjlist().size(); i++) {
        if (i != previousPivot) {
            unordered_set<int> incidentActives = activeEdges[i];
            refineBySet(forest, incidentActives, i, leftNodes[i]);
        }
    }
}

/**
* Executes the third step of the algorithm, the promotion. For more details on the promotion - process,
* see the algorithm description.
*
* @param forest The MD_Forest to execute promotion on.
* @return The promoted forest as a vector of MD_Trees
*/
vector<MD_Tree> promotion(TreeList& forest) {
    vector<MD_Tree> newTreeList;
    MD_Tree* current = forest.getStart();
    while (current != nullptr) {
        vector<MD_Tree> promotedList = getPromotedTree(current->root);
        newTreeList.insert(newTreeList.end(), promotedList.begin(), promotedList.end());
        current = current->right;
    }
    cleanUp(newTreeList);
    return newTreeList;
}

/**
* Executes the fourth (and final) step of the algorithm, the assembly. For more details on the
* assembly - process, see the algorithm description.
*
* @param forest The MD_Forest resulting from step 3: promotion.
* @param graph The graph.
* @param pivot The previously chosen pivot element.
* @return The final assembled MD-tree.
*/
MD_Tree assembly(vector<MD_Tree>& forest, const Graph& graph, int pivot) {

    MD_Tree pivotTree = MD_Tree(new TreeNode(pivot));
    int pivotIndex = 1;
    while (pivotIndex < forest.size() && isConnectedToPivot(forest[pivotIndex].root, 
        graph.getAdjlist()[pivot])) {
        pivotIndex++;
    }
    forest.insert(forest.begin() + pivotIndex, pivotTree);
    insertLeftRightPointers(forest, graph);

    TreeNode* lastModule = new TreeNode(pivot);

    int currentLeft = pivotIndex;
    int currentRight = pivotIndex + 1;
    int includedLeft = pivotIndex;
    int includedRight = pivotIndex + 1;

    do {
        bool addedRight = false;
        bool addedLeft = false;
        queue<int> maxModuleIndices;

        if (checkForParallel(forest, currentLeft, currentRight)) {
            maxModuleIndices.push(currentRight);
            addedRight = true;
            currentRight++;
        }
        else {
            maxModuleIndices.push(currentLeft - 1);
            addedLeft = true;
            currentLeft--;
        }

        do {
            int currentMaxModule = maxModuleIndices.front();
            maxModuleIndices.pop();

            int leftPointer = forest[currentMaxModule].leftIndex;
            int rightPointer = forest[currentMaxModule].rightIndex;

            if (leftPointer < currentLeft) {
                currentLeft = leftPointer;
                maxModuleIndices.push(leftPointer);
                addedLeft = true;
            }
            if (rightPointer > currentRight) {
                currentRight = rightPointer;
                maxModuleIndices.push(rightPointer - 1);
                addedRight = true;
            }
        } while (!maxModuleIndices.empty());

        Label moduleType = PARALLEL;
        if (addedLeft && addedRight) {
            moduleType = PRIME;
        }
        else if (addedLeft) {
            moduleType = SERIES;
        }
        TreeNode* moduleNode = new TreeNode(moduleType);
        setChild(moduleNode, lastModule);
        TreeNode* lastNode = lastModule;

        for (int i = currentLeft; i < includedLeft; i++) {
            addToMDTree(forest[i].root, moduleType, lastNode);
        }
        for (int i = includedRight; i < currentRight; i++) {
            addToMDTree(forest[i].root, moduleType, lastNode);
        }

        includedLeft = currentLeft;
        includedRight = currentRight;
        lastModule = moduleNode;

    } while (currentLeft > 0 || currentRight < forest.size());

    return MD_Tree(lastModule);
}

/**
* Returns the Modular Decomposition for a disconnected Graph, by creating one PARALLEL node
* and setting the recursively computed MD_Trees of the graphs components as its children.
*
* @param graph The graph
* @return The recursively computed MD_Tree
*/
MD_Tree getModularDecompositionDisconnectedGraph(const Graph& graph) {
    vector<unordered_set<int>> components = graph.getConnectedComponents();
    TreeNode* rootNode = new TreeNode(PARALLEL);
    TreeNode* lastChild = nullptr;

    for (int i = 0; i < components.size(); i++) {
        vector<int> subgraphIndexMapping;
        Graph subgraph = graph.getSubGraph(components[i], subgraphIndexMapping);
        MD_Tree tree = getModularDecomposition(subgraph);
        updateTreeValues(tree.root, subgraphIndexMapping);
        if (i == 0) {
            setChild(rootNode, tree.root);
            lastChild = tree.root;
        }
        else {
            setSibling(lastChild, tree.root);
            lastChild = tree.root;
        }
    }
    return MD_Tree(rootNode);
}

/**
* Returns the modular decomposition tree for a given graph. This is done by using the four steps:
*   - Recursion
*   - Refinement
*   - Promotion
*   - Assembly
* For more details on all of these steps and the algorithm itself, see the algorithm description.
*
* @param graph The graph to be modular decomposed.
*/
MD_Tree getModularDecomposition(const Graph& graph) {

    if (graph.getAdjlist().size() == 0) {
        cerr << "You cannot get a modular decomposition of a graph with no vertices!" << endl;
    }
    else if (graph.getAdjlist().size() == 1) {
        TreeNode* root = new TreeNode(0);
        // Create a vector that maps value -> TreeNode*
        return MD_Tree(root);
    }

    if (graph.isConnected()) {
        int pivot = graph.getAdjlist().size() >= 14 ? 14 : 0;
        vector<unordered_set<int>> activeEdges;
        vector<bool> leftNodes;

        TreeList forest = recursion(graph, pivot, activeEdges, leftNodes);

        /*cout << "After Recursion: " << endl;
        printForest(&forest);
        cout << endl;*/

        refinement(graph, pivot, forest, activeEdges, leftNodes);

        /*cout << "After Refinement: " << endl;
        printForest(&forest);
        cout << endl;*/

        vector<MD_Tree> forestVec = promotion(forest);

        /*cout << "After Promotion: " << endl;
        printForest(forestVec);
        cout << endl;*/

        MD_Tree finalResult = assembly(forestVec, graph, pivot);
        resetTimestemps(finalResult.root);

        return finalResult;
    }
    else {
        MD_Tree finalTree = getModularDecompositionDisconnectedGraph(graph);
        return finalTree;
    }
}

/**
* The main method. Doesn't do that much.
*/
int main() {
    string filePath;
    cout << "Input the path to the file containing the graph as adjacency list: ";
    cin >> filePath;
    string adjList = readFile(filePath);
    vector<int> indexMapping;

    if (!Util::isConsecutivelyOrdered(adjList)) {
        vector<int> indexMapping = Util::rewriteAdjacencyList(adjList);
    }

    Graph graph = Graph(adjList);

    auto start = chrono::high_resolution_clock::now();
    MD_Tree mdTree = getModularDecomposition(graph);
    auto end = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

    Util::sortTree(mdTree);
    if (indexMapping.size() > 0) {
        updateTreeValues(mdTree.root, indexMapping);
    }
    cout << "The final MD-tree: " << endl << endl;
    printTree(mdTree.root);
    if (Util::testModularDecompositionTree(graph, mdTree)) {
        cout << endl << "This result appears to be correct" << endl;
    }
    else {
        cout << endl << "Unfortunately, this result appears to be incorrect!" << endl;
    }
    cout << endl << "Time needed: " << duration.count() << "microseconds" << endl;
}

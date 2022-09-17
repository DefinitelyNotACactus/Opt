//
//  BranchBound.cpp
//  Opt
//
//  Created by David Galvao on 28/08/22.
//

#define TIMEOUT 60

#include "BranchBound.hpp"

#include <iostream>
#include <queue>
#include <stack>

std::vector<std::vector<int>> getSubtours(hungarian_problem_t &p) {
    // Obter os subtours da solução do algoritmo hungaro
    std::vector<std::vector<int>> subtours;
    std::vector<bool> visited(p.num_rows, false);
    for (int i = 0; i < p.num_rows; i++) {
        if (!visited[i]) {
            int next = i;
            std::vector<int> subtour;
            subtour.push_back(i);
            do {
                for (int j = 0; j < p.num_cols; j++) {
                    if(p.assignment[next][j]) {
                        next = j;
                        visited[j] = true;
                        subtour.push_back(j);
                        break;
                    }
                }
            } while(next != i);
            visited[i] = true;
            subtours.push_back(subtour);
        }
    }
    return subtours;
}

void printSubtours(const std::vector<std::vector<int>> &subtours) {
    int i = 0;
    for (std::vector<int> subtour : subtours) {
        std::cout << "Subtour " << i++ << ":";
        printSubtour(subtour);
    }
}

void printSubtour(const std::vector<int> &subtour) {
    for(int v : subtour) {
        std::cout << " " << v;
    }
    std::cout << "\n";
}

Tree buildTree(Data &data, double upperBound, int strategy, bool lagragean) {
    int dimension = data.getDimension();
    Tree tree;
    tree.start = std::chrono::steady_clock::now();
    Node root;
    // Computar um upper bound, utilizando o vizinho mais próximo
    if (upperBound < 0) {
        nearestNeighbor(tree.solution, data.getMatrixCost(), dimension);
        tree.upperBound = tree.solution.cost;
        rvnd(data, tree.solution.tour, tree.upperBound);
    } else {
        tree.upperBound = upperBound;
    }
    std::cout << "Starting UB: " << tree.upperBound << "\n";
    // Chamar o algoritmo hungaro na raiz
    if (!lagragean) {
        computeSolution(root, data, dimension, tree.upperBound);
    } else {
        Kruskal kruskal(data, tree.upperBound, 1);
        root.multipliers = kruskal.getMultipliers();
        root.multipliersSum = kruskal.getMultipliersSum();
        root.forbiddenEdges = kruskal.getNewForbiddenEdges();
        root.lowerBound = kruskal.getCost();
    }
    // Percorrer a árvore de acordo com a estratégia escolhida
    switch (strategy) {
        case 0:
            traverseTreeDFS(data, tree, root, lagragean);
            break;
        case 1:
            traverseTreeBFS(data, tree, root, lagragean);
            break;
        case 2:
            traverseTreeBestBound(data, tree, root, lagragean);
            break;
        default:
            traverseTreeDFS(data, tree, root, lagragean);
    }
    
    return tree;
}

void traverseTreeBestBound(Data &data, Tree &tree, Node &root, bool lagragean) {
    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> nodes;
    nodes.push(root);
    while (!nodes.empty()) {
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(now - tree.start).count() > TIMEOUT) {
            std::cout << "Timeout\n";
            break;
        }
        Node best = nodes.top();
        if (best.lowerBound > tree.upperBound) {
            break;
        }
        nodes.pop();
        for (int i = 0; i < best.subtours[best.chosen].size() - 1; i++) {
            std::cout << "Tree Size: " << nodes.size() << " | ";
            Node n;
            n.forbiddenEdges = best.forbiddenEdges;
            std::pair<int, int> forbiddenEdge;
            forbiddenEdge.first = best.subtours[best.chosen][i];
            forbiddenEdge.second = best.subtours[best.chosen][i + 1];
            n.forbiddenEdges.push_back(forbiddenEdge);
            //n.forbiddenEdgesWeights.push_back(data.getDistance(forbiddenEdge.first, forbiddenEdge.second));
            computeSolution(n, data, data.getDimension(), tree.upperBound);
            if (n.subtours.size() == 1) { // Encontramos uma solução viável
                n.prune = true;
                // Tentar melhorar uma solução viável por meio de movimentos de vizinhança
                rvnd(data, n.subtours[0], n.lowerBound);
                if (n.lowerBound < tree.upperBound) {
                    tree.upperBound = n.lowerBound;
                    std::cout << "New UB: " << tree.upperBound << "\n";
                    if(!nodes.empty()) {
                        std::priority_queue<Node, std::vector<Node>, std::greater<Node>> nodes_aux; // Podar a árvore
                        nodes_aux.push(nodes.top());
                        nodes.pop();
                        while(!nodes.empty()) {
                            Node s = nodes.top();
                            nodes.pop();
                            if (s.lowerBound < tree.upperBound) {
                                nodes_aux.push(s);
                            }
                        }
                        nodes = nodes_aux;
                    }
                } else {
                    std::cout << "Encontrada solução viável (Custo = " << n.lowerBound << ")\n";
                }
            }
            if (!n.prune) {
                nodes.push(n);
            }
        }
    }
}

void traverseTreeBFS(Data &data, Tree &tree, Node &root, bool lagragean) {
    std::queue<Node> nodes;
    nodes.push(root);
    while (!nodes.empty()) {
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(now - tree.start).count() > TIMEOUT) {
            std::cout << "Timeout\n";
            break;
        }
        Node front = nodes.front();
        nodes.pop();
        if (front.lowerBound > tree.upperBound) {
            continue;
        }
        for (int i = 0; i < front.subtours[front.chosen].size() - 1; i++) {
            std::cout << "Tree Size: " << nodes.size() << " | ";
            Node n;
            n.forbiddenEdges = front.forbiddenEdges;
            std::pair<int, int> forbiddenEdge;
            forbiddenEdge.first = front.subtours[front.chosen][i];
            forbiddenEdge.second = front.subtours[front.chosen][i + 1];
            n.forbiddenEdges.push_back(forbiddenEdge);
            //n.forbiddenEdgesWeights.push_back(data.getDistance(forbiddenEdge.first, forbiddenEdge.second));
            computeSolution(n, data, data.getDimension(), tree.upperBound);
            if (n.subtours.size() == 1) { // Encontramos uma solução viável
                n.prune = true;
                // Tentar melhorar uma solução viável por meio de movimentos de vizinhança
                rvnd(data, n.subtours[0], n.lowerBound);
                if (n.lowerBound < tree.upperBound) {
                    tree.upperBound = n.lowerBound;
                    std::cout << "New UB: " << tree.upperBound << "\n";
                    if(!nodes.empty()) {
                        std::queue<Node> nodes_aux; // Podar a árvore
                        nodes_aux.push(nodes.front());
                        nodes.pop();
                        while(!nodes.empty()) {
                            Node s = nodes.front();
                            nodes.pop();
                            if (s.lowerBound < tree.upperBound) {
                                nodes_aux.push(s);
                            }
                        }
                        nodes = nodes_aux;
                    }
                } else {
                    std::cout << "Encontrada solução viável (Custo = " << n.lowerBound << ")\n";
                }
            }
            if (!n.prune) {
                nodes.push(n);
            }
        }
    }
}

void traverseTreeDFS(Data &data, Tree &tree, Node &root, bool lagragean) {
    std::stack<Node> nodes;
    nodes.push(root);
    while (!nodes.empty()) {
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(now - tree.start).count() > TIMEOUT) {
            std::cout << "Timeout\n";
            break;
        }
        Node top = nodes.top();
        nodes.pop();
        if (top.lowerBound > tree.upperBound) {
            continue;
        }
        for (int i = 0; i < top.subtours[top.chosen].size() - 1; i++) {
            std::cout << "Tree Size: " << nodes.size() << " | ";
            Node n;
            n.forbiddenEdges = top.forbiddenEdges;
            std::pair<int, int> forbiddenEdge;
            forbiddenEdge.first = top.subtours[top.chosen][i];
            forbiddenEdge.second = top.subtours[top.chosen][i + 1];
            n.forbiddenEdges.push_back(forbiddenEdge);
            //n.forbiddenEdgesWeights.push_back(data.getDistance(forbiddenEdge.first, forbiddenEdge.second));
            computeSolution(n, data, data.getDimension(), tree.upperBound);
            if (n.subtours.size() == 1) { // Encontramos uma solução viável
                n.prune = true;
                // Tentar melhorar uma solução viável por meio de movimentos de vizinhança
                rvnd(data, n.subtours[0], n.lowerBound);
                if (n.lowerBound < tree.upperBound) {
                    tree.upperBound = n.lowerBound;
                    std::cout << "New UB: " << tree.upperBound << "\n";
                    if(!nodes.empty()) {
                        std::stack<Node> nodes_aux; // Podar a árvore
                        nodes_aux.push(nodes.top());
                        nodes.pop();
                        while(!nodes.empty()) {
                            Node s = nodes.top();
                            nodes.pop();
                            if (s.lowerBound < tree.upperBound) {
                                nodes_aux.push(s);
                            }
                        }
                        nodes = nodes_aux;
                    }
                } else {
                    std::cout << "Encontrada solução viável (Custo = " << n.lowerBound << ")\n";
                }
            }
            if (!n.prune) {
                nodes.push(n);
            }
        }
    }
}

void computeSolution(Node &node, Data &data, int dimension, double upperBound) {
    double **matrix = new double*[dimension];
    for (int i = 0; i < dimension; i++){
        matrix[i] = new double[dimension];
        for (int j = 0; j < dimension; j++){
            matrix[i][j] = data.getDistance(i,j);
        }
    }
    // Definir as arestas proibidas como custo infinito na matriz
    for (std::pair<int, int> edge : node.forbiddenEdges) {
        matrix[edge.first][edge.second] = INFINITE;
    }
    // Resolver o problema utilizando o algoritmo hungaro
    hungarian_problem_t p;
    hungarian_init(&p, matrix, dimension, dimension, HUNGARIAN_MODE_MINIMIZE_COST); // Carregando o problema

    double objValue = hungarian_solve(&p);
    std::cout << "LB: " << objValue << " UB: " << upperBound << "\n";
    if (objValue > upperBound) {
        node.prune = true;
    } else {
        node.prune = false;
        node.lowerBound = objValue;
        node.subtours = getSubtours(p);
        node.chosen = chooseSubtour(node.subtours);
    }
    // Desfazer os valores das arestas proibidas
//    for (int i = 0; i < node.forbiddenEdges.size(); i++) {
//        data.getMatrixCost()[node.forbiddenEdges[i].first][node.forbiddenEdges[i].second] = node.forbiddenEdgesWeights[i];
//    }
    // Desalocar
    hungarian_free(&p);
    for (int i = 0; i < dimension; i++) delete [] matrix[i];
    delete [] matrix;
}

std::vector<std::pair<int, int>> computeSolutionLagrangean(Node &node, Data &data, int dimension, double upperBound) {
    Kruskal kruskal(data.getMatrixCost(), dimension, upperBound, 1, node.multipliersSum, node.multipliers, node.forbiddenEdges);
    double objValue = kruskal.getCost();
    std::cout << "LB: " << objValue << " UB: " << upperBound << "\n";
    if (objValue > upperBound) {
        node.prune = true;
    } else {
        node.prune = false;
        node.lowerBound = objValue;
        node.multipliers = kruskal.getMultipliers();
        node.multipliersSum = kruskal.getMultipliersSum();
    }
    return kruskal.getNewForbiddenEdges();
}

int chooseSubtour(const std::vector<std::vector<int>> &subtours) {
    int smallestSubtour = std::numeric_limits<int>::max(), initialIndex = std::numeric_limits<int>::max(), chosen = 0;
    std::vector<int> subtourSize = computeSubtourSize(subtours);
    // Escolher o índice do subtour a ser avaliado
    for (int i = 0; i < subtours.size(); i++) {
        if (subtourSize[i] < smallestSubtour) {
            smallestSubtour = subtourSize[i];
            initialIndex = subtours[i][0];
            chosen = i;
        } else if (subtourSize[i] == smallestSubtour and subtours[i][0] < initialIndex) {
            initialIndex = subtours[i][0];
            chosen = i;
        }
    }
    
    return chosen;
}

std::vector<int> computeSubtourSize(const std::vector<std::vector<int>> &subtours) {
    std::vector<int> size(subtours.size(), 0);
    for (int i = 0; i < subtours.size(); i++) {
        size[i] = (int) subtours[i].size();
    }
    return size;
}

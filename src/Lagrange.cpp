//
//  Kruskal.cpp
//  Opt
//
//  Created by David Galvao on 13/09/22.
//
#include "Lagrange.hpp"
#include "Kruskal.hpp"

#define EPSILON 0.000001

Lagrange::Lagrange(int dimension) : Lagrange(std::vector<std::pair<int, int>>(), std::vector<double>(dimension, 0)) { }

Lagrange::Lagrange(std::vector<std::pair<int, int>> forbiddenEdges, std::vector<double> multipliers) : forbiddenEdges(forbiddenEdges), multipliers(multipliers), multipliersSum(0), dimension((int) multipliers.size()), prune(false) {
    subgradient = std::vector<double>(dimension);
    for(int i = 0; i < multipliers.size(); i++) {
        multipliersSum += multipliers[i];
    }
}

void Lagrange::addForbiddenEdge(std::pair<int, int> edge) {
    forbiddenEdges.push_back(edge);
}

void Lagrange::solve(double **costMatrix, const double w, const double e, const int maxItr) {
    // Alocar a matriz
    double **lagrangeMatrix = new double*[dimension];
    for(int i = 0; i < dimension; i++) lagrangeMatrix[i] = new double[dimension];
    // Variáveis para guardar os dados da melhor iteração
    double cost = 0, bestCost = INFINITE * -1;
    std::vector<double> bestMultipliers = multipliers;
    Vertex chosenVertex = NULL;
    // Loop principal
    for(int itr = 0; abs(bestCost - w) > EPSILON; itr++) {
        if (itr >= maxItr) break;
        for(int i = 0; i < dimension; i++) {
            for(int j = i + 1; j < dimension; j++) {
                lagrangeMatrix[i][j] = costMatrix[i][j] - multipliers[i] - multipliers[j];
            }
        }
        for(auto edge : forbiddenEdges) {
            if(edge.first < edge.second) {
                lagrangeMatrix[edge.first][edge.second] = INFINITE;
            } else {
                lagrangeMatrix[edge.second][edge.first] = INFINITE;
            }
        }
        Kruskal kruskal = Kruskal(lagrangeMatrix, dimension);
        cost = kruskal.getCost();
        if(cost == INFINITE) {
            prune = true;
            break;
        }
//        std::cout << "Itr: " << itr << " Cost: " << cost << "\n";
        bool best = false;
        if(cost > bestCost) {
            itr = 0;
            best = true;
            bestCost = cost;
            bestMultipliers = multipliers;
            branchEdges.clear();
        }
        // Atualizar subgradiente e o custo
        subgradientSum = 0;
        feasible = true;
        int highestDegree = -1;
        for(int i = 0; i < dimension; i++) {
            Vertex vertex = kruskal.getVertex(i);
            subgradient[i] = 2 - vertex.degree;
            subgradientSum += subgradient[i] * subgradient[i];
            cost += 2 * multipliers[i];
            if(vertex.degree > highestDegree && best) {
                highestDegree = vertex.degree;
                chosenVertex = vertex;
            }
            if(subgradient[i] != 0 && feasible) feasible = false;
        }
//        std::cout << subgradientSum << "\n";
        // Checar se a solução é viável
        if(feasible) {
//            std::cout << "Feasible\n";
            prune = true;
            break;
        }
        // Atualizar multiplicadores
        multipliersSum = 0;
        double constantTerm = e * ((w - cost) / subgradientSum);
        std::vector<double> newMultipliers;
        for(int i = 0; i < dimension; i++) {
            newMultipliers.push_back(multipliers[i] + constantTerm * subgradient[i]);
            multipliersSum += multipliers[i] + constantTerm * subgradient[i];
        }
        multipliers = newMultipliers;
    }
    // Atualizar as variáveis de acordo com a melhor solução encontrada
    multipliers = bestMultipliers;
    lowerBound = bestCost;
    for(auto adj : chosenVertex.adjacencyList) branchEdges.push_back(std::make_pair(chosenVertex.index, adj));
    if(abs(bestCost - w) <= EPSILON) {
        prune = true;
    }
    // Desalocar a matriz
    for(int i = 0; i < dimension; i++) delete [] lagrangeMatrix[i];
    delete [] lagrangeMatrix;
}

/*
#include "Lagrange.hpp"

#include "Kruskal.hpp"
#include "Util.hpp"

Kruskal::Kruskal(Data &data, double w, double e) : Kruskal(data.getMatrixCost(), data.getDimension(), w, e, 0, std::vector<double>(data.getDimension(), 0), std::vector<std::pair<int, int>>()) { }

Kruskal::Kruskal(double **distMatrix, int dimension, double w, double e, double multipliersSum, std::vector<double> multipliers, const std::vector<std::pair<int, int>> &forbiddenEdges) : costMatrix(new double*[dimension]), cost(0), subgradientSum(0), dimension(dimension), w(w), e(e), multipliersSum(multipliersSum), multipliers(multipliers), validSolution(false) {
    for(int i = 0; i < dimension; i++) {
        costMatrix[i] = new double[dimension];
        for(int j = i + 1; j < dimension; j++) {
            costMatrix[i][j] = distMatrix[i][j];
        }
    }
    for(std::pair<int, int> edge : forbiddenEdges) {
        costMatrix[edge.first][edge.second] = INFINITE;
        costMatrix[edge.second][edge.first] = INFINITE;
    }
    for(int itr = 0; itr < 1; itr++) {
        vertices.clear();
        sets.clear();
        mst.clear();
        for(int i = 0; i < dimension; i++) {
            vertices.push_back(Vertex(i, 0));
            sets.push_back(Set(i, 0));
            if(itr == 0) {
                subgradient.push_back(0);
            }
        }
        build1Tree();
        connectFirstVertex();
        updateSubgradient();
        updateMultipliers();
        std::cout << "(Iteração " << itr << ") Custo: " << getCost() << "\n";
        if(itr > 0 and getCost() < w) { break; }
//        updateCostMatrix();
//        printTreeInfo();
    }
}

Kruskal::Kruskal(double **distMatrix, int dimension, const std::vector<std::pair<int, int>> &forbiddenEdges) : costMatrix(new double*[dimension]), cost(0), subgradientSum(0), dimension(dimension), multipliers(std::vector<double>(dimension, 0)), multipliersSum(0), validSolution(false) {
    for(int i = 0; i < dimension; i++) {
        costMatrix[i] = new double[dimension];
        vertices.push_back(Vertex(i, 0));
        sets.push_back(Set(i, 0));
        for(int j = 0; j < dimension; j++) {
            costMatrix[i][j] = distMatrix[i][j];
        }
    }
    std::cout << "Forbidden: ";
    for(std::pair<int, int> edge : forbiddenEdges) {
        std::cout << "(" << edge.first << "," << edge.second << ") ";
        costMatrix[edge.first][edge.second] = INFINITE;
//        costMatrix[edge.second][edge.first] = INFINITE;
    }
    std::cout << "\n";
    build1Tree();
    if(cost != INFINITE) connectFirstVertex();
    //printTreeInfo();
}

Kruskal::~Kruskal() {
    for(int i = 0; i < dimension; i++) delete [] costMatrix[i];
    delete [] costMatrix;
}

int Kruskal::findSet(int index) {
    if(sets[index].parent != index) {
        sets[index].parent = findSet(sets[index].parent);
    }
    return sets[index].parent;
}

void Kruskal::unionSet(int a, int b) {
    int rootA = findSet(a), rootB = findSet(b);
    if(sets[rootA].value < sets[rootB].value) {
        sets[rootA].parent = rootB;
    } else if(sets[rootA].value > sets[rootB].value) {
        sets[rootB].parent = rootA;
    } else {
        sets[rootB].parent = rootA;
        sets[rootA].value++;
    }
}

std::vector<Edge> Kruskal::buildEdges() {
    std::vector<Edge> edges;
    for(int i = 1; i < dimension; i++) {
        for(int j = i + 1; j < dimension; j++) {
            if(costMatrix[i][j] < INFINITE) {
                edges.push_back(Edge(i, j, getEdgeWeight(i, j)));
            }
        }
    }
    return edges;
}

void Kruskal::printTreeInfo() {
    for(Vertex vertex : vertices) {
        std::cout << vertex.index << ":";
        for(int adj : vertex.adjacencyList) {
            std::cout << " " << adj;
        }
        std::cout << " (Degree = " << vertex.degree << ")\n";
    }
    std::cout << "Cost: " << getCost() << "\n";
//    std::cout << "U = (";
//    for(int i = 0; i < dimension - 1; i++) {
//        std::cout << multipliers[i] << ", ";
//    }
//    std::cout << multipliers[dimension - 1] << ")\n";
//    std::cout << "G = (";
//    for(int i = 0; i < dimension - 1; i++) {
//        std::cout << subgradient[i] << ", ";
//    }
//    std::cout << subgradient[dimension - 1] << ")\n";
}

// Atualizar os subgradientes
void Kruskal::updateSubgradient() {
    subgradientSum = 0;
    for(int i = 0; i < dimension; i++) {
        subgradient[i] = 2 - vertices[i].degree;
        subgradientSum += subgradient[i] * subgradient[i];
    }
}

void Kruskal::updateMultipliers() {
    std::vector<double> newMultipliers;
    multipliersSum = 0;
    double constantTerm = e * ((w - cost) / subgradientSum);
    for(int i = 0; i < dimension; i++) {
        newMultipliers.push_back(multipliers[i] + constantTerm * subgradient[i]);
        multipliersSum += multipliers[i] + constantTerm * subgradient[i];
    }
    multipliers = newMultipliers;
}

void Kruskal::build1Tree() {
    std::vector<Edge> edges = buildEdges();
    // Construir a 1-Arvore
    std::sort(edges.begin(), edges.end());
    for(Edge edge : edges) {
        if (findSet(edge.a) != findSet(edge.b)) {
            mst.push_back(edge);
            cost += edge.weight;
            unionSet(edge.a, edge.b);
            vertices[edge.a].degree++;
            vertices[edge.a].adjacencyList.push_back(edge.b);
            vertices[edge.b].degree++;
            vertices[edge.b].adjacencyList.push_back(edge.a);
            if(mst.size() == dimension - 2) return;
        }
    }
    cost = INFINITE;
    //isTree(dimension, mst);
}

// Conectar o primeiro vértice a arvore
void Kruskal::connectFirstVertex() {
    std::vector<Edge> edges1;
    for(int i = 1; i < dimension; i++) {
        if(costMatrix[0][i] < INFINITE) {
            edges1.push_back(Edge(0, i, getEdgeWeight(0, i)));
        }
    }
    if(edges1.size() < 2) {
        cost = INFINITE;
        return;
    }
    // Adicionar as duas arestas menos custosas que saem do vértice 1
    std::sort(edges1.begin(), edges1.end());
    // Adicionar a aresta mais próxima
    mst.push_back(edges1[0]);
    cost += edges1[0].weight;
    vertices[edges1[0].a].degree++;
    vertices[edges1[0].a].adjacencyList.push_back(edges1[0].b);
    vertices[edges1[0].b].degree++;
    vertices[edges1[0].b].adjacencyList.push_back(edges1[0].a);
    // Adicionar a segunda aresta mais próxima
    mst.push_back(edges1[1]);
    cost += edges1[1].weight;
    vertices[edges1[1].a].degree++;
    vertices[edges1[1].a].adjacencyList.push_back(edges1[1].b);
    vertices[edges1[1].b].degree++;
    vertices[edges1[1].b].adjacencyList.push_back(edges1[1].a);
}

//void Kruskal::updateCostMatrix() {
//    for(int i = 0; i < dimension; i++) {
//        for(int j = 0; j < dimension; j++) {
//            if(i != j and costMatrix[i][j] < INFINITE) {
//                costMatrix[i][j] = costMatrix[i][j] - multipliers[i] - multipliers[j];
//            }
//        }
//    }
//}

double Kruskal::getCost() {
    return cost + 2 * multipliersSum;
}

double Kruskal::getMultipliersSum() {
    return multipliersSum;
}

std::vector<double> Kruskal::getMultipliers() {
    return multipliers;
}

std::vector<std::pair<int, int>> Kruskal::getNewForbiddenEdges() {
    Vertex *highestDegree = &vertices[0];
    double lowestDegree = dimension + 1;
    std::vector<std::pair<int, int>> forbiddenEdges;
    for(int i = 1; i < dimension; i++) {
        if(vertices[i].degree > highestDegree->degree) {
            highestDegree = &vertices[i];
        } else if(vertices[i].degree < lowestDegree) {
            lowestDegree = vertices[i].degree;
        }
    }
//    std::cout << lowestDegree << " " << highestDegree->degree << "\n";
    for(int adj : highestDegree->adjacencyList) {
        std::pair<int, int> forbiddenEdge;
        if(highestDegree->index > adj) {
            forbiddenEdge.first = adj;
            forbiddenEdge.second = highestDegree->index;
        } else {
            forbiddenEdge.first = highestDegree->index;
            forbiddenEdge.second = adj;
        }
        forbiddenEdges.push_back(forbiddenEdge);
    }
    return forbiddenEdges;
}

double Kruskal::getEdgeWeight(int a, int b) {
    if(costMatrix[a][b] >= INFINITE or a == b) return INFINITE;
    return costMatrix[a][b] - multipliers[a] - multipliers[b];
//    return costMatrix[a][b];
}

bool Kruskal::isValidSolution() {
    if(vertices[0].degree != 2) return false;
    std::vector<bool> visited(dimension, false);
    std::stack<int> visitStack;
    visitStack.push(0);
    while(!visitStack.empty()) {
        int current = visitStack.top();
        if(vertices[current].degree != 2) return false;
        visitStack.pop();
        visited[current] = true;
        for(int adj: vertices[current].adjacencyList) {
            if(!visited[adj]) {
                visitStack.push(adj);
            }
        }
    }
    for(bool visit : visited) {
        if(!visit) return false;
    }
    printTreeInfo();
    return true;
}
*/

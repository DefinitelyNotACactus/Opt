//
//  Kruskal.cpp
//  Opt
//
//  Created by David Galvao on 13/09/22.
//

#include "Kruskal.hpp"

Kruskal::Kruskal(Data &data, double w, double e) : Kruskal(data.getMatrixCost(), data.getDimension(), w, e, 0, std::vector<double>(data.getDimension(), 0), std::vector<std::pair<int, int>>()) { }

Kruskal::Kruskal(double **distMatrix, int dimension, double w, double e, double multipliersSum, std::vector<double> multipliers, const std::vector<std::pair<int, int>> &forbiddenEdges) : costMatrix(new double*[dimension]), cost(0), subgradientSum(0), dimension(dimension), w(w), e(e), multipliersSum(multipliersSum), multipliers(multipliers) {
//    this->costMatrix = new double*[dimension];
    for(int i = 0; i < dimension; i++) {
        costMatrix[i] = new double[dimension];
        for(int j = 0; j < dimension; j++) {
            costMatrix[i][j] = distMatrix[i][j] - multipliers[i] - multipliers[j];
        }
    }
    for(std::pair<int, int> edge : forbiddenEdges) {
        costMatrix[edge.first][edge.second] = INFINITE;
    }
    for(int itr = 0; getCost() < w; itr++) {
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
        updateCostMatrix();
        std::cout << "(Iteração " << itr << ")\n";
        printTreeInfo();
    }
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
                edges.push_back(Edge(i, j, costMatrix[i][j]));
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
    std::cout << "U = (";
    for(int i = 0; i < dimension - 1; i++) {
        std::cout << multipliers[i] << ", ";
    }
    std::cout << multipliers[dimension - 1] << ")\n";
    std::cout << "G = (";
    for(int i = 0; i < dimension - 1; i++) {
        std::cout << subgradient[i] << ", ";
    }
    std::cout << subgradient[dimension - 1] << ")\n";
}

// Atualizar os subgradientes
void Kruskal::updateSubgradient() {
    subgradientSum = 0;
    for(int i = 0; i < dimension; i++) {
        subgradient[i] = 2 - vertices[i].degree;
        subgradientSum += (2 - vertices[i].degree) * (2 - vertices[i].degree);
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
    for(int i = 0; mst.size() < dimension - 2; i++) {
        Edge edge = edges[i];
        int a = findSet(edge.a), b = findSet(edge.b);
        if (a != b) {
            mst.push_back(edge);
            cost += edge.weight;
            unionSet(a, b);
            vertices[a].degree++;
            vertices[a].adjacencyList.push_back(b);
            vertices[b].degree++;
            vertices[b].adjacencyList.push_back(a);
        }
    }
}

// Conectar o primeiro vértice a arvore
void Kruskal::connectFirstVertex() {
    std::vector<Edge> edges1;
    for(int i = 1; i < dimension; i++) {
        if(costMatrix[0][i] < INFINITE) {
            edges1.push_back(Edge(0, i, costMatrix[0][i]));
        }
    }
    // Adicionar as duas arestas menos custosas que saem do vértice 1
    std::sort(edges1.begin(), edges1.end());
    mst.push_back(edges1[0]);
    cost += edges1[0].weight;
    vertices[edges1[0].a].degree++;
    vertices[edges1[0].a].adjacencyList.push_back(edges1[0].b);
    vertices[edges1[0].b].degree++;
    vertices[edges1[0].b].adjacencyList.push_back(edges1[0].a);
    mst.push_back(edges1[1]);
    cost += edges1[1].weight;
    vertices[edges1[1].a].degree++;
    vertices[edges1[1].a].adjacencyList.push_back(edges1[1].b);
    vertices[edges1[1].b].degree++;
    vertices[edges1[1].b].adjacencyList.push_back(edges1[1].a);
}

void Kruskal::updateCostMatrix() {
    for(int i = 0; i < dimension; i++) {
        for(int j = 0; j < dimension; j++) {
            if(i != j and costMatrix[i][j] < INFINITE) {
                costMatrix[i][j] = costMatrix[i][j] - multipliers[i] - multipliers[j];
            }
        }
    }
}

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
    std::vector<std::pair<int, int>> forbiddenEdges;
    for(int i = 1; i < dimension; i++) {
        if(vertices[i].degree > highestDegree->degree) {
            highestDegree = &vertices[i];
        }
    }
    for(int adj : highestDegree->adjacencyList) {
        std::pair<int, int> forbiddenEdge;
        forbiddenEdge.first = highestDegree->index;
        forbiddenEdge.second = adj;
        forbiddenEdges.push_back(forbiddenEdge);
    }
    return forbiddenEdges;
}

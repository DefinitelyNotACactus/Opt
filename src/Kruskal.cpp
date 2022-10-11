//
//  Kruskal.cpp
//  Opt
//
//  Created by David Galvao on 13/09/22.
//

#include "Kruskal.hpp"
#include "Util.hpp"

Kruskal::Kruskal(double **distMatrix, int dimension) : dimension(dimension), cost(0) {
    for(int i = 0; i < dimension; i++) {
        vertices.push_back(Vertex(i));
        sets.push_back(i);
    }
    // Obter as arestas do grafo
    std::vector<Edge> edges;
    for(int i = 1; i < dimension; i++) {
        for(int j = i + 1; j < dimension; j++) {
            if(distMatrix[i][j] < INFINITE) {
                edges.push_back(Edge(i, j, distMatrix[i][j]));
            }
        }
    }
    // Construir a 1-Arvore
    std::sort(edges.begin(), edges.end(), std::less<Edge>());
    for(Edge edge : edges) {
        if (findSet(edge.a) != findSet(edge.b)) {
            unionSet(edge.a, edge.b);
            insertEdge(edge);
            if(mst.size() == dimension - 2) break;
        }
    }
    if(mst.size() != dimension - 2) {
        cost = INFINITE;
        return;
    }
    std::vector<Edge> edges1;
    for(int i = 1; i < dimension; i++) {
        if(distMatrix[0][i] < INFINITE) {
            edges1.push_back(Edge(0, i, distMatrix[0][i]));
        }
    }
    if(edges1.size() < 2) {
        cost = INFINITE;
        return;
    }
    // Adicionar as duas arestas menos custosas que saem do vértice 1
    std::sort(edges1.begin(), edges1.end(), std::less<Edge>());
    insertEdge(edges1[0]);
    insertEdge(edges1[1]);
}

void Kruskal::insertEdge(Edge edge) {
    mst.push_back(edge);
    cost += edge.weight;
    vertices[edge.a].degree++;
    vertices[edge.a].adjacencyList.push_back(edge.b);
    vertices[edge.b].degree++;
    vertices[edge.b].adjacencyList.push_back(edge.a);
}

int Kruskal::findSet(int index) {
    if(sets[index] != index) {
        sets[index] = findSet(sets[index]);
    }
    return sets[index];
}

void Kruskal::unionSet(int a, int b) {
    sets[findSet(a)] = findSet(b);
}

double Kruskal::getCost() {
    return cost;
}

Vertex Kruskal::getVertex(int index) {
    return vertices[index];
}

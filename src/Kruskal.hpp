//
//  Kruskal.hpp
//  Opt
//
//  Created by David Galvao on 13/09/22.
//

#ifndef Kruskal_hpp
#define Kruskal_hpp

#include "Heuristics.hpp"
#include "data.hpp"

struct Vertex {
    int index;
    int degree;
    std::vector<int> adjacencyList;
    
    Vertex(int index): index(index), degree(0) { };
};

struct Edge {
    int a, b;
    double weight;
    
    Edge(int a, int b, double weight): a(a), b(b), weight(weight) { };
    bool operator <(const Edge &b) const { return weight < b.weight; }
    bool operator >(const Edge &b) const { return weight > b.weight; }
};

class Kruskal {
public:
    Kruskal(double **, int);
    
    void printTreeInfo();
    double getCost();
    Vertex getVertex(int);
    
private:
    std::vector<Vertex> vertices;
    std::vector<Edge> mst;
    std::vector<int> sets;
    double cost;
    const int dimension;
    
    int findSet(int);
    void unionSet(int, int);
    void insertEdge(Edge);
};

#endif /* Kruskal_hpp */

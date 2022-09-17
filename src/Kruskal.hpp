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
    
    Vertex(int index, int degree): index(index), degree(degree) { };
};

struct Set {
    int parent;
    int value;
    
    Set(int parent, int value): parent(parent), value(value) { };
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
    Kruskal(Data &, double, double);
    Kruskal(double **, int, double, double, double, std::vector<double>, const std::vector<std::pair<int, int>> &);
    ~Kruskal();
    
    void printTreeInfo();
    double getCost();
    double getMultipliersSum();
    double **getCostMatrix() { return costMatrix; };
    std::vector<double> getMultipliers();
    std::vector<std::pair<int, int>> getNewForbiddenEdges();
    
private:
    std::vector<Vertex> vertices;
    std::vector<Set> sets;
    std::vector<Edge> mst;
    std::vector<double> multipliers, subgradient;
    double e, w, cost, subgradientSum, multipliersSum;
    double **costMatrix;
    int dimension;
    
    int findSet(int);
    void unionSet(int, int);
    void build1Tree();
    void connectFirstVertex();
    void updateSubgradient();
    void updateMultipliers();
    void updateCostMatrix();
    std::vector<Edge> buildEdges();
};

#endif /* Kruskal_hpp */

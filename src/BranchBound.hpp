//
//  BranchBound.hpp
//  Opt
//
//  Created by David Galvao on 28/08/22.
//

#ifndef BranchBound_hpp
#define BranchBound_hpp

#include <vector>
#include <chrono>

#include "hungarian.hpp"
#include "Heuristics.hpp"
#include "data.hpp"
#include "Lagrange.hpp"

struct Node {
    std::vector<std::pair<int, int>> forbiddenEdges;
    std::vector<std::vector<int>> subtours;
    double lowerBound;
    int chosen;
    bool prune;
    
    bool operator <(const Node &b) const { return lowerBound < b.lowerBound; }
    bool operator >(const Node &b) const { return lowerBound > b.lowerBound; }
};

struct Tree {
    std::vector<Node> nodes;
    Solution solution;
    double upperBound;
    std::chrono::high_resolution_clock::time_point start;
};

std::vector<std::vector<int>> getSubtours(hungarian_problem_t &);
void printSubtours(const std::vector<std::vector<int>> &);
void printSubtour(const std::vector<int> &);
Tree buildTree(Data &, double = -1, int = 0);
void traverseTreeDFS(Data &, Tree &, Node &);
void traverseTreeDFSLagrange(Data &, Tree &, Lagrange &);
void traverseTreeBestBound(Data &, Tree &, Node &);
void traverseTreeBestBoundLagrange(Data &, Tree &, Lagrange &);
void traverseTreeBFS(Data &, Tree &, Node &);
void computeSolution(Node &, Data &, int, double);
std::vector<std::pair<int, int>> computeSolutionLagrangean(Node &, Data &, int, double);
int chooseSubtour(const std::vector<std::vector<int>> &);
std::vector<int> computeSubtourSize(const std::vector<std::vector<int>> &);
#endif /* BranchBound_hpp */

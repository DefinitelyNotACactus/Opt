//
//  Heuristics.hpp
//  Opt
//
//  Created by David Galvao on 28/08/22.
//

#ifndef Heuristics_hpp
#define Heuristics_hpp

#include "data.hpp"

#include <vector>

struct Solution {
    std::vector<int> tour;
    double cost;
};

void nearestNeighbor(Solution &solution, double **matrix, int dimension);
bool swap(Data &data, std::vector<int> &st, double &cost);
bool twoOpt(Data &data, std::vector<int> &st, double &cost);
bool reinsertion(Data &data, std::vector<int> &st, int subtourSize, double &cost);
void rvnd(Data &data, std::vector<int> &s, double &cost);
#endif /* Heuristics_hpp */

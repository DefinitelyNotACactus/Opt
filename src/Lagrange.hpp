//
//  Lagrange.hpp
//  Opt
//
//  Created by David Galvao on 10/10/22.
//

#ifndef Lagrange_hpp
#define Lagrange_hpp

#include <vector>

class Lagrange {
public:
    Lagrange(int);
    Lagrange(std::vector<std::pair<int, int>>, std::vector<double>);
    
    void addForbiddenEdge(std::pair<int, int>);
    void solve(double **, const double, const double = 1, const int = 10);
    
    inline double getLB() { return lowerBound; }
    inline double getMultipliersSum() { return multipliersSum; }
    
    inline std::vector<int> getTour() { return tour; }
    inline std::vector<double> getMultipliers() { return multipliers; }
    
    inline bool isFeasible() { return feasible; }
    inline bool shouldPrune() { return prune; }
    
    inline std::vector<std::pair<int, int>> getBranches() { return branchEdges; }
    inline std::vector<std::pair<int, int>> getForbiddenEdges() { return forbiddenEdges; }
    
    inline bool operator <(const Lagrange &b) const { return lowerBound < b.lowerBound; }
    inline bool operator >(const Lagrange &b) const { return lowerBound > b.lowerBound; }
    
private:
    std::vector<int> tour;
    std::vector<std::pair<int, int>> forbiddenEdges, branchEdges;
    std::vector<double> multipliers, subgradient;
    double lowerBound, multipliersSum, subgradientSum;
    int dimension;
    bool prune, feasible;
    
};
#endif /* Lagrange_hpp */

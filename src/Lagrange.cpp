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

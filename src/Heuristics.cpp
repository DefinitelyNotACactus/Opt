//
//  Heuristics.cpp
//  Opt
//
//  Created by David Galvao on 28/08/22.
//

#include "Heuristics.hpp"
#include "data.hpp"

#include <iostream>

// Heurística construtiva, vizinho mais próximo
void nearestNeighbor(Solution &solution, double **matrix, int dimension) {
    std::vector<bool> visited(dimension, false);
    double nearestDistance, totalCost = 0;
    int nearestNeighbor, currentPoint;
    // Iniciamos a solução no vértice 1
    solution.tour.clear();
    solution.tour.push_back(0);
    visited[0] = true;
    currentPoint = 0;
//    std::cout << "Construindo solução (vizinho mais próximo)\n0->";
    for(int i = 0; i < dimension - 1; i++) {
        nearestDistance = std::numeric_limits<double>::infinity();
        nearestNeighbor = currentPoint;
        for(int j = 0; j < dimension; j++) {
            if (visited[j] == false and j != currentPoint and matrix[currentPoint][j] < nearestDistance) {
                nearestDistance = matrix[currentPoint][j];
                nearestNeighbor = j;
            }
        }
        if (nearestNeighbor != currentPoint) {
//            std::cout << nearestNeighbor << "->";
            solution.tour.push_back(nearestNeighbor);
            totalCost += nearestDistance;
            visited[nearestNeighbor] = true;
            currentPoint = nearestNeighbor;
        }
    }
//    for (int i = 0; i < dimension; i++) {
//        if (!visited[i]) { std::cout << i << " não visitado!\n"; }
//    }
    // Fechar o tour
//    std::cout << "0\n";
    totalCost += matrix[currentPoint][0];
    solution.tour.push_back(0);
    solution.cost = totalCost;
//    std::cout << "Custo total: " << totalCost << "\n";
}

// Movimento de vizinhança: troca
bool swap(Data &data, std::vector<int> &st, double &cost) {
    double delta = std::numeric_limits<double>::infinity(), nDelta = 0;
    int a = 0, b = 0, dimension = data.getDimension();
    for(int i = 1; i < dimension - 1; i++) {
        for(int j = i + 1; j < dimension - 1; j++) {
            if(i == j - 1) {
                nDelta = data.getDistance(st.at(i - 1), st.at(j)) + data.getDistance(st.at(j + 1), st.at(i)) - data.getDistance(st.at(i - 1), st.at(i)) - data.getDistance(st.at(j + 1), st.at(j));
            } else {
                nDelta = data.getDistance(st.at(i), st.at(j - 1)) + data.getDistance(st.at(i), st.at(j + 1)) + data.getDistance(st.at(j), st.at(i - 1)) + data.getDistance(st.at(j), st.at(i + 1)) - data.getDistance(st.at(j), st.at(j - 1)) - data.getDistance(st.at(j), st.at(j + 1)) - data.getDistance(st.at(i), st.at(i - 1)) - data.getDistance(st.at(i), st.at(i + 1));
            }
            if(nDelta < delta) {
                a = i;
                b = j;
                delta = nDelta;
            }
        }
    }
    
    if(delta < 0) {
        int aux = st.at(a);
        
        st.at(a) = st.at(b);
        st.at(b) = aux;
        
        cost += delta;
        return true;
    }
    return false;
}

bool twoOpt(Data &data, std::vector<int> &st, double &cost) {
    double delta = std::numeric_limits<double>::infinity(), nDelta = 0;
    int a = 0, b = 0;
    for(int i = 1; i < st.size() - 1; i++) {
        for(int j = i + 1; j < st.size() - 1; j++) {
            nDelta = data.getDistance(st.at(i - 1), st.at(j)) + data.getDistance(st.at(j + 1), st.at(i)) - data.getDistance(st.at(j + 1), st.at(j)) - data.getDistance(st.at(i - 1), st.at(i));
            if(nDelta < delta) {
                a = i;
                b = j;
                delta = nDelta;
            }
        }
    }
    
    if(delta < 0) {
        vector<int> aux; // assumindo que a > b
        for(int i = b; i >= a; i--) {
            aux.push_back(st.at(i));
        }
        st.erase(st.begin() + a, st.begin() + b + 1);
        st.insert(st.begin() + a, aux.begin(), aux.end());
        
        cost += delta;
        return true;
    }
    return false;
}

bool reinsertion(Data &data, std::vector<int> &st, int subtourSize, double &cost) {
    double delta = std::numeric_limits<double>::infinity(), nDelta = 0;
    int a = 0, b = 0;
    for(int i = 1; i < st.size() - subtourSize; i++) {
        for(int j = 1; j < st.size() - 1; j++) {
            if(j >= i && j <= i + subtourSize) {
                continue;
            }
            nDelta = data.getDistance(st.at(j - 1), st.at(i)) + data.getDistance(st.at(i + subtourSize - 1), st.at(j)) + data.getDistance(st.at(i - 1), st.at(i + subtourSize)) - data.getDistance(st.at(i + subtourSize - 1), st.at(i + subtourSize)) - data.getDistance(st.at(j - 1), st.at(j)) - data.getDistance(st.at(i - 1), st.at(i));
            if(nDelta < delta) {
                a = i;
                b = j;
                delta = nDelta;
            }
        }
    }
    
    if(delta < 0) {
        std::vector<int> aux;
        aux.insert(aux.end(), st.begin() + a, st.begin() + a + subtourSize);
        st.erase(st.begin() + a, st.begin() + a + subtourSize);
        if(a >= b) {
            st.insert(st.begin() + b, aux.begin(), aux.end());
        } else {
            st.insert(st.begin() + b - subtourSize, aux.begin(), aux.end());
        }
        
        cost += delta;
        return true;
    }
    return false;
}

void rvnd(Data &data, std::vector<int> &s, double &cost) {
    std::vector<int> nl = {0, 1, 2, 3, 4};
    double tempCost = cost;
    
    int n = 0;
    bool improve = false;
    while(!nl.empty()) {
        n = rand() % nl.size();
        switch(nl.at(n)) {
            case 0: // swap
                improve = swap(data, s, tempCost);
                break;
            case 1: // 2-opt
                improve = twoOpt(data, s, tempCost);
                break;
            case 2: // reinsercao 1
                improve = reinsertion(data, s, 1, tempCost);
                break;
            case 3: // reinsercao 2
                improve = reinsertion(data, s, 2, tempCost);
                break;
            case 4: // reinsercao 3
                improve = reinsertion(data, s, 3, tempCost);
                break;
        }
        if(improve) {
            nl = {0, 1, 2, 3, 4};
        } else {
            nl.erase(nl.begin() + n);
        }
    }

    cost = tempCost;
}

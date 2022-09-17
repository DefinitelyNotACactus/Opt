#include <iostream>
#include <vector>
#include <chrono>

using namespace std;

#include "data.hpp"
#include "hungarian.hpp"
#include "BranchBound.hpp"
#include "Kruskal.hpp"

int main(int argc, char** argv) {
    int strategy = 0;
    double upperBound = -1;
    if (argc > 2) {
        upperBound = atof(argv[2]);
    }
    if (argc > 3) {
        strategy = atoi(argv[3]);
    }
    
    Data data(argc, argv[1]);
    data.readData();
    
    std::vector<std::vector<double>> costVector;
    double **cost = new double*[data.getDimension()];
    for (int i = 0; i < data.getDimension(); i++){
        std::vector<double> costVectorLine;
        cost[i] = new double[data.getDimension()];
        for (int j = 0; j < data.getDimension(); j++){
            cost[i][j] = data.getDistance(i,j);
            costVectorLine.push_back(data.getDistance(i, j));
        }
        costVector.push_back(costVectorLine);
    }

//    hungarian_problem_t p;
//    int mode = HUNGARIAN_MODE_MINIMIZE_COST;
//    hungarian_init(&p, cost, data.getDimension(), data.getDimension(), mode); // Carregando o problema
//
//    double obj_value = hungarian_solve(&p);
//    cout << "Obj. value: " << obj_value << endl;
//
//    cout << "Assignment" << endl;
//    hungarian_print_status(&p);
//    std::vector<std::vector<int>> subtours = getSubtours(p);
//    printSubtours(subtours);
    auto start = std::chrono::steady_clock::now();
    Tree tree = buildTree(data, upperBound, strategy);
    auto end = std::chrono::steady_clock::now();
//    std::cout << chrono::duration_cast<chrono::milliseconds>(end - start).count() << ";" << tree.upperBound << "\n";
//    hungarian_free(&p);
    
//    Kruskal *kruskal = new Kruskal(data, 7000, 1);
//    std::vector<std::pair<int, int>> forbidden = kruskal->getNewForbiddenEdges();
//    for(auto pair : forbidden) {
//        std::cout << "(" << pair.first << ", " << pair.second << ")\n";
//    }
    for (int i = 0; i < data.getDimension(); i++) delete [] cost[i];
    delete [] cost;
//    delete kruskal;
    return 0;
}

#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <vector>
#include "Instance.hpp"

using namespace std;

struct Route {
    vector<int>stations;
    int capacity;
};

class Solver {
    private:
        Route route;        
    public:
        Solver() : route{} {};


        void Solve(Instance *instance);
};

vector<Route> clarkeWright(const Instance *instance);
bool validateRoutes(const vector<Route> &routes, const Instance *instance);
void printRoute(const Route &route);
void printAllRoutes(const vector<Route> &routes, const Instance *instance);

#endif
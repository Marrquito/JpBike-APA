#include "Solver.hpp"
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <chrono>

#define QTD_VIZINHANCAS 3

typedef enum {
    SWAP = 1,
    TWO_OPT,
    REINSERTION
} Vizinhança;

vector<Route> VND(vector<Route> routes, const Instance *instance);
vector<Route> performSwap(vector<Route> routes, const Instance *instance);
vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance);
vector<Route> performReinsertion(vector<Route> routes, const Instance *instance);

vector<Route> GRASPConstruction(const Instance *instance);

double calculateTotalCost(const vector<Route> &routes, const Instance *instance) {
    if (routes.empty()) return 0.0;
    
    const auto &distancias  = instance->get_distancias();
    double totalCost        = 0.0;
    
    for (const auto &route : routes) {
        const auto &stations = route.stations;
        const int stationsSize = (int)stations.size();
        
        // [0, station, 0]
        if      (stationsSize == 3) totalCost += distancias[stations[0]][stations[1]] + distancias[stations[1]][stations[2]];
        else if (stationsSize == 4) { // [0, s1, s2, 0]
            totalCost += distancias[stations[0]][stations[1]] 
                       + distancias[stations[1]][stations[2]]
                       + distancias[stations[2]][stations[3]];
        } else {
            // rotas maiores
            for (int i = 0; i < stationsSize - 1; i++) {
                totalCost += distancias[stations[i]][stations[i + 1]];
            }
        }
    }

    return totalCost;
}

int calculateRouteCapacity(const vector<int> &stations, const vector<int> &demandas) {
    if (stations.size() <= 2) return 0; // rota só com depósito

    int currentLoad = 0;
    int minPrefix   = 0;
    int maxPrefix   = 0;

    // percorre apenas as estações (ignora depósito no início e fim)
    for (int i = 1; i < (int)stations.size() - 1; i++) {
        int node    = stations[i];
        int demand  = demandas[node - 1]; // node começa em 1

        currentLoad += demand;

        if (currentLoad < minPrefix) minPrefix = currentLoad;
        if (currentLoad > maxPrefix) maxPrefix = currentLoad;
    }

    int initialLoad = (minPrefix < 0) ? -minPrefix : 0;
    int maxLoad     = initialLoad + maxPrefix;

    return maxLoad;
}

bool validateRoutes(const vector<Route> &routes, const vector<int> &demandas, int maxCapacity) {
    for (const auto &route : routes) {
        int requiredCapacity = calculateRouteCapacity(route.stations, demandas);

        if (requiredCapacity > maxCapacity) return false;
        if (route.stations.size() < 3) return false;
        if (route.stations.front() != route.stations.back()) return false;
    }
    
    return true;
}

void printAllRoutes(const vector<Route> &routes, const Instance *instance) {
    cout << (calculateTotalCost(routes, instance)) << endl;
    cout << routes.size() << endl;

    for (const auto &route : routes) {
        for (int i = 0; i < (int)route.stations.size(); i++) {
            cout << route.stations[i];

            if (i < (int)route.stations.size() - 1) cout << " ";
        }

        cout << endl;
    }
}

void Solver::Solve(Instance *instance) {
    srand(static_cast<unsigned int>(time(nullptr)));
    auto startTime = chrono::high_resolution_clock::now();

    vector<Route> bestSolutions;
    double bestCost     = numeric_limits<double>::max();
    double valorOtimo   = instance->valorOtimo;

    for (int i = 0; i < instance->repetitions; i++) {
        cout << "\n--- Iteração " << (i + 1) << " ---\n";

        // fase 1 - construção GRASP
        vector<Route> routes = GRASPConstruction(instance);

        // fase 2 - VND  
        vector<Route>   improvedRoutes  = VND(routes, instance);
        double          improvedCost    = calculateTotalCost(improvedRoutes, instance);

        // validar a solução antes de aceitar
        if (validateRoutes(improvedRoutes, instance->get_demandas(), instance->get_capVeiculos())) {
            if (improvedCost < bestCost) {
                bestCost        = improvedCost;
                bestSolutions   = improvedRoutes;
            }
        }
        
        // se encontrou o otimo, para execução
        if (abs(bestCost - valorOtimo) < 0.01) {
            break;
        }
    }

    auto endTime = chrono::high_resolution_clock::now();
    chrono::duration<double> timeDiff = endTime - startTime;

    cout << "\n=== RESULTADO FINAL GRASP ===\n";
    cout << fixed << setprecision(3);
    cout << "--------------------------------------------------------------------------\n";
    cout << "| Instancia    | Custo GRASP + VND | Valor Otimo |   Gap (%)   | Tempo (s) |\n";
    cout << "--------------------------------------------------------------------------\n";
    cout << "|"
              << setw(11) << instance->instanceName << " |"
              << setw(20) << bestCost << " |"
              << setw(12) << valorOtimo << " |"
              << setw(11) << (((bestCost - valorOtimo) / valorOtimo) * 100.0) << " |"
              << setw(10) << timeDiff.count() << " |\n";
    cout << "--------------------------------------------------------------------------\n";
    
    cout << "\n=== MELHOR SOLUÇÃO ENCONTRADA ===\n";
    printAllRoutes(bestSolutions, instance);
}

vector<Route> VND(vector<Route> routes, const Instance *instance) {
    vector<Vizinhança> neighborhoodOrder = {REINSERTION, TWO_OPT, SWAP};
    
    double currentCost = calculateTotalCost(routes, instance);
    
    int k = 0; // vizinhança atual
    while (k < QTD_VIZINHANCAS) {
        vector<Route> newRoutes = routes;

        switch (neighborhoodOrder[k]) {
            case SWAP:
                newRoutes = performSwap(newRoutes, instance);
                break;
            case TWO_OPT:
                newRoutes = performTwoOpt(newRoutes, instance);
                break;
            case REINSERTION:
                newRoutes = performReinsertion(newRoutes, instance);
                break;
        }

        double newCost = calculateTotalCost(newRoutes, instance);

        if (newCost < currentCost) {
            routes = newRoutes;
            currentCost = newCost;
            k = 0;
        } else k++;
    }

    return routes;
}

vector<Route> performSwap(vector<Route> routes, const Instance *instance) {
    const auto &distancias  = instance->get_distancias();
    const int maxCapacity   = instance->get_capVeiculos();

    // guardar o melhor movimento
    double bestDelta        = 0.0;
    bool foundImprovement   = false;
    
    vector<int> bestIntraStations;
    int bestIntraRouteIdx = -1;
    
    vector<int> bestInterR1Stations, bestInterR2Stations;
    int bestInterR1      = -1;
    int bestInterR2      = -1;
    bool isBestMoveIntra = true;

    // SWAP INTRA-ROTA
    for (int i = 0; i < (int)(routes.size()); i++) {
        if (routes[i].stations.size() <= 3) continue;

        for (int j = 1; j < (int)routes[i].stations.size() - 1; j++) {
            for (int k = j + 1; k < (int)routes[i].stations.size() - 1; k++) {
                
                vector<int> tempStations = routes[i].stations;
                swap(tempStations[j], tempStations[k]);

                double delta = 0;
                if (k == j + 1) { // Adjacentes
                    delta = - distancias[tempStations[j-1]][routes[i].stations[j]] 
                            - distancias[routes[i].stations[j]][routes[i].stations[k]] 
                            - distancias[routes[i].stations[k]][tempStations[k+1]]
                            + distancias[tempStations[j-1]][tempStations[j]] 
                            + distancias[tempStations[j]][tempStations[k]] 
                            + distancias[tempStations[k]][tempStations[k+1]];
                } else { // Não adjacentes
                    delta = - distancias[tempStations[j-1]][routes[i].stations[j]] 
                            - distancias[routes[i].stations[j]][tempStations[j+1]]
                            - distancias[tempStations[k-1]][routes[i].stations[k]] 
                            - distancias[routes[i].stations[k]][tempStations[k+1]]
                            + distancias[tempStations[j-1]][tempStations[j]] 
                            + distancias[tempStations[j]][tempStations[j+1]]
                            + distancias[tempStations[k-1]][tempStations[k]] 
                            + distancias[tempStations[k]][tempStations[k+1]];
                }

                // verifica se tem melhoria e a melhor até agora
                if (delta < bestDelta - 1e-9) {
                    int newCapacity = calculateRouteCapacity(tempStations, instance->get_demandas());
                    
                    if (newCapacity <= maxCapacity) {
                        bestDelta           = delta;
                        isBestMoveIntra     = true;
                        foundImprovement    = true;
                        bestIntraRouteIdx   = i;
                        bestIntraStations   = tempStations;
                    }
                }
            }
        }
    }

    // SWAP INTER-ROTAS
    for (int r1_idx = 0; r1_idx < (int)(routes.size()); r1_idx++) {
        for (int r2_idx = r1_idx + 1; r2_idx < (int)(routes.size()); r2_idx++) {
            for (int s1_idx = 1; s1_idx < (int)routes[r1_idx].stations.size() - 1; s1_idx++) {
                for (int s2_idx = 1; s2_idx < (int)routes[r2_idx].stations.size() - 1; s2_idx++) {
                    
                    const auto &stations1 = routes[r1_idx].stations;
                    const auto &stations2 = routes[r2_idx].stations;

                    double delta = - (distancias[stations1[s1_idx-1]][stations1[s1_idx]] + distancias[stations1[s1_idx]][stations1[s1_idx+1]])
                                   - (distancias[stations2[s2_idx-1]][stations2[s2_idx]] + distancias[stations2[s2_idx]][stations2[s2_idx+1]])
                                   + (distancias[stations1[s1_idx-1]][stations2[s2_idx]] + distancias[stations2[s2_idx]][stations1[s1_idx+1]])
                                   + (distancias[stations2[s2_idx-1]][stations1[s1_idx]] + distancias[stations1[s1_idx]][stations2[s2_idx+1]]);


                    if (delta < bestDelta - 1e-9) {
                        vector<int> r1_stations = stations1;
                        vector<int> r2_stations = stations2;
                        swap(r1_stations[s1_idx], r2_stations[s2_idx]);

                        int newCap1 = calculateRouteCapacity(r1_stations, instance->get_demandas());
                        int newCap2 = calculateRouteCapacity(r2_stations, instance->get_demandas());

                        if (newCap1 <= maxCapacity && newCap2 <= maxCapacity) {
                            foundImprovement    = true;
                            isBestMoveIntra     = false;
                            
                            bestDelta   = delta;
                            bestInterR1 = r1_idx;
                            bestInterR2 = r2_idx;
                            bestInterR1Stations = r1_stations;
                            bestInterR2Stations = r2_stations;
                        }
                    }
                }
            }
        }
    }

    // aplicar o melhor swap encontrado
    if (foundImprovement) {
        if (isBestMoveIntra) {  // melhor swap intra-rota
            routes[bestIntraRouteIdx].stations = bestIntraStations;
            routes[bestIntraRouteIdx].capacity = calculateRouteCapacity(bestIntraStations, instance->get_demandas());
        } else {                // melhor swap inter-rotas
            routes[bestInterR1].stations = bestInterR1Stations;
            routes[bestInterR1].capacity = calculateRouteCapacity(bestInterR1Stations, instance->get_demandas());
            routes[bestInterR2].stations = bestInterR2Stations;
            routes[bestInterR2].capacity = calculateRouteCapacity(bestInterR2Stations, instance->get_demandas());
        }
    }
    
    return routes;
}

vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance) {
    const auto &distancias  = instance->get_distancias();
    const int maxCapacity   = instance->get_capVeiculos();
    
    // guardar melhor movimento
    vector<int> bestStations;
    bool    foundImprovement = false;
    int     bestRouteIdx     = -1;
    double  bestDelta        = 0.0;

    // 2-OPT
    for (int i = 0; i < (int)(routes.size()); i++) {
        if (routes[i].stations.size() <= 4) continue; // poucas estações na rota

        for (int j = 1; j < (int)routes[i].stations.size() - 2; j++) {
            for (int k = j + 1; k < (int)routes[i].stations.size() - 1; k++) {
                const auto &stations = routes[i].stations;
                
                double delta = - distancias[stations[j-1]][stations[j]]
                               - distancias[stations[k]][stations[k+1]]
                               + distancias[stations[j-1]][stations[k]]
                               + distancias[stations[j]][stations[k+1]];

                // verifica se é melhoria
                if (delta < bestDelta - 1e-9) {
                    vector<int> newStations = stations;
                    reverse(newStations.begin() + j, newStations.begin() + k + 1);

                    int newCapacity = calculateRouteCapacity(newStations, instance->get_demandas());
                    if (newCapacity <= maxCapacity) {
                        foundImprovement = true;
                        bestDelta       = delta;
                        bestRouteIdx    = i;
                        bestStations    = newStations;
                    }
                }
            }
        }
    }

    // aplica o melhor 2-opt encontrado
    if (foundImprovement) {
        routes[bestRouteIdx].stations = bestStations;
        routes[bestRouteIdx].capacity = calculateRouteCapacity(bestStations, instance->get_demandas());
    }
    
    return routes;
}

vector<Route> performReinsertion(vector<Route> routes, const Instance *instance) {
    const auto &distancias  = instance->get_distancias();
    const int maxCapacity   = instance->get_capVeiculos();

    // guardar o melhor movimento
    bool foundImprovement   = false;
    bool isBestMoveIntra    = false;
    double bestDelta        = 0.0;
    
    // intra-rota
    int bestIntraRouteIdx = -1;
    vector<int> bestIntraStations;
    
    // inter-rotas
    vector<int> bestInterR1Stations, bestInterR2Stations;
    int bestInterR1             = -1;
    int bestInterR2             = -1;
    bool bestInterDeletesRoute  = false;

    // REINSERTION INTRA-ROTA
    for (int i = 0; i < (int)(routes.size()); i++) {
        if (routes[i].stations.size() <= 3) continue;

        for (int j = 1; j < (int)(routes[i].stations.size() - 1); j++) {
            for (int k = 1; k < (int)(routes[i].stations.size()); k++) {
                if (j == k || j == k - 1) continue;

                const auto &stations    = routes[i].stations;
                int stationToMove       = stations[j];

                double delta = - (distancias[stations[j-1]][stationToMove] + distancias[stationToMove][stations[j+1]])
                               - distancias[stations[k-1]][stations[k]]
                               + distancias[stations[j-1]][stations[j+1]]
                               + distancias[stations[k-1]][stationToMove] + distancias[stationToMove][stations[k]];

                // verifica se tem melhoria
                if (delta < bestDelta - 1e-9) {
                    vector<int> tempStations = stations;
                    tempStations.erase(tempStations.begin() + j);
                    tempStations.insert(tempStations.begin() + (k > j ? k - 1 : k), stationToMove);

                    int newCapacity = calculateRouteCapacity(tempStations, instance->get_demandas());
                    if (newCapacity <= maxCapacity) {
                        foundImprovement    = true;
                        isBestMoveIntra     = true;
                        bestIntraRouteIdx   = i;
                        bestIntraStations   = tempStations;
                        bestDelta           = delta;
                    }
                }
            }
        }
    }

    // REINSERTION INTER-ROTAS
    for (int r1_idx = 0; r1_idx < (int)(routes.size()); r1_idx++) {
        if (routes[r1_idx].stations.size() <= 2) continue;

        for (int r2_idx = 0; r2_idx < (int)(routes.size()); r2_idx++) {
            if (r1_idx == r2_idx) continue;

            for (int s1_idx = 1; s1_idx < (int)routes[r1_idx].stations.size() - 1; s1_idx++) {
                for (int s2_idx = 1; s2_idx < (int)routes[r2_idx].stations.size(); s2_idx++) {
                    
                    const auto &stations1   = routes[r1_idx].stations;
                    const auto &stations2   = routes[r2_idx].stations;

                    int stationToMove       = stations1[s1_idx];

                    double delta = - (distancias[stations1[s1_idx-1]][stationToMove] + distancias[stationToMove][stations1[s1_idx+1]])
                                   - distancias[stations2[s2_idx-1]][stations2[s2_idx]]
                                   + distancias[stations1[s1_idx-1]][stations1[s1_idx+1]]
                                   + distancias[stations2[s2_idx-1]][stationToMove] + distancias[stationToMove][stations2[s2_idx]];

                    // verifica melhoria
                    if (delta < bestDelta - 1e-9) {
                        vector<int> r1_stations = stations1;
                        vector<int> r2_stations = stations2;
                        
                        r1_stations.erase(r1_stations.begin() + s1_idx);
                        r2_stations.insert(r2_stations.begin() + s2_idx, stationToMove);

                        int newCap1 = (r1_stations.size() > 2) ? calculateRouteCapacity(r1_stations, instance->get_demandas()) : 0;
                        int newCap2 = calculateRouteCapacity(r2_stations, instance->get_demandas());

                        if (newCap1 <= maxCapacity && newCap2 <= maxCapacity) {
                            foundImprovement        = true;
                            isBestMoveIntra         = false;
                            bestDelta               = delta;
                            bestInterR1             = r1_idx;
                            bestInterR2             = r2_idx;
                            bestInterR1Stations     = r1_stations;
                            bestInterR2Stations     = r2_stations;
                            bestInterDeletesRoute   = (r1_stations.size() <= 2);
                        }
                    }
                }
            }
        }
    }

    // aplica mlelhor reinsertion encontrada
    if (foundImprovement) {
        if (isBestMoveIntra) { // intra-rota
            routes[bestIntraRouteIdx].stations = bestIntraStations;
            routes[bestIntraRouteIdx].capacity = calculateRouteCapacity(bestIntraStations, instance->get_demandas());
        } else { // inter-rotas
            routes[bestInterR2] = {bestInterR2Stations, calculateRouteCapacity(bestInterR2Stations, instance->get_demandas())};
            
            if (bestInterDeletesRoute)  routes.erase(routes.begin() + bestInterR1);
            else                        routes[bestInterR1] = {bestInterR1Stations, calculateRouteCapacity(bestInterR1Stations, instance->get_demandas())};
        }
    }
    
    return routes;
}

vector<Route> GRASPConstruction(const Instance *instance) {
    const vector<int> &demandas             = instance->get_demandas();
    const vector<vector<int>> &distancias   = instance->get_distancias();
    
    double alpha = instance->alpha;
    
    int maxCapacity     = instance->get_capVeiculos();
    int totalStations   = instance->get_qtdEstacoes();
    int unvisitedCount  = totalStations;

    vector<bool> visitedStations(totalStations + 1, false);
    vector<Route> routes;
    routes.reserve(max(1, totalStations / 5));


    while (unvisitedCount > 0) {
        // seed inicial = estação mais próxima do depósito
        int seed        = -1;
        int bestDepDist = numeric_limits<int>::max();
        
        for (int s = 1; s <= totalStations; ++s) {
            if (visitedStations[s]) continue;
            
            int d0s = distancias[0][s];
            if (d0s < bestDepDist) { // acha a estação mais próxima do depósito
                bestDepDist = d0s;
                seed        = s;
            }
        }

        if (seed == -1) break;

        // cria nova rota começando com seed (a mais proxima do depósito)
        vector<int> routeStations = {0, seed, 0};
        
        visitedStations[seed]   = true;
        unvisitedCount--;

        bool extended = true;
        while (extended && unvisitedCount > 0) {
            extended = false;
            int last = routeStations[routeStations.size() - 2]; // último nó antes do depósito

            // possiveis candidatos viáveis
            vector<pair<int,int>> candidates; // (station, dist)
            int dmin = numeric_limits<int>::max();
            int dmax = 0;

            for (int station = 1; station <= totalStations; station++) {
                if (visitedStations[station]) continue;

                int distance = distancias[last][station];

                // simula adicionar estação station
                vector<int> temp = routeStations;
                temp.insert(temp.end() - 1, station);

                int requiredCap = calculateRouteCapacity(temp, demandas); // valida a capacidade da rota
                if (requiredCap <= maxCapacity) {
                    candidates.push_back({station, distance});
                    dmin = min(dmin, distance);
                    dmax = max(dmax, distance);
                }
            }

            if (!candidates.empty()) {
                // monta RCL - lista restrita de candidatos
                vector<int> RCL;
                double limit = dmin + (alpha * (dmax - dmin)); // limit representa o custo máximo aceitável na RCL
                for (auto &c : candidates) {
                    // (s, d)
                    if (c.second <= limit) RCL.push_back(c.first);
                }

                if (!RCL.empty()) {
                    // escolhe aleatoriamente da RCL
                    int chosen = RCL[rand() % RCL.size()];
                    
                    routeStations.insert(routeStations.end() - 1, chosen); // insere antes do depósito
                    unvisitedCount--;
                    
                    visitedStations[chosen] = true;
                    extended                = true; // conseguiu estender a rota, continua no loop
                }
            }
        }

        int finalCap = calculateRouteCapacity(routeStations, demandas);
        routes.push_back({routeStations, finalCap});
    }

    return routes;
}
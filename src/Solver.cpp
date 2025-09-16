#include "Solver.hpp"
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <numeric>

#define QTD_VIZINHANCAS 3

typedef enum {
    SWAP = 1,
    TWO_OPT,
    REINSERTION
} Vizinhança;

struct Saving {
    double value;
    int fromStation;
    int toStation;
};

vector<Route> VND(vector<Route> routes, const Instance *instance);
vector<Route> performSwap(vector<Route> routes, const Instance *instance);
vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance);
vector<Route> performReinsertion(vector<Route> routes, const Instance *instance);

vector<Route> GRASPConstruction(const Instance *instance);

double calculateTotalCost(const vector<Route> &routes, const Instance *instance) {
    if (routes.empty()) return 0.0;
    
    const auto &distancias = instance->get_distancias(); // cache
    double totalCost = 0.0;
    
    for (const auto &route : routes) {
        const auto &stations = route.stations; // cache
        const size_t stationsSize = stations.size();
        
        if (stationsSize == 3) { // [0, station, 0]
            totalCost += distancias[stations[0]][stations[1]] + distancias[stations[1]][stations[2]];
        } else if (stationsSize == 4) { // [0, s1, s2, 0]
            totalCost += distancias[stations[0]][stations[1]] 
                       + distancias[stations[1]][stations[2]]
                       + distancias[stations[2]][stations[3]];
        } else {
            // rotas maioresclear
            for (size_t i = 0; i < stationsSize - 1; ++i) {
                totalCost += distancias[stations[i]][stations[i + 1]];
            }
        }
    }
    return totalCost;
}

int calculateRouteCapacity(const std::vector<int>& stations,
                           const std::vector<int>& demandas) {
    if (stations.size() <= 2) return 0; // rota só com depósito

    int currentLoad = 0;
    int minPrefix = 0;
    int maxPrefix = 0;

    // percorre apenas as estações (ignora depósito no início e fim)
    for (size_t i = 1; i < stations.size() - 1; i++) {
        int node = stations[i];
        int demand = demandas[node - 1]; // atenção: node começa em 1

        currentLoad += demand;

        if (currentLoad < minPrefix) minPrefix = currentLoad;
        if (currentLoad > maxPrefix) maxPrefix = currentLoad;
    }

    int initialLoad = (minPrefix < 0) ? -minPrefix : 0;
    int maxLoad = initialLoad + maxPrefix;

    return maxLoad;
}

bool validateRoutesQuiet(const std::vector<Route>& routes,
                         const std::vector<int>& demandas,
                         int maxCapacity) {
    for (const auto& route : routes) {
        int requiredCapacity = calculateRouteCapacity(route.stations, demandas);
        if (requiredCapacity > maxCapacity) {
            return false;
        }
    }
    return true;
}

void printAllRoutes(const vector<Route> &routes, const Instance *instance) {
    // Formato especificado: valor da solução, número de veículos, depois cada rota
    cout << static_cast<int>(calculateTotalCost(routes, instance)) << endl;
    cout << routes.size() << endl;
    
    for (const auto& route : routes) {
        for (size_t i = 0; i < route.stations.size(); i++) {
            cout << route.stations[i];
            if (i < route.stations.size() - 1) cout << " ";
        }
        cout << endl;
    }
}

void Solver::Solve(Instance *instance) {
    // Inicializar semente aleatória
    srand(static_cast<unsigned int>(time(nullptr)));

    vector<Route> bestSolutions;
    double bestCost = numeric_limits<double>::max();
    double valorOtimo = instance->valorOtimo;

    for (int i = 0; i < instance->repetitions; i++) {
        cout << "\n--- Iteração " << (i + 1) << " ---\n";

        // fase 1 - construção GRASP
        vector<Route> routes = GRASPConstruction(instance);

        // fase 2 - VND  
        vector<Route> improvedRoutes = VND(routes, instance);
        double improvedCost = calculateTotalCost(improvedRoutes, instance);

        // validar a solução antes de aceitar
        if (validateRoutesQuiet(improvedRoutes, instance->get_demandas(), instance->get_capVeiculos())) {
            if (improvedCost < bestCost) {
                bestCost = improvedCost;
                bestSolutions = improvedRoutes;
            }
        }
        
        // se encontrou o otimo, para execução
        if (abs(bestCost - valorOtimo) < 0.01) {
            break;
        }
    }

    // Tabela final de comparação, agora com a melhor solução do GRASP
    cout << "\n=== RESULTADO FINAL GRASP ===\n";
    cout << fixed << setprecision(2);
    cout << "--------------------------------------------------------------------------\n";
    cout << "| Instancia | Custo GRASP + VND | Valor Otimo |   Gap (%)   |\n";
    cout << "--------------------------------------------------------------------------\n";
    cout << "|"
              << setw(10) << instance->instanceName << " |"
              << setw(22) << bestCost << " |"
              << setw(12) << valorOtimo << " |"
              << setw(11) << (((bestCost - valorOtimo) / valorOtimo) * 100.0) << " |\n";
    cout << "--------------------------------------------------------------------------\n";
    
    // Opcionalmente, você pode adicionar a melhor rota final
    cout << "\n=== MELHOR SOLUÇÃO ENCONTRADA ===\n";
    printAllRoutes(bestSolutions, instance);
}

vector<Route> VND(vector<Route> routes, const Instance *instance) {
    vector<Vizinhança> neighborhoodOrder = {SWAP, TWO_OPT, REINSERTION};
    
    int k = 0; // vizinhança atual
    double currentCost = calculateTotalCost(routes, instance);

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
        } else {
            k++;
        }
    }

    return routes;
}

vector<Route> performSwap(vector<Route> routes, const Instance *instance) {
    const auto& distancias = instance->get_distancias();
    const int maxCapacity = instance->get_capVeiculos();

    // Variáveis para memorizar o melhor movimento
    double bestDelta = 0.0;
    bool foundImprovement = false;
    
    // Para movimentos intra-rota
    int bestIntraRouteIdx = -1;
    size_t bestIntraJ = 0, bestIntraK = 0;
    vector<int> bestIntraStations;
    int bestIntraCapacity = 0;
    
    // Para movimentos inter-rotas
    int bestInterR1 = -1, bestInterR2 = -1;
    size_t bestInterS1 = 0, bestInterS2 = 0;
    vector<int> bestInterR1Stations, bestInterR2Stations;
    int bestInterCap1 = 0, bestInterCap2 = 0;
    bool isBestMoveIntra = true;

    // AVALIAÇÃO DE TODOS OS SWAPS INTRA-ROTA
    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        if (routes[i].stations.size() <= 3) continue;

        for (size_t j = 1; j < routes[i].stations.size() - 1; j++) {
            for (size_t k = j + 1; k < routes[i].stations.size() - 1; k++) {
                
                vector<int> tempStations = routes[i].stations;
                swap(tempStations[j], tempStations[k]);

                double delta;
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

                // Verificar se é uma melhoria válida e a melhor até agora
                if (delta < bestDelta - 1e-9) {
                    int newCapacity = calculateRouteCapacity(tempStations, instance->get_demandas());
                    if (newCapacity <= maxCapacity) {
                        bestDelta = delta;
                        foundImprovement = true;
                        isBestMoveIntra = true;
                        bestIntraRouteIdx = i;
                        bestIntraJ = j;
                        bestIntraK = k;
                        bestIntraStations = tempStations;
                        bestIntraCapacity = newCapacity;
                    }
                }
            }
        }
    }

    // AVALIAÇÃO DE TODOS OS SWAPS INTER-ROTAS
    for (int r1_idx = 0; r1_idx < static_cast<int>(routes.size()); r1_idx++) {
        for (int r2_idx = r1_idx + 1; r2_idx < static_cast<int>(routes.size()); r2_idx++) {
            for (size_t s1_idx = 1; s1_idx < routes[r1_idx].stations.size() - 1; s1_idx++) {
                for (size_t s2_idx = 1; s2_idx < routes[r2_idx].stations.size() - 1; s2_idx++) {
                    
                    const auto& stations1 = routes[r1_idx].stations;
                    const auto& stations2 = routes[r2_idx].stations;

                    double delta = - (distancias[stations1[s1_idx-1]][stations1[s1_idx]] + distancias[stations1[s1_idx]][stations1[s1_idx+1]])
                                   - (distancias[stations2[s2_idx-1]][stations2[s2_idx]] + distancias[stations2[s2_idx]][stations2[s2_idx+1]])
                                   + (distancias[stations1[s1_idx-1]][stations2[s2_idx]] + distancias[stations2[s2_idx]][stations1[s1_idx+1]])
                                   + (distancias[stations2[s2_idx-1]][stations1[s1_idx]] + distancias[stations1[s1_idx]][stations2[s2_idx+1]]);

                    // Verificar se é uma melhoria válida e a melhor até agora
                    if (delta < bestDelta - 1e-9) {
                        vector<int> r1_stations = stations1;
                        vector<int> r2_stations = stations2;
                        swap(r1_stations[s1_idx], r2_stations[s2_idx]);

                        int newCap1 = calculateRouteCapacity(r1_stations,   instance->get_demandas());
                        int newCap2 = calculateRouteCapacity(r2_stations, instance->get_demandas());

                        if (newCap1 <= maxCapacity && newCap2 <= maxCapacity) {
                            bestDelta = delta;
                            foundImprovement = true;
                            isBestMoveIntra = false;
                            bestInterR1 = r1_idx;
                            bestInterR2 = r2_idx;
                            bestInterS1 = s1_idx;
                            bestInterS2 = s2_idx;
                            bestInterR1Stations = r1_stations;
                            bestInterR2Stations = r2_stations;
                            bestInterCap1 = newCap1;
                            bestInterCap2 = newCap2;
                        }
                    }
                }
            }
        }
    }

    // APLICAR O MELHOR MOVIMENTO ENCONTRADO
    if (foundImprovement) {
        if (isBestMoveIntra) {
            // Aplicar o melhor swap intra-rota
            routes[bestIntraRouteIdx].stations = bestIntraStations;
            routes[bestIntraRouteIdx].capacity = calculateRouteCapacity(bestIntraStations, instance->get_demandas());
        } else {
            // Aplicar o melhor swap inter-rotas
            routes[bestInterR1].stations = bestInterR1Stations;
            routes[bestInterR1].capacity = calculateRouteCapacity(bestInterR1Stations, instance->get_demandas());
            routes[bestInterR2].stations = bestInterR2Stations;
            routes[bestInterR2].capacity = calculateRouteCapacity(bestInterR2Stations, instance->get_demandas());
        }
    }
    
    return routes;
}

vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance) {
    const auto& distancias = instance->get_distancias();
    const int maxCapacity = instance->get_capVeiculos();
    
    // Variáveis para memorizar o melhor movimento
    bool foundImprovement = false;
    double bestDelta = 0.0;
    int bestRouteIdx = -1;
    vector<int> bestStations;
    int bestCapacity = 0;

    // AVALIAR TODOS OS MOVIMENTOS 2-OPT POSSÍVEIS
    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        if (routes[i].stations.size() <= 4) continue;

        for (size_t j = 1; j < routes[i].stations.size() - 2; j++) {
            for (size_t k = j + 1; k < routes[i].stations.size() - 1; k++) {
                
                const auto& stations = routes[i].stations;
                double delta = - distancias[stations[j-1]][stations[j]]
                               - distancias[stations[k]][stations[k+1]]
                               + distancias[stations[j-1]][stations[k]]
                               + distancias[stations[j]][stations[k+1]];

                // Verificar se é uma melhoria válida e a melhor até agora
                if (delta < bestDelta - 1e-9) {
                    vector<int> newStations = stations;
                    reverse(newStations.begin() + j, newStations.begin() + k + 1);

                    int newCapacity = calculateRouteCapacity(newStations, instance->get_demandas());
                    if (newCapacity <= maxCapacity) {
                        bestDelta = delta;
                        foundImprovement = true;
                        bestRouteIdx = i;
                        bestStations = newStations;
                        bestCapacity = newCapacity;
                    }
                }
            }
        }
    }

    // APLICAR O MELHOR MOVIMENTO ENCONTRADO
    if (foundImprovement) {
        routes[bestRouteIdx].stations = bestStations;
        routes[bestRouteIdx].capacity = calculateRouteCapacity(bestStations, instance->get_demandas());
    }
    
    return routes;
}

vector<Route> performReinsertion(vector<Route> routes, const Instance *instance) {
    const auto& distancias = instance->get_distancias();
    const int maxCapacity = instance->get_capVeiculos();

    // Variáveis para memorizar o melhor movimento
    bool foundImprovement = false;
    double bestDelta = 0.0;
    bool isBestMoveIntra = false;
    
    // Para movimentos intra-rota
    int bestIntraRouteIdx = -1;
    vector<int> bestIntraStations;
    int bestIntraCapacity = 0;
    
    // Para movimentos inter-rotas
    int bestInterR1 = -1, bestInterR2 = -1;
    vector<int> bestInterR1Stations, bestInterR2Stations;
    int bestInterCap1 = 0, bestInterCap2 = 0;
    bool bestInterDeletesRoute = false;

    // AVALIAR REINSERTION INTRA-ROTA
    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        if (routes[i].stations.size() <= 3) continue;

        for (size_t j = 1; j < routes[i].stations.size() - 1; j++) {
            for (size_t k = 1; k < routes[i].stations.size(); k++) {
                if (j == k || j == k - 1) continue;

                const auto& stations = routes[i].stations;
                int stationToMove = stations[j];

                double delta = - (distancias[stations[j-1]][stationToMove] + distancias[stationToMove][stations[j+1]])
                               - distancias[stations[k-1]][stations[k]]
                               + distancias[stations[j-1]][stations[j+1]]
                               + distancias[stations[k-1]][stationToMove] + distancias[stationToMove][stations[k]];

                // Verificar se é uma melhoria válida e a melhor até agora
                if (delta < bestDelta - 1e-9) {
                    vector<int> tempStations = stations;
                    tempStations.erase(tempStations.begin() + j);
                    tempStations.insert(tempStations.begin() + (k > j ? k - 1 : k), stationToMove);

                    int newCapacity = calculateRouteCapacity(tempStations, instance->get_demandas());
                    if (newCapacity <= maxCapacity) {
                        bestDelta = delta;
                        foundImprovement = true;
                        isBestMoveIntra = true;
                        bestIntraRouteIdx = i;
                        bestIntraStations = tempStations;
                        bestIntraCapacity = newCapacity;
                    }
                }
            }
        }
    }

    // AVALIAR REINSERTION INTER-ROTAS
    for (int r1_idx = 0; r1_idx < static_cast<int>(routes.size()); r1_idx++) {
        if (routes[r1_idx].stations.size() <= 2) continue;

        for (int r2_idx = 0; r2_idx < static_cast<int>(routes.size()); r2_idx++) {
            if (r1_idx == r2_idx) continue;

            for (size_t s1_idx = 1; s1_idx < routes[r1_idx].stations.size() - 1; s1_idx++) {
                for (size_t s2_idx = 1; s2_idx < routes[r2_idx].stations.size(); s2_idx++) {
                    
                    const auto& stations1 = routes[r1_idx].stations;
                    const auto& stations2 = routes[r2_idx].stations;
                    int stationToMove = stations1[s1_idx];

                    double delta = - (distancias[stations1[s1_idx-1]][stationToMove] + distancias[stationToMove][stations1[s1_idx+1]])
                                   - distancias[stations2[s2_idx-1]][stations2[s2_idx]]
                                   + distancias[stations1[s1_idx-1]][stations1[s1_idx+1]]
                                   + distancias[stations2[s2_idx-1]][stationToMove] + distancias[stationToMove][stations2[s2_idx]];

                    // Verificar se é uma melhoria válida e a melhor até agora
                    if (delta < bestDelta - 1e-9) {
                        vector<int> r1_stations = stations1;
                        vector<int> r2_stations = stations2;
                        
                        r1_stations.erase(r1_stations.begin() + s1_idx);
                        r2_stations.insert(r2_stations.begin() + s2_idx, stationToMove);

                        int newCap1 = (r1_stations.size() > 2) ? calculateRouteCapacity(r1_stations, instance->get_demandas()) : 0;
                        int newCap2 = calculateRouteCapacity(r2_stations, instance->get_demandas());

                        if (newCap1 <= maxCapacity && newCap2 <= maxCapacity) {
                            bestDelta = delta;
                            foundImprovement = true;
                            isBestMoveIntra = false;
                            bestInterR1 = r1_idx;
                            bestInterR2 = r2_idx;
                            bestInterR1Stations = r1_stations;
                            bestInterR2Stations = r2_stations;
                            bestInterCap1 = newCap1;
                            bestInterCap2 = newCap2;
                            bestInterDeletesRoute = (r1_stations.size() <= 2);
                        }
                    }
                }
            }
        }
    }

    // APLICAR O MELHOR MOVIMENTO ENCONTRADO
    if (foundImprovement) {
        if (isBestMoveIntra) {
            // Aplicar o melhor movimento intra-rota
            routes[bestIntraRouteIdx].stations = bestIntraStations;
            routes[bestIntraRouteIdx].capacity = calculateRouteCapacity(bestIntraStations, instance->get_demandas());
        } else {
            // Aplicar o melhor movimento inter-rotas
            routes[bestInterR2] = {bestInterR2Stations, calculateRouteCapacity(bestInterR2Stations, instance->get_demandas())};
            if (bestInterDeletesRoute) {
                routes.erase(routes.begin() + bestInterR1);
            } else {
                routes[bestInterR1] = {bestInterR1Stations, calculateRouteCapacity(bestInterR1Stations, instance->get_demandas())};
            }
        }
    }
    
    return routes;
}

// Estrutura para representar uma possível inserção
struct Insertion {
    double cost;
    int station;
    int route_idx;
    int position;

    bool operator<(const Insertion &other) const {
        return cost < other.cost;
    }
};

vector<Route> GRASPConstruction(const Instance *instance) {
    // IMPLEMENTAÇÃO: Nearest Neighbor construction (substitui GRASP)
    const vector<int>& demandas = instance->get_demandas(); 
    const vector<vector<int>>& distancias = instance->get_distancias();
    
    int maxCapacity = instance->get_capVeiculos();
    int totalStations = instance->get_qtdEstacoes();
    
    vector<bool> visitedStations(totalStations + 1, false);
    int unvisitedCount = totalStations;
    vector<Route> routes;
    routes.reserve( std::max(1, totalStations / 5) );

    // Precompute distances from depot for seeding
    while (unvisitedCount > 0) {
        // find seed: nearest unvisited station to depot (0)
        int seed = -1;
        int bestDepDist = std::numeric_limits<int>::max();
        for (int s = 1; s <= totalStations; ++s) {
            if (visitedStations[s]) continue;
            int d0s = distancias[0][s];
            if (d0s < bestDepDist) {
                bestDepDist = d0s;
                seed = s;
            }
        }
        if (seed == -1) break; // safety

        // create new route starting with seed
        vector<int> routeStations;
        routeStations.push_back(0);
        routeStations.push_back(seed);
        routeStations.push_back(0);
        visitedStations[seed] = true;
        unvisitedCount--;

        // Update capacity for this route
        int currCap = calculateRouteCapacity(routeStations, instance->get_demandas());

        bool extended = true;
        while (extended && unvisitedCount > 0) {
            extended = false;
            // last real node before depot
            int last = routeStations[routeStations.size() - 2];

            // find nearest feasible unvisited station to 'last'
            int bestStation = -1;
            int bestDist = std::numeric_limits<int>::max();
            for (int s = 1; s <= totalStations; ++s) {
                if (visitedStations[s]) continue;
                int d = distancias[last][s];
                if (d >= bestDist) continue;

                // simulate appending s at end (before depot)
                vector<int> temp = routeStations;
                temp.insert(temp.end() - 1, s); // before final 0
                int requiredCap = calculateRouteCapacity(temp, instance->get_demandas());
                if (requiredCap <= maxCapacity) {
                    bestDist = d;
                    bestStation = s;
                }
            }

            if (bestStation != -1) {
                // append bestStation
                routeStations.insert(routeStations.end() - 1, bestStation);
                visitedStations[bestStation] = true;
                unvisitedCount--;
                currCap = calculateRouteCapacity(routeStations, instance->get_demandas());
                extended = true;
            }
        }

        // finalize route
        int finalCap = calculateRouteCapacity(routeStations, instance->get_demandas());
        routes.push_back({routeStations, finalCap});
    }

    return routes;
}
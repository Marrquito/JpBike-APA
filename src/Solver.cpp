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
    
    const auto &distancias = instance->get_distancias(); // Cache da referência
    double totalCost = 0.0;
    
    for (const auto &route : routes) {
        const auto &stations = route.stations; // Cache da referência
        const size_t stationsSize = stations.size();
        
        // Loop desenrolado para as rotas pequenas (mais comuns)
        if (stationsSize == 3) { // [0, station, 0]
            totalCost += distancias[stations[0]][stations[1]] + distancias[stations[1]][stations[2]];
        } else if (stationsSize == 4) { // [0, s1, s2, 0]
            totalCost += distancias[stations[0]][stations[1]] 
                       + distancias[stations[1]][stations[2]]
                       + distancias[stations[2]][stations[3]];
        } else {
            // Loop otimizado para rotas maiores
            for (size_t i = 0; i < stationsSize - 1; ++i) {
                totalCost += distancias[stations[i]][stations[i + 1]];
            }
        }
    }
    return totalCost;
}

bool validateRoutesQuiet(const vector<Route> &routes, const Instance *instance) {
    vector<bool> visitedStations(instance->get_qtdEstacoes() + 1, false);
    visitedStations[0] = true; // deposito
    vector<int> demandas = instance->get_demandas();

    for (int rota = 0; rota < static_cast<int>(routes.size()); rota++) {
        const Route &route = routes[rota];
        
        // verifica se a rota começa e termina no depósito
        if (route.stations.empty() || route.stations[0] != 0 || route.stations.back() != 0) {
            return false;
        }
        
        // Verificar se há estações suficientes na rota
        if (route.stations.size() < 3) {
            return false;
        }
        
        //  carga inicial necessária (total de demandas negativas na rota)
        int initialLoad = 0;
        for (int i = 1; i < static_cast<int>(route.stations.size()) - 1; i++) {
            int station = route.stations[i];
            if (demandas[station - 1] < 0) { // -1 porque demandas não inclui depósito
                initialLoad += abs(demandas[station - 1]); // demandas negativas
            }
        }
        
        int currentLoad = initialLoad; // Começar com a carga necessária
        
        // rota estação por estação (- depósito inicial e final)
        for (int i = 1; i < static_cast<int>(route.stations.size()) - 1; i++) {
            int station = route.stations[i];
            
            // validando estação
            if (station == 0) {
                return false;
            }
            
            // validando unica visita na estação
            if (visitedStations[station]) {
                return false;
            }
            visitedStations[station] = true;
            
            // aplicar demanda da estação
            currentLoad += demandas[station - 1]; // -1 porque demandas não inclui depósito
            
            // Verificar limites de capacidade
            if (currentLoad < 0) {
                return false;
            }
            
            if (currentLoad > instance->get_capVeiculos()) {
                return false;
            }
        }
    }
    
    // todas estações visitadas
    for (int i = 1; i <= instance->get_qtdEstacoes(); i++) {
        if (!visitedStations[i]) {
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
        if (validateRoutesQuiet(improvedRoutes, instance)) {
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
    cout << "| Instancia | Custo GRASP + VND | Valor Otimo |   Gap (%)   |\n";
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

int calculateRouteCapacity(const vector<int> &stations, const Instance *instance) {
    if (stations.size() <= 2) return 0;
    
    const vector<int> &demandas = instance->get_demandas();
    
    // 1. Calcular a carga inicial necessária (soma de todas as coletas)
    int initialLoad = 0;
    for (size_t i = 1; i < stations.size() - 1; ++i) {
        int demand = demandas[stations[i] - 1];
        if (demand < 0) {
            initialLoad -= demand; // Somar o valor absoluto da demanda negativa
        }
    }
    
    // 2. Simular a rota para encontrar a carga máxima transportada
    int currentLoad = initialLoad;
    int maxLoad = initialLoad;
    
    for (size_t i = 1; i < stations.size() - 1; ++i) {
        // Aplica a demanda da estação (positiva ou negativa)
        currentLoad += demandas[stations[i] - 1];
        if (currentLoad > maxLoad) {
            maxLoad = currentLoad;
        }
    }
    
    return maxLoad;
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
                    int newCapacity = calculateRouteCapacity(tempStations, instance);
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

                        int newCap1 = calculateRouteCapacity(r1_stations, instance);
                        int newCap2 = calculateRouteCapacity(r2_stations, instance);

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
            routes[bestIntraRouteIdx].capacity = bestIntraCapacity;
        } else {
            // Aplicar o melhor swap inter-rotas
            routes[bestInterR1].stations = bestInterR1Stations;
            routes[bestInterR1].capacity = bestInterCap1;
            routes[bestInterR2].stations = bestInterR2Stations;
            routes[bestInterR2].capacity = bestInterCap2;
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
                    
                    int newCapacity = calculateRouteCapacity(newStations, instance);
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
        routes[bestRouteIdx].capacity = bestCapacity;
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
                    
                    int newCapacity = calculateRouteCapacity(tempStations, instance);
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

                        int newCap1 = (r1_stations.size() > 2) ? calculateRouteCapacity(r1_stations, instance) : 0;
                        int newCap2 = calculateRouteCapacity(r2_stations, instance);

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
            routes[bestIntraRouteIdx].capacity = bestIntraCapacity;
        } else {
            // Aplicar o melhor movimento inter-rotas
            routes[bestInterR2] = {bestInterR2Stations, bestInterCap2};
            if (bestInterDeletesRoute) {
                routes.erase(routes.begin() + bestInterR1);
            } else {
                routes[bestInterR1] = {bestInterR1Stations, bestInterCap1};
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
    vector<Route> routes;
    const vector<int>& demandas = instance->get_demandas(); 
    const vector<vector<int>>& distancias = instance->get_distancias();
    
    float alpha = instance->alpha;
    int maxCapacity = instance->get_capVeiculos();
    int totalStations = instance->get_qtdEstacoes();
    
    vector<bool> visitedStations(totalStations + 1, false);
    int unvisitedCount = totalStations;

    while (unvisitedCount > 0) {
        vector<Insertion> candidates;

        // 1. Gerar todos os candidatos de inserção com heurística melhorada
        for (int station = 1; station <= totalStations; station++) {
            if (visitedStations[station]) continue;

            int stationDemand = abs(demandas[station - 1]);

            // a) Inserir em rotas existentes (com prioridade para rotas com espaço)
            for (int r_idx = 0; r_idx < static_cast<int>(routes.size()); r_idx++) {
                // Checar capacidade primeiro
                if (routes[r_idx].capacity + stationDemand > maxCapacity) continue;

                // Calcular porcentagem de capacidade utilizada
                double capacityUsage = static_cast<double>(routes[r_idx].capacity + stationDemand) / maxCapacity;
                
                for (int pos = 1; pos < static_cast<int>(routes[r_idx].stations.size()); pos++) {
                    int prev = routes[r_idx].stations[pos - 1];
                    int next = routes[r_idx].stations[pos];
                    double insertionCost = distancias[prev][station] + distancias[station][next] - distancias[prev][next];
                    
                    // Penalizar levemente inserções que deixam rota pouco utilizada
                    if (capacityUsage < 0.6) {
                        insertionCost += insertionCost * 0.1; // 10% de penalidade
                    }
                    
                    candidates.push_back({insertionCost, station, r_idx, pos});
                }
            }

            // b) Criar uma nova rota para a estação (com penalidade maior para desencorajar rotas desnecessárias)
            double newRouteCost = distancias[0][station] + distancias[station][0];
            newRouteCost *= 1.2; // 20% de penalidade para novas rotas
            candidates.push_back({newRouteCost, station, -1, -1});
        }

        if (candidates.empty()) break;

        // 2. Construir a RCL (Lista Restrita de Candidatos) mais restritiva
        sort(candidates.begin(), candidates.end());

        double minCost = candidates.front().cost;
        double maxCost = candidates.back().cost;
        double threshold = minCost + alpha * (maxCost - minCost);

        vector<Insertion> rcl;
        for (const auto& candidate : candidates) {
            if (candidate.cost <= threshold) {
                rcl.push_back(candidate);
            }
        }

        // Limitar o tamanho da RCL para melhor qualidade
        if (rcl.size() > 5) {
            rcl.resize(5);
        }

        if (rcl.empty()) {
            rcl.push_back(candidates.front());
        }

        // 3. Selecionar da RCL com leve preferência pelos melhores candidatos
        int chosenIndex;
        if (rcl.size() == 1) {
            chosenIndex = 0;
        } else {
            // 60% de chance para os melhores 50% da RCL
            if ((rand() % 100) < 60) {
                chosenIndex = rand() % ((rcl.size() + 1) / 2);
            } else {
                chosenIndex = rand() % rcl.size();
            }
        }
        
        const Insertion& chosen = rcl[chosenIndex];

        if (chosen.route_idx == -1) {
            // Criar nova rota
            vector<int> newRouteStations = {0, chosen.station, 0};
            int newCapacity = calculateRouteCapacity(newRouteStations, instance);
            routes.push_back({newRouteStations, newCapacity});
        } else {
            // Inserir em rota existente
            routes[chosen.route_idx].stations.insert(routes[chosen.route_idx].stations.begin() + chosen.position, chosen.station);
            routes[chosen.route_idx].capacity = calculateRouteCapacity(routes[chosen.route_idx].stations, instance);
        }

        visitedStations[chosen.station] = true;
        unvisitedCount--;
    }
    
    return routes;
}
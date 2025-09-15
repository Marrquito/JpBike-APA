#include "Solver.hpp"
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <ctime>

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

void printRoute(const Route &route) {
    cout << "[";
    for (int i = 0; i < static_cast<int>(route.stations.size()); i++) {
        cout << route.stations[i];
        if (i < static_cast<int>(route.stations.size()) - 1) cout << " -> ";
    }
    cout << "] capacidade: " << route.capacity << ")\n";
}

void printAllRoutes(const vector<Route> &routes, const Instance *instance) {
    cout << "=== ROTAS ENCONTRADAS ===\n";
    cout << "número total de rotas: " << routes.size() << "\n";
    cout << "custo total: " << calculateTotalCost(routes, instance) << "\n\n";

    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        cout << "rota " << (i + 1) << ": ";
        printRoute(routes[i]);
    }
    cout << endl;
}

void Solver::Solve(Instance *instance) {
    cout << "\n=== Iniciando GRASP ===\n";
    
    // Inicializar semente aleatória
    srand(static_cast<unsigned int>(time(nullptr)));

    vector<Route> bestSolutions;
    double bestCost = numeric_limits<double>::max();
    double valorOtimo = instance->valorOtimo;

    for (int i = 0; i < instance->repetitions; i++) {
        cout << "\n--- Iteração " << (i + 1) << " ---\n";

        // fase 1 - construção GRASP
        vector<Route> routes = GRASPConstruction(instance);
        double constructionCost = calculateTotalCost(routes, instance);
        cout << "Custo após construção: " << constructionCost << "\n";

        // fase 2 - VND  
        vector<Route> improvedRoutes = VND(routes, instance);
        double improvedCost = calculateTotalCost(improvedRoutes, instance);
        cout << "Custo após VND: " << improvedCost << "\n";

        // validar a solução antes de aceitar
        if (validateRoutesQuiet(improvedRoutes, instance)) {
            if (improvedCost < bestCost) {
                bestCost = improvedCost;
                bestSolutions = improvedRoutes;
                cout << "*** NOVA MELHOR SOLUÇÃO VÁLIDA! (Gap: " 
                     << fixed << setprecision(2) 
                     << (((bestCost - valorOtimo) / valorOtimo) * 100.0) << "%) ***\n";
            } 
        }  
        // se encontrou o otimo, para execução
        if (abs(bestCost - valorOtimo) < 0.01) {
            cout << "*** ÓTIMO ENCONTRADO! Parando execução. ***\n";
            break;
        }
    }

    // Tabela final de comparação, agora com a melhor solução do GRASP
    cout << "\n=== RESULTADO FINAL GRASP ===\n";
    cout << fixed << setprecision(2);
    cout << "--------------------------------------------------------------------------\n";
    cout << "| Instancia | Custo GRASP + VND     | Valor Otimo |   Gap (%)   |\n";
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

// comparar economias em ordem decrescente
bool compareSavings(const Saving &a, const Saving &b) {
    return a.value > b.value;
}

int calculateRouteCapacity(const vector<int> &stations, const Instance *instance) {
    if (stations.size() <= 2) return 0;
    
    const vector<int> &demandas = instance->get_demandas();
    int maxLoad = 0;
    int currentLoad = 0;
    
    //carga inicial necessária (somatório de demandas negativas)
    for (size_t i = 1; i < stations.size() - 1; ++i) {
        int demand = demandas[stations[i] - 1];
        if (demand < 0) {
            currentLoad -= demand; // demand é negativo, então -= torna positivo
        }
    }
    
    maxLoad = currentLoad; // carga inicial
    
    // simulando a rota pra encontrar a carga máxima - otimizado
    for (size_t i = 1; i < stations.size() - 1; ++i) {
        currentLoad += demandas[stations[i] - 1]; // demanda negativa diminui, positiva aumenta
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
            default:
                break;
        }

        double newCost = calculateTotalCost(newRoutes, instance);

        if (newCost < currentCost - 1e-9) {
            routes = newRoutes;
            currentCost = newCost;
            k = 0; // reinicia a vizinhança com a nova solução melhor
        } else {
            k++; // incrementa a vizinhança
        }
    }

    return routes;
}

vector<Route> performSwap(vector<Route> routes, const Instance *instance) {
    const auto& distancias = instance->get_distancias();
    
    // SWAP INTRA-ROTA - First Improvement com cálculo incremental
    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        if (routes[i].stations.size() <= 3) continue;

        for (size_t j = 1; j < routes[i].stations.size() - 1; j++) {
            for (size_t k = j + 1; k < routes[i].stations.size() - 1; k++) {
                

                double delta = 0;
                int s_j_prev = routes[i].stations[j - 1];
                int s_j = routes[i].stations[j];
                int s_j_next = routes[i].stations[j + 1];
                int s_k_prev = routes[i].stations[k - 1];
                int s_k = routes[i].stations[k];
                int s_k_next = routes[i].stations[k + 1];

                if (k == j + 1) { // Estações adjacentes - caso especial otimizado
                    delta = -distancias[s_j_prev][s_j] - distancias[s_j][s_k] - distancias[s_k][s_k_next]
                            + distancias[s_j_prev][s_k] + distancias[s_k][s_j] + distancias[s_j][s_k_next];
                } else { // Estações não adjacentes
                    delta = -distancias[s_j_prev][s_j] - distancias[s_j][s_j_next] 
                            - distancias[s_k_prev][s_k] - distancias[s_k][s_k_next]
                            + distancias[s_j_prev][s_k] + distancias[s_k][s_j_next]
                            + distancias[s_k_prev][s_j] + distancias[s_j][s_k_next];
                }

                if (delta < -1e-9) {
                    swap(routes[i].stations[j], routes[i].stations[k]);
                    int newCapacity = calculateRouteCapacity(routes[i].stations, instance);
                    if (newCapacity <= instance->get_capVeiculos()) {
                        routes[i].capacity = newCapacity;
                        return routes; // First Improvement - retorna imediatamente
                    } else {
                        swap(routes[i].stations[j], routes[i].stations[k]); // Reverte
                    }
                }
            }
        }
    }

    // SWAP INTER-ROTAS - First Improvement com cálculo incremental
    for (int r1_idx = 0; r1_idx < static_cast<int>(routes.size()); r1_idx++) {
        for (int r2_idx = r1_idx + 1; r2_idx < static_cast<int>(routes.size()); r2_idx++) {
            for (size_t s1_idx = 1; s1_idx < routes[r1_idx].stations.size() - 1; s1_idx++) {
                for (size_t s2_idx = 1; s2_idx < routes[r2_idx].stations.size() - 1; s2_idx++) {
                    
                    int s1_prev = routes[r1_idx].stations[s1_idx - 1];
                    int s1 = routes[r1_idx].stations[s1_idx];
                    int s1_next = routes[r1_idx].stations[s1_idx + 1];
                    
                    int s2_prev = routes[r2_idx].stations[s2_idx - 1];
                    int s2 = routes[r2_idx].stations[s2_idx];
                    int s2_next = routes[r2_idx].stations[s2_idx + 1];

                    double delta = - (distancias[s1_prev][s1] + distancias[s1][s1_next])
                                   - (distancias[s2_prev][s2] + distancias[s2][s2_next])
                                   + (distancias[s1_prev][s2] + distancias[s2][s1_next])
                                   + (distancias[s2_prev][s1] + distancias[s1][s2_next]);

                    if (delta < -1e-9) {
                        swap(routes[r1_idx].stations[s1_idx], routes[r2_idx].stations[s2_idx]);

                        int newCap1 = calculateRouteCapacity(routes[r1_idx].stations, instance);
                        int newCap2 = calculateRouteCapacity(routes[r2_idx].stations, instance);

                        if (newCap1 <= instance->get_capVeiculos() && newCap2 <= instance->get_capVeiculos()) {
                            routes[r1_idx].capacity = newCap1;
                            routes[r2_idx].capacity = newCap2;
                            return routes;
                        } else {
                            swap(routes[r1_idx].stations[s1_idx], routes[r2_idx].stations[s2_idx]); // Reverte
                        }
                    }
                }
            }
        }
    }
    
    return routes; // Nenhuma melhoria encontrada
}

vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance) {
    // First Improvement Strategy
    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        if (routes[i].stations.size() <= 4) continue;

        for (size_t j = 1; j < routes[i].stations.size() - 2; j++) {
            for (size_t k = j + 1; k < routes[i].stations.size() - 1; k++) {
                
                double delta = -instance->get_distancias()[routes[i].stations[j - 1]][routes[i].stations[j]]
                               - instance->get_distancias()[routes[i].stations[k]][routes[i].stations[k + 1]]
                               + instance->get_distancias()[routes[i].stations[j - 1]][routes[i].stations[k]]
                               + instance->get_distancias()[routes[i].stations[j]][routes[i].stations[k + 1]];

                if (delta < -1e-9) {
                    
                    vector<int> newStations = routes[i].stations;
                    reverse(newStations.begin() + j, newStations.begin() + k + 1);
                    
                    int newCapacity = calculateRouteCapacity(newStations, instance);

                    if (newCapacity <= instance->get_capVeiculos()) {
                        routes[i].stations = newStations;
                        routes[i].capacity = newCapacity;
                        return routes;
                    }
                }
            }
        }
    }
    return routes; // Nenhuma melhoria encontrada
}

vector<Route> performReinsertion(vector<Route> routes, const Instance *instance) {
    const auto& distancias = instance->get_distancias();

    // REINSERTION INTRA-ROTA - First Improvement
    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        if (routes[i].stations.size() <= 3) continue;

        for (size_t j = 1; j < routes[i].stations.size() - 1; j++) {
            for (size_t k = 1; k < routes[i].stations.size() - 1; k++) {
                if (j == k || j == k-1) continue;

                int s_j_prev = routes[i].stations[j - 1];
                int s_j = routes[i].stations[j];
                int s_j_next = routes[i].stations[j + 1];
                int s_k_prev = routes[i].stations[k - 1];
                int s_k = routes[i].stations[k];

                double delta = - (distancias[s_j_prev][s_j] + distancias[s_j][s_j_next])
                               - distancias[s_k_prev][s_k]
                               + distancias[s_j_prev][s_j_next]
                               + distancias[s_k_prev][s_j] + distancias[s_j][s_k];

                if (delta < -1e-9) {
                    vector<int> tempStations = routes[i].stations;
                    tempStations.erase(tempStations.begin() + j);
                    tempStations.insert(tempStations.begin() + (k > j ? k-1 : k), s_j);
                    
                    int newCapacity = calculateRouteCapacity(tempStations, instance);
                    if (newCapacity <= instance->get_capVeiculos()) {
                        routes[i].stations = tempStations;
                        routes[i].capacity = newCapacity;
                        return routes;
                    }
                }
            }
        }
    }

    // REINSERTION INTER-ROTAS
    for (int r1_idx = 0; r1_idx < static_cast<int>(routes.size()); r1_idx++) {
        if (routes[r1_idx].stations.size() <= 2) continue;

        for (int r2_idx = 0; r2_idx < static_cast<int>(routes.size()); r2_idx++) {
            if (r1_idx == r2_idx) continue;

            for (size_t s1_idx = 1; s1_idx < routes[r1_idx].stations.size() - 1; s1_idx++) {
                for (size_t s2_idx = 1; s2_idx < routes[r2_idx].stations.size(); s2_idx++) {
                    
                    int s1_prev = routes[r1_idx].stations[s1_idx - 1];
                    int s1 = routes[r1_idx].stations[s1_idx];
                    int s1_next = routes[r1_idx].stations[s1_idx + 1];
                    int s2_prev = routes[r2_idx].stations[s2_idx - 1];
                    int s2 = routes[r2_idx].stations[s2_idx];

                    double delta = - (distancias[s1_prev][s1] + distancias[s1][s1_next])
                                   - distancias[s2_prev][s2]
                                   + distancias[s1_prev][s1_next]
                                   + distancias[s2_prev][s1] + distancias[s1][s2];

                    if (delta < -1e-9) {
                        vector<int> r1_stations = routes[r1_idx].stations;
                        vector<int> r2_stations = routes[r2_idx].stations;
                        
                        r1_stations.erase(r1_stations.begin() + s1_idx);
                        r2_stations.insert(r2_stations.begin() + s2_idx, s1);

                        int newCap1 = (r1_stations.size() > 2) ? calculateRouteCapacity(r1_stations, instance) : 0;
                        int newCap2 = calculateRouteCapacity(r2_stations, instance);

                        if (newCap1 <= instance->get_capVeiculos() && newCap2 <= instance->get_capVeiculos()) {
                            routes[r2_idx] = {r2_stations, newCap2};
                            if (r1_stations.size() <= 2) {
                                routes.erase(routes.begin() + r1_idx);
                            } else {
                                routes[r1_idx] = {r1_stations, newCap1};
                            }
                            return routes;
                        }
                    }
                }
            }
        }
    }
    
    return routes; // Nenhuma melhoria encontrada
}

vector<Route> GRASPConstruction(const Instance *instance) {
    vector<Route> routes;
    vector<int> demandas = instance->get_demandas(); 
    vector<vector<int>> distancias = instance->get_distancias();
    
    // 0 = guloso, 1 = totalmente aleatório
    float alpha = instance->alpha;
    
    // ordem aleatória de criação das rotas iniciais
    vector<int> stationOrder;
    for (int i = 1; i <= instance->get_qtdEstacoes(); i++) {
        stationOrder.push_back(i);
    }
    
    // embaralhar a ordem das estações
    for (int i = 0; i < static_cast<int>(stationOrder.size()); i++) {
        int j = rand() % stationOrder.size();
        swap(stationOrder[i], stationOrder[j]);
    }
    
    // novas rotas na ordem embaralhada
    for (int stationId : stationOrder) {
        int capacityNeeded = abs(demandas[stationId - 1]); // -1 porque demandas não inclui depósito
        routes.push_back({ {0, stationId, 0}, capacityNeeded});
    }

    // fazer até não conseguir mais fazer uniões viáveis
    bool foundImprovement = true;
    while (foundImprovement && routes.size() > 1) {
        foundImprovement = false;
        vector<Saving> candidateSavings;

        // todas as economias viáveis possíveis
        for (int i = 0; i < static_cast<int>(routes.size()); i++) {
            for (int j = i + 1; j < static_cast<int>(routes.size()); j++) {
                // testa se as rotas podem ser unidas
                if (routes[i].capacity + routes[j].capacity <= instance->get_capVeiculos()) {
                    // estações que podem ser conectadas
                    vector<int> connectableStations_i, connectableStations_j;
                    
                    // estações conectáveis da rota i (primeira e última estação não-depósito)
                    if (routes[i].stations.size() > 2) {
                        connectableStations_i.push_back(routes[i].stations[1]); // primeira estação
                        connectableStations_i.push_back(routes[i].stations[routes[i].stations.size()-2]); // última estação
                    }
                    
                    // estações conectáveis da rota j
                    if (routes[j].stations.size() > 2) {
                        connectableStations_j.push_back(routes[j].stations[1]); // primeira estação
                        connectableStations_j.push_back(routes[j].stations[routes[j].stations.size()-2]); // última estação
                    }
                    
                    // savings para todas as combinações possíveis
                    for (int station_i : connectableStations_i) {
                        for (int station_j : connectableStations_j) {
                            double savingValue = distancias[0][station_i] + distancias[0][station_j] - distancias[station_i][station_j];
                            candidateSavings.push_back({savingValue, station_i, station_j});
                        }
                    }
                }
            }
        }

        if (candidateSavings.empty()) break;

        sort(candidateSavings.begin(), candidateSavings.end(), compareSavings);

        // Lista Restrita de Candidatos (RCL)
        double bestSaving = candidateSavings[0].value;
        double worstSaving = candidateSavings.back().value;
        double threshold = worstSaving + alpha * (bestSaving - worstSaving);
        
        vector<Saving> rcl;
        for (const auto &saving : candidateSavings) {
            if (saving.value >= threshold) {
                rcl.push_back(saving);
            }
        }

        // pegando aleatoriamente um saving da RCL
        if (!rcl.empty()) {
            int randomIndex = rand() % rcl.size();
            Saving selectedSaving = rcl[randomIndex];
            
            int route_i_idx = -1, route_j_idx = -1;
            
            for (int k = 0; k < static_cast<int>(routes.size()); k++) {
                // valida se a rota contém fromStation
                for (int station : routes[k].stations) {
                    if (station == selectedSaving.fromStation) {
                        route_i_idx = k;
                        break;
                    }
                }
                // valida se a rota contém toStation
                for (int station : routes[k].stations) {
                    if (station == selectedSaving.toStation) {
                        route_j_idx = k;
                        break;
                    }
                }
            }

            // valida se encontrou ambas as rotas e se são diferentes
            if (route_i_idx != -1 && route_j_idx != -1 && route_i_idx != route_j_idx) {
                // valida restrição de capacidade
                if (routes[route_i_idx].capacity + routes[route_j_idx].capacity <= instance->get_capVeiculos()) {
                    // une das rotas
                    vector<int> newStations = routes[route_i_idx].stations;
                    newStations.pop_back(); // remove o depósito final
                    
                    for (int j = 1; j < static_cast<int>(routes[route_j_idx].stations.size()); j++) {
                        newStations.push_back(routes[route_j_idx].stations[j]);
                    }

                    int newCapacity = routes[route_i_idx].capacity + routes[route_j_idx].capacity;
                    
                    // remove as rotas antigas (maior índice primeiro)
                    routes.erase(routes.begin() + max(route_i_idx, route_j_idx));
                    routes.erase(routes.begin() + min(route_i_idx, route_j_idx));
                    
                    // nova rota
                    routes.push_back({newStations, newCapacity});
                    foundImprovement = true;
                }
            }
        }
    }

    return routes;
}
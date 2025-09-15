#include "Solver.hpp"
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <ctime>

#define QTD_VIZINHANCAS 3
// #define VALIDATE_ROUTES

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
    double totalCost = 0.0;
    for (const auto &route : routes) {
        for (size_t i = 0; i < route.stations.size() - 1; ++i) {
            totalCost += instance->get_distancias()[route.stations[i]][route.stations[i + 1]];
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

bool validateRoutes(const vector<Route> &routes, const Instance *instance) {
    vector<bool> visitedStations(instance->get_qtdEstacoes() + 1, false);
    visitedStations[0] = true; // deposito
    vector<int> demandas = instance->get_demandas();
    
    cout << "\n=== VALIDAÇÃO DAS ROTAS ===\n";

    for (int rota = 0; rota < static_cast<int>(routes.size()); rota++) {
        const Route &route = routes[rota];
        cout << "\nvalidando Rota " << (rota + 1) << ": ";
        printRoute(route);
        
        // verifica se a rota começa e termina no depósito
        if (route.stations.empty() || route.stations[0] != 0 || route.stations.back() != 0) {
            cout << "❌ rota não começa/termina no depósito!\n";
            return false;
        }
        
        // Verificar se há estações suficientes na rota
        if (route.stations.size() < 3) {
            cout << "❌ rota.stations.size < 3: " << route.stations.size() << "\n";
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
        cout << "   carga inicial necessária: " << initialLoad << "\n";
        
        // rota estação por estação (- depósito inicial e final)
        for (int i = 1; i < static_cast<int>(route.stations.size()) - 1; i++) {
            int station = route.stations[i];
            
            // validando estação
            if (station == 0) {
                cout << "❌ depósito aparece no meio da rota\n";
                return false;
            }
            
            // validando unica visita na estação
            if (visitedStations[station]) {
                cout << "❌ estação " << station << " visitada dnv\n";
                return false;
            }
            visitedStations[station] = true;
            
            // aplicar demanda da estação
            currentLoad += demandas[station - 1]; // -1 porque demandas não inclui depósito
            
            // Verificar limites de capacidade
            if (currentLoad < 0) {
                cout << "❌ capacidade negativa (" << currentLoad << ") após visitar estação " << station << "!\n";
                cout << "   demanda da estação: " << demandas[station - 1] << "\n";
                return false;
            }
            
            if (currentLoad > instance->get_capVeiculos()) {
                cout << "❌ capacidade excedida (" << currentLoad << "/" << instance->get_capVeiculos() 
                     << ") após visitar estação " << station << "!\n";
                return false;
            }
            
            string action = (demandas[station - 1] < 0) ? "retirou" : "entregou";
            cout << "   Estação " << station << " (" << action << " " << abs(demandas[station - 1]) 
                 << " bicicletas) -> Carga atual: " << currentLoad << "/" << instance->get_capVeiculos() << "\n";
        }

        cout << "✅ Rota " << (rota + 1) << " válida! Carga final: " << currentLoad 
             << " (inicial: " << initialLoad << ")\n";
    }
    
    // todas estações visitadas
    for (int i = 1; i <= instance->get_qtdEstacoes(); i++) {
        if (!visitedStations[i]) {
            cout << "❌ estação " << i << " não foi visitada!\n";
            return false;
        }
    }
    
    cout << "\n✅ TODAS AS ROTAS SÃO VÁLIDAS!\n";
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
    vector<int> demandas = instance->get_demandas();
    int maxLoad = 0;
    int currentLoad = 0;
    
    // 1° calcular carga inicial necessária (somatório de demandas negativas)
    for (int i = 1; i < static_cast<int>(stations.size()) - 1; i++) {
        int station = stations[i];
        if (demandas[station - 1] < 0) {
            currentLoad += abs(demandas[station - 1]);
        }
    }
    
    maxLoad = currentLoad; // carga inicial
    
    // simulando a rota pra encontrar a carga máxima
    for (int i = 1; i < static_cast<int>(stations.size()) - 1; i++) {
        int station = stations[i];
        currentLoad += demandas[station - 1]; // demanda negativa diminui, positiva aumenta
        maxLoad = max(maxLoad, currentLoad);
    }
    
    return maxLoad;
}

vector<Route> VND(vector<Route> routes, const Instance *instance) {
    vector<Vizinhança> neighborhoodOrder = {SWAP, TWO_OPT, REINSERTION};
    
    int k = 0; // vizinhança atual

    while (k < QTD_VIZINHANCAS) {
        vector<Route> newRoutes;
        double currentCost = calculateTotalCost(routes, instance);

        switch (neighborhoodOrder[k]) {
            case SWAP:
                newRoutes = performSwap(routes, instance);
                break;
            case TWO_OPT:
                newRoutes = performTwoOpt(routes, instance);
                break;
            case REINSERTION:
                newRoutes = performReinsertion(routes, instance);
                break;
            default:
                break;
        }

        double newCost = calculateTotalCost(newRoutes, instance);

        if (newCost < currentCost) {
            routes = newRoutes;
            k = 0; // reinicia a vizinhança com a nova solução melhor
        } else {
            k++; // incrementa a vizinhança
        }
    }

    return routes;
}

vector<Route> performSwap(vector<Route> routes, const Instance *instance) {
    vector<Route> bestRoutes = routes;
    double bestCost = calculateTotalCost(bestRoutes, instance);
    vector<int> demandas = instance->get_demandas();

    // troca estações dentro da mesma rota - SWAP INTRA-ROTA
    for (int i = 0; i < static_cast<int>(routes.size()); i++) {
        vector<int> &currentStations = routes[i].stations;

        if (currentStations.size() <= 3) continue; // não há estações para trocar

        // percorre todas combinações ignorando deposito (0 e n-1)
        for (int j = 1; j < static_cast<int>(currentStations.size()) - 1; j++) {
            for (int k = j + 1; k < static_cast<int>(currentStations.size()) - 1; k++) {
                vector<Route> tempRoutes = routes;
                swap(tempRoutes[i].stations[j], tempRoutes[i].stations[k]);

                double tempCost = calculateTotalCost(tempRoutes, instance);

                if (tempCost < bestCost) {
                    bestCost = tempCost;
                    bestRoutes = tempRoutes;
                }
            }
        }
    }

    // trocar estações entre rotas - SWAP INTER-ROTA
    for (int route1 = 0; route1 < static_cast<int>(routes.size()); route1++) {
        for (int route2 = route1 + 1; route2 < static_cast<int>(routes.size()); route2++) {
            
            if (routes[route1].stations.size() <= 2 || routes[route2].stations.size() <= 2) continue;
            
            // tenta trocar cada estação da rota1 com cada estação da rota2
            for (int pos1 = 1; pos1 < static_cast<int>(routes[route1].stations.size()) - 1; pos1++) {
                for (int pos2 = 1; pos2 < static_cast<int>(routes[route2].stations.size()) - 1; pos2++) {
                    
                    int station1 = routes[route1].stations[pos1];
                    int station2 = routes[route2].stations[pos2];
                    
                    // calcula mudança de capacidade
                    int demand1 = abs(demandas[station1 - 1]);
                    int demand2 = abs(demandas[station2 - 1]);
                    
                    int newCapacity1 = routes[route1].capacity - demand1 + demand2;
                    int newCapacity2 = routes[route2].capacity - demand2 + demand1;
                    
                    // ve a viabilidade da troca
                    if (newCapacity1 <= instance->get_capVeiculos() && 
                        newCapacity2 <= instance->get_capVeiculos()) {
                        
                        vector<Route> tempRoutes = routes;
                        
                        // troca
                        tempRoutes[route1].stations[pos1] = station2;
                        tempRoutes[route2].stations[pos2] = station1;
                        tempRoutes[route1].capacity = newCapacity1;
                        tempRoutes[route2].capacity = newCapacity2;
                        
                        double tempCost = calculateTotalCost(tempRoutes, instance);
                        
                        if (tempCost < bestCost) {
                            bestCost = tempCost;
                            bestRoutes = tempRoutes;
                        }
                    }
                }
            }
        }
    }

    return bestRoutes;
}

vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance) {
    vector<Route> bestRoutes = routes;
    double bestCost = calculateTotalCost(bestRoutes, instance);
 
    // melhorar ordem dentro de cada rota - 2-OPT INTRA-ROTA
    for (int route = 0; route < static_cast<int>(routes.size()); route++) {
        vector<int> &currentStations = routes[route].stations;

        if (currentStations.size() <= 4) continue; // tem q ter no minimo 4 estações pro 2-opt

        for (int i = 1; i < static_cast<int>(currentStations.size()) - 2; i++) {
            for (int j = i + 1; j < static_cast<int>(currentStations.size()) - 1; j++) {
                //  delta antes da modificação
                double delta = -instance->get_distancias()[currentStations[i - 1]][currentStations[i]]
                               - instance->get_distancias()[currentStations[j]][currentStations[j + 1]]
                               + instance->get_distancias()[currentStations[i - 1]][currentStations[j]]
                               + instance->get_distancias()[currentStations[i]][currentStations[j + 1]];

                if (delta < 0) { // só fazer a troca se melhorar
                    vector<int> tempStations = currentStations;
                    reverse(tempStations.begin() + i, tempStations.begin() + j + 1);

                    vector<Route> tempRoutes = routes;
                    tempRoutes[route].stations = tempStations;
                    
                    // capacidade real após a mudança
                    tempRoutes[route].capacity = calculateRouteCapacity(tempStations, instance);
                    
                    // valida capacidade
                    if (tempRoutes[route].capacity <= instance->get_capVeiculos()) {
                        double tempCost = calculateTotalCost(tempRoutes, instance);

                        if (tempCost < bestCost) {
                            bestCost = tempCost;
                            bestRoutes = tempRoutes;
                        }
                    }
                }
            }
        }
    }

    return bestRoutes;
}

vector<Route> performReinsertion(vector<Route> routes, const Instance *instance) {
    vector<Route> bestRoutes = routes;
    double bestCost = calculateTotalCost(bestRoutes, instance);
    vector<int> demandas = instance->get_demandas();

    // mover estação dentro da mesma rota - REINSERTION INTRA-ROTA
    for (int route = 0; route < static_cast<int>(routes.size()); route++) {
        vector<int> &currentStations = routes[route].stations;

        if (currentStations.size() <= 3) continue; // sem estações suficientes pra realocar

        for (int i = 1; i < static_cast<int>(currentStations.size()) - 1; i++) {
            int stationToMove = currentStations[i];
            
            for (int j = 1; j < static_cast<int>(currentStations.size()) - 1; j++) {
                if (i == j) continue;

                vector<int> tempStations = currentStations;
                tempStations.erase(tempStations.begin() + i); // removendo a estação

                if (j < i)  tempStations.insert(tempStations.begin() + j, stationToMove); // inserindo na nova posição
                else        tempStations.insert(tempStations.begin() + j - 1, stationToMove); // ajusta o índice se a remoção foi antes da inserção

                // capacidade não muda em reinsertion (mesmas estações, só ordem diferente)
                vector<Route> tempRoutes = routes;
                tempRoutes[route].stations = tempStations;

                double tempCost = calculateTotalCost(tempRoutes, instance);

                if (tempCost < bestCost) {
                    bestCost = tempCost;
                    bestRoutes = tempRoutes;
                }
            }
        }
    }

    // mover estação de uma rota para outra - REINSERTION INTER-ROTA
    for (int fromRoute = 0; fromRoute < static_cast<int>(routes.size()); fromRoute++) {
        for (int toRoute = 0; toRoute < static_cast<int>(routes.size()); toRoute++) {
            if (fromRoute == toRoute) continue;

            // tenta mover cada estação da rota origem (mesmo se a rota tiver apenas uma estação)
            for (int stationPos = 1; stationPos < static_cast<int>(routes[fromRoute].stations.size()) - 1; stationPos++) {
                int stationToMove = routes[fromRoute].stations[stationPos];

                // tenta inserir em diferentes posições da rota destino
                for (int insertPos = 1; insertPos < static_cast<int>(routes[toRoute].stations.size()); insertPos++) {
                    vector<Route> tempRoutes = routes;

                    // remove estação da rota origem
                    tempRoutes[fromRoute].stations.erase(tempRoutes[fromRoute].stations.begin() + stationPos);

                    // insere estação na rota destino
                    tempRoutes[toRoute].stations.insert(tempRoutes[toRoute].stations.begin() + insertPos, stationToMove);

                    // valida capacidade
                    int newToCap = calculateRouteCapacity(tempRoutes[toRoute].stations, instance);
                    if (newToCap > instance->get_capVeiculos()) {
                        continue; // destino inválido
                    }

                    // remove rota se ficar vazia
                    bool removedFromRoute = false;
                    if (static_cast<int>(tempRoutes[fromRoute].stations.size()) == 2) { // [0, 0]
                        // Se removemos a rota com índice menor que o destino, o índice do destino diminui
                        tempRoutes.erase(tempRoutes.begin() + fromRoute);
                    
                        removedFromRoute = true;
                    } else {
                        // nova capacidade da rota origem
                        int newFromCap = calculateRouteCapacity(tempRoutes[fromRoute].stations, instance);
                        if (newFromCap > instance->get_capVeiculos()) {
                            continue; // origem inválida
                        }
                        tempRoutes[fromRoute].capacity = newFromCap;
                    }

                    // atualizar capacidade da rota destino
                    if (!removedFromRoute) {
                        tempRoutes[toRoute].capacity = newToCap;
                    } else {
                        // Se removido e fromRoute < toRoute, o índice do destino reduziu 1
                        int adjustedToIndex = toRoute - (fromRoute < toRoute ? 1 : 0);
                        if (adjustedToIndex >= 0 && adjustedToIndex < static_cast<int>(tempRoutes.size())) {
                            tempRoutes[adjustedToIndex].capacity = newToCap;
                        }
                    }

                    double tempCost = calculateTotalCost(tempRoutes, instance);

                    if (tempCost < bestCost) {
                        bestCost = tempCost;
                        bestRoutes = tempRoutes;
                    }
                }
            }
        }
    }

    return bestRoutes;
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
        double threshold = bestSaving - alpha * (bestSaving - worstSaving);
        
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
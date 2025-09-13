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

bool validateRoutes(const vector<Route> &routes, const Instance *instance) {
    vector<bool> visitedStations(instance->get_qtdEstacoes(), false);
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
            if (demandas[station] < 0) {
                initialLoad += abs(demandas[station]); // demandas negativas
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
            currentLoad += demandas[station];
            
            // Verificar limites de capacidade
            if (currentLoad < 0) {
                cout << "❌ capacidade negativa (" << currentLoad << ") após visitar estação " << station << "!\n";
                cout << "   demanda da estação: " << demandas[station] << "\n";
                return false;
            }
            
            if (currentLoad > instance->get_capVeiculos()) {
                cout << "❌ capacidade excedida (" << currentLoad << "/" << instance->get_capVeiculos() 
                     << ") após visitar estação " << station << "!\n";
                return false;
            }
            
            string action = (demandas[station] < 0) ? "retirou" : "entregou";
            cout << "   Estação " << station << " (" << action << " " << abs(demandas[station]) 
                 << " bicicletas) -> Carga atual: " << currentLoad << "/" << instance->get_capVeiculos() << "\n";
        }

        cout << "✅ Rota " << (rota + 1) << " válida! Carga final: " << currentLoad 
             << " (inicial: " << initialLoad << ")\n";
    }
    
    // todas estações visitadas
    for (int i = 1; i < instance->get_qtdEstacoes(); i++) {
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

        // fase 1
        vector<Route> routes = GRASPConstruction(instance);

        // fase 2
        vector<Route> improvedRoutes = VND(routes, instance);
        double improvedCost = calculateTotalCost(improvedRoutes, instance);

        if (improvedCost < bestCost) {
            bestCost = improvedCost;
            bestSolutions = improvedRoutes;
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

// algoritmo de clarke-wright
vector<Route> clarkeWright(const Instance *instance) {
    vector<Route> routes;
    vector<Saving> savings;

    vector<int> demandas = instance->get_demandas(); 
    vector<vector<int>> distancias = instance->get_distancias();
    
    // cada estação tem sua propria rota (exceto o depósito)
    for (int i = 1; i < instance->get_qtdEstacoes(); i++) { // começa em 1 para ignorar o deposito
        int capacityNeeded = abs(demandas[i]); // Capacidade necessária é o valor absoluto da demanda
        routes.push_back({ {0, i, 0}, capacityNeeded}); // deposito -> estacao -> deposito
    }

    // calcular economias: S(i,j) = d(0,i) + d(0,j) - d(i,j)
    for (int i = 1; i < instance->get_qtdEstacoes(); i++) {
        for (int j = i + 1; j < instance->get_qtdEstacoes(); j++) {
            double savingValue = distancias[0][i] + distancias[0][j] - distancias[i][j];
            savings.push_back({ savingValue, i, j});
        }
    }

    // ordenar economias
    sort(savings.begin(), savings.end(), compareSavings);

    // processamento guloso: tenta unir as rotas com a maior economia
    for (const auto &saving : savings) {
        // encontrar as rotas que contêm fromStation e toStation
        // onde uma estação está no final da rota e a outra no início
        int route_i_idx = -1; 
        int route_j_idx = -1;

        for ( int k = 0; k < static_cast<int>(routes.size()); k++) {
            if (routes[k].stations.back() == saving.fromStation && routes[k].stations.size() > 2) {
                route_i_idx = k;
            } else if (routes[k].stations.back() == saving.toStation && routes[k].stations.size() > 2) {
                route_j_idx = k;
            } else if (routes[k].stations[1] == saving.fromStation) {
                route_i_idx = k;
            } else if (routes[k].stations[1] == saving.toStation) {
                route_j_idx = k;
            }
        }

        if (route_i_idx == route_j_idx || route_i_idx == -1 || route_j_idx == -1) {
            continue; // mesma rota ou rota não encontrada
        }

        if (routes[route_i_idx].capacity + routes[route_j_idx].capacity <= instance->get_capVeiculos()) {
            // faz a união das rotas
            vector<int> newStations = routes[route_i_idx].stations;
            newStations.pop_back(); // remove o deposito final
            for ( int j = 1; j < static_cast<int>(routes[route_j_idx].stations.size()); j++) {
                newStations.push_back(routes[route_j_idx].stations[j]); // adiciona todas as estações da segunda rota, incluindo o depósito final
            }

            int newCapacity = routes[route_i_idx].capacity + routes[route_j_idx].capacity;
            routes.erase(routes.begin() + max(route_i_idx, route_j_idx)); // remove a rota com o índice maior primeiro
            routes.erase(routes.begin() + min(route_i_idx, route_j_idx)); // depois remove a outra
            routes.push_back({ newStations, newCapacity }); 
        }
    }

    return routes;
}

vector<Route> VND(vector<Route> routes, const Instance *instance) {
    int k = 1;

    // se encontrar uma solução melhor, atualiza routes e reinicia a busca com mesmo k
    // caso contrário, incrementar k
    while (k <= QTD_VIZINHANCAS) {
        vector<Route> newRoutes;
        double currentCost = calculateTotalCost(routes, instance);

        switch (k) {
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
            k = 1; // reinicia a vizinhança com a nova solução melhor
        } else k++; // incrementa a vizinhança
    }

    return routes;
}

vector<Route> performSwap(vector<Route> routes, const Instance *instance) {
    vector<Route> bestRoutes = routes;
    double bestCost = calculateTotalCost(bestRoutes, instance);

    for (int i = 0; i < static_cast<int>(routes.size()); i++) { // passando por todas as rotas
        vector<int> &currentStations = routes[i].stations;

        if (currentStations.size() <= 3) continue; // não há estações para trocar

        for (int j = 1; j < static_cast<int>(currentStations.size()) -1; j++) {
            for (int k = j + 1; k < static_cast<int>(currentStations.size()) -1; k++) {
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

    return bestRoutes;
}

vector<Route> performTwoOpt(vector<Route> routes, const Instance *instance) {
    vector<Route> bestRoutes = routes;
    double bestCost = calculateTotalCost(bestRoutes, instance);
 
    for (int route = 0; route < static_cast<int>(routes.size()); route++) { // passando por todas as rotas
        vector<int> &currentStations = routes[route].stations;

        if (currentStations.size() <= 4) continue; // minimo de 4 estações para fazer o 2-opt

        for (int i = 1; i < static_cast<int>(currentStations.size()) -2; i++) {
            for (int j = i + 1; j < static_cast<int>(currentStations.size()) -1; j++) {
                double delta = -instance->get_distancias()[currentStations[i - 1]][currentStations[i]]
                               - instance->get_distancias()[currentStations[j]][currentStations[j + 1]]
                               + instance->get_distancias()[currentStations[i - 1]][currentStations[j]]
                               + instance->get_distancias()[currentStations[i]][currentStations[j + 1]];

                if (delta < 0) { // encontrou melhoria
                    vector<int> tempStations = currentStations;
                    reverse(tempStations.begin() + i, tempStations.begin() + j + 1);

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
    }

    return bestRoutes;
}

vector<Route> performReinsertion(vector<Route> routes, const Instance *instance) {
    vector<Route> bestRoutes = routes;
    double bestCost = calculateTotalCost(bestRoutes, instance);

    for (int route = 0; route < static_cast<int>(routes.size()); route++) {
        vector<int> &currentStations = routes[route].stations;

        if (currentStations.size() <= 3 ) continue; // sem estações suficientes pra realocar

        for (int i = 1; i < static_cast<int>(currentStations.size()) - 1; i++) {
            int stationToMove = currentStations[i];
            
            for (int j = 1; j < static_cast<int>(currentStations.size()) - 1; j++) {
                if (i == j) continue;

                vector<int> tempStations = currentStations;
                tempStations.erase(tempStations.begin() + i); // removendo a estação

                if ( j < i ) tempStations.insert(tempStations.begin() + j, stationToMove); // inserindo na nova posição
                else tempStations.insert(tempStations.begin() + j - 1, stationToMove); // ajusta o índice se a remoção foi antes da inserção

                // capacidade não muda em reinsertion (mesmas estações, só ordem diferente)
                // tempRoutes[route].capacity permanece a mesma
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

    return bestRoutes;
}

vector<Route> GRASPConstruction(const Instance *instance) {
    vector<Route> routes;
    vector<int> demandas = instance->get_demandas(); 
    vector<vector<int>> distancias = instance->get_distancias();
    
    // Parâmetro alpha para controlar a aleatoriedade (0 = guloso, 1 = totalmente aleatório)
    float alpha = instance->alpha;
    
    // Cada estação tem sua própria rota inicial (exceto o depósito)
    for (int i = 1; i < instance->get_qtdEstacoes(); i++) {
        int capacityNeeded = abs(demandas[i]);
        routes.push_back({ {0, i, 0}, capacityNeeded});
    }

    // Continua até não conseguir mais fazer uniões viáveis
    bool foundImprovement = true;
    while (foundImprovement && routes.size() > 1) {
        foundImprovement = false;
        vector<Saving> candidateSavings;

        // Calcular todas as economias viáveis possíveis
        for (int i = 0; i < static_cast<int>(routes.size()); i++) {
            for (int j = i + 1; j < static_cast<int>(routes.size()); j++) {
                // Verificar se as rotas podem ser unidas
                if (routes[i].capacity + routes[j].capacity <= instance->get_capVeiculos()) {
                    // Encontrar as estações que podem ser conectadas
                    vector<int> connectableStations_i, connectableStations_j;
                    
                    // Estações conectáveis da rota i (primeira e última estação não-depósito)
                    if (routes[i].stations.size() > 2) {
                        connectableStations_i.push_back(routes[i].stations[1]); // primeira estação
                        connectableStations_i.push_back(routes[i].stations[routes[i].stations.size()-2]); // última estação
                    }
                    
                    // Estações conectáveis da rota j
                    if (routes[j].stations.size() > 2) {
                        connectableStations_j.push_back(routes[j].stations[1]); // primeira estação
                        connectableStations_j.push_back(routes[j].stations[routes[j].stations.size()-2]); // última estação
                    }
                    
                    // Calcular savings para todas as combinações possíveis
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

        // Ordenar savings em ordem decrescente
        sort(candidateSavings.begin(), candidateSavings.end(), compareSavings);

        // Criar Lista Restrita de Candidatos (RCL)
        double bestSaving = candidateSavings[0].value;
        double worstSaving = candidateSavings.back().value;
        double threshold = bestSaving - alpha * (bestSaving - worstSaving);
        
        vector<Saving> rcl;
        for (const auto &saving : candidateSavings) {
            if (saving.value >= threshold) {
                rcl.push_back(saving);
            }
        }

        // Escolher aleatoriamente um saving da RCL
        if (!rcl.empty()) {
            int randomIndex = rand() % rcl.size();
            Saving selectedSaving = rcl[randomIndex];
            
            // Encontrar as rotas que contêm as estações do saving selecionado
            int route_i_idx = -1, route_j_idx = -1;
            
            for (int k = 0; k < static_cast<int>(routes.size()); k++) {
                // Verificar se a rota contém fromStation
                for (int station : routes[k].stations) {
                    if (station == selectedSaving.fromStation) {
                        route_i_idx = k;
                        break;
                    }
                }
                // Verificar se a rota contém toStation
                for (int station : routes[k].stations) {
                    if (station == selectedSaving.toStation) {
                        route_j_idx = k;
                        break;
                    }
                }
            }

            // Verificar se encontrou ambas as rotas e se são diferentes
            if (route_i_idx != -1 && route_j_idx != -1 && route_i_idx != route_j_idx) {
                // Verificar restrição de capacidade
                if (routes[route_i_idx].capacity + routes[route_j_idx].capacity <= instance->get_capVeiculos()) {
                    // Fazer a união das rotas
                    vector<int> newStations = routes[route_i_idx].stations;
                    newStations.pop_back(); // remove o depósito final
                    
                    for (int j = 1; j < static_cast<int>(routes[route_j_idx].stations.size()); j++) {
                        newStations.push_back(routes[route_j_idx].stations[j]);
                    }

                    int newCapacity = routes[route_i_idx].capacity + routes[route_j_idx].capacity;
                    
                    // Remover as rotas antigas (maior índice primeiro)
                    routes.erase(routes.begin() + max(route_i_idx, route_j_idx));
                    routes.erase(routes.begin() + min(route_i_idx, route_j_idx));
                    
                    // Adicionar a nova rota
                    routes.push_back({newStations, newCapacity});
                    foundImprovement = true;
                }
            }
        }
    }

    return routes;
}
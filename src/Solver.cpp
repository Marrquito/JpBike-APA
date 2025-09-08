#include "Solver.hpp"
#include <iomanip>
#include <algorithm>

struct Saving {
    double value;
    int fromStation;
    int toStation;
};

void Solver::Solve(Instance *instance) {
    // primeiro fazer o algoritmo guloso
    vector<Route>routes = clarkeWright(instance);
    double valorOtimo = instance->valorOtimo;

    // custo total das rotas
    double custoHeuristica = 0.0;
    for (const auto& route : routes) {
        for (size_t i = 0; i < route.stations.size() - 1; ++i) {
            custoHeuristica += instance->get_distancias()[route.stations[i]][route.stations[i+1]];
        }
    }

    // tabela de comparação
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "--------------------------------------------------------------------------\n";
    std::cout << "| Instancia | Heuristica Construtiva | Valor Otimo |   Gap (%)   |\n";
    std::cout << "--------------------------------------------------------------------------\n";
    std::cout << "|"
              << std::setw(10) << instance->instanceName << " |"
              << std::setw(22) << custoHeuristica << " |"
              << std::setw(12) << valorOtimo << " |"
              << std::setw(11) << ((custoHeuristica - valorOtimo) / valorOtimo * 100.0) << " |\n";
    std::cout << "--------------------------------------------------------------------------\n";

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
    
    // cada estação tem sua propria rota
    for (int i = 0; i < instance->get_qtdEstacoes(); i++) {
        routes.push_back({ {0, i, 0}, demandas[i]}); // deposito -> estacao -> deposito
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
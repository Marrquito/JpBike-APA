#include "Instance.hpp"
#include <iomanip>

//  operador de sobrecarga >>
istream& operator>>(istream& is, Instance& instance) {
    is >> instance.qtdEstacoes >> instance.qtdVeiculos >> instance.capVeiculos;

    if (is.fail()) {
        return is;
    }

    // demandas
    instance.demandas.resize(instance.qtdEstacoes);
    for (int i = 0; i < instance.qtdEstacoes; ++i) {
        is >> instance.demandas[i];
    }
    if (is.fail()) return is;

    // matriz de custos
    instance.distancias.resize(instance.qtdEstacoes, vector<int>(instance.qtdEstacoes));
    for (int i = 0; i < instance.qtdEstacoes; ++i) {
        for (int j = 0; j < instance.qtdEstacoes; ++j) {
            is >> instance.distancias[i][j];
        }
    }
    return is;
}

void Instance::print() const {
    cout << "Quantidade de Estações: " << this->qtdEstacoes << endl;
    cout << "Quantidade de Veículos: " << this->qtdVeiculos << endl;
    cout << "Capacidade dos Veículos: " << this->capVeiculos << endl;

    cout << "Demandas: ";
    for (int i = 0; i < this->qtdEstacoes; ++i) {
        cout << this->demandas[i] << " ";
    }
    cout << endl;

    cout << "Distâncias:" << endl;
    for (int i = 0; i < this->qtdEstacoes; ++i) {
        for (int j = 0; j < this->qtdEstacoes; ++j) {
            cout << setw(4) << this->distancias[i][j] << " ";
        }
        cout << endl;
    }

}

int Instance::get_qtdEstacoes() const {
    return this->qtdEstacoes;
}

vector<int> Instance::get_demandas() const {
    return this->demandas;
}

vector<vector<int>> Instance::get_distancias() const {
    return this->distancias;
}

int Instance::get_capVeiculos() const {
    return this->capVeiculos;
}
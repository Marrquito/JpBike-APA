#ifndef INSTANCE_H
#define INSTANCE_H

#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

class Instance {
    private:
        int qtdEstacoes;
        int qtdVeiculos;
        int capVeiculos;
        std::vector<int> demandas;
        std::vector<std::vector<int>> distancias;

    public:
        Instance() : qtdEstacoes(0), qtdVeiculos(0), capVeiculos(0) {}

        // Sobrecarga do operador de entrada (>>), declarada como friend
        friend istream &operator>>(istream &is, Instance &instance);

        int get_qtdEstacoes() const;
        void print() const;
};

#endif
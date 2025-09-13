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
        vector<int> demandas;
        vector<vector<int>> distancias;


    public:
        int valorOtimo;
        string instanceName;
        float alpha;
        int repetitions;

        Instance(int valor_otimo, string instance_name, float alpha, int repetitions) : qtdEstacoes(0), qtdVeiculos(0), capVeiculos(0), valorOtimo(valor_otimo), instanceName(instance_name), alpha(alpha), repetitions(repetitions) {}

        // Sobrecarga do operador de entrada (>>), declarada como friend
        friend istream &operator>>(istream &is, Instance &instance);

        int get_qtdEstacoes() const;
        int get_capVeiculos() const;
        vector<int> get_demandas() const;
        vector<vector<int>> get_distancias() const;

        void print() const;
};

#endif
#include <iostream>
#include "Instance.hpp"
#include "Solver.hpp"
#include <fstream>

using namespace std;

int main() 
{
    string filename = "";
    cout << "Digite o nome do arquivo de instância: ";
    cin >> filename;

    ifstream file("instancias/" + filename);

    int valor_otimo = 0;
    cout << "Digite o valor ótimo da instância: ";
    cin >> valor_otimo;

    float alpha = 0.0;
    cout << "Digite o valor de alpha (0.0 a 1.0): ";
    cin >> alpha;

    int repetitions = 0;
    cout << "Digite o número de repetições para o GRASP: ";
    cin >> repetitions;
    
    if (!file) 
    {
        cerr << "Erro ao abrir arquivo!" << endl;
        return -1;
    }

    Instance instancia(valor_otimo, filename, alpha, repetitions);

    file >> instancia;
    // instancia.print();

    Solver *solver = new Solver();
    solver->Solve(&instancia);


    delete solver;
    return 0;
}


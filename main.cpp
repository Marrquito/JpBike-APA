#include <iostream>
#include "Instance.hpp"
#include "Solver.hpp"
#include <fstream>

using namespace std;

int main() 
{
    cout << "Hello, karai!" << endl;
    
    string filename;
    cout << "Digite o nome do arquivo de instância: ";
    cin >> filename;

    ifstream file("instancias/" + filename);

    int valor_otimo;
    cout << "Digite o valor ótimo da instância: ";
    cin >> valor_otimo;
    
    if (!file) 
    {
        cerr << "Erro ao abrir arquivo!" << endl;
        return -1;
    }

    Instance instancia(valor_otimo, filename);

    file >> instancia;
    // instancia.print();

    Solver solver;
    solver.Solve(&instancia);

    return 0;
}


#include <iostream>
#include "Instance.hpp"
#include <fstream>

using namespace std;

int main() 
{
    cout << "Hello, karai!" << endl;

    ifstream file("./instancias/n12_q20.txt");
    
    if (!file) 
    {
        cerr << "Erro ao abrir arquivo!" << endl;
    }

    Instance instancia;


    file >> instancia;
    instancia.print();

    return 0;
}


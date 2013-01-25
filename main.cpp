#include <iostream>
#include <fstream>
#include "crystal.h"
#include <armadillo>
#include <ctime>

using namespace arma;
using namespace std;

void printing(int nc, double b){
    ofstream test;
    test.open("/home/jonathan/projectsFSAP/project1/project1/locationatoms.xyz");

    Crystal crystal(nc, b);
    test << crystal << endl;
    test.close();

}

int main()
{
    //cubic lattice with nc x nc x nc cells
    int nc=8;
    //latice parameter in unit Angstrom
    double b=5.28;
    cout << "this is a test!" << endl;
    cout << "Hello World!" << endl;

    printing(nc, b);
    return 0;
}



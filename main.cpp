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
    //UNIT SYSTEM!!
    //Distances = Angstrom
    //Time = picosecond
    //Energy = electronvolt
    //Mass = enter in amu, wrong unit system is addressed for in value of k
    //Temperature = Kelvin
    //velocity = angstrom/ps
    //k = 0.8314766196505026 when masses are expressed in amu and T in K to get velocities in A/ps

    //cubic lattice with nc x nc x nc cells
    int nc=8;
    //latice parameter in unit Angstrom
    double b=5.28;
    cout << "this is a test!" << endl;
    cout << "Hello World!" << endl;

    printing(nc, b);
    return 0;
}



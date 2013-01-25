#include <iostream>
#include <fstream>
#include "crystal.h"
#include <armadillo>
#include <ctime>

using namespace arma;
using namespace std;

void printing(Crystal &crystal){
    ofstream output;
    output.open("/home/jonathan/projectsFSAP/project1/project1/output/locationatoms.xyz");
    output << crystal << endl;
    output.close();
}

void printvelocities(Crystal &crystal){
    ofstream outputx, outputy, outputz, outputtot;
    outputx.open("/home/jonathan/projectsFSAP/project1/project1/output/velocity-x.dat");
    outputy.open("/home/jonathan/projectsFSAP/project1/project1/output/velocity-y.dat");
    outputz.open("/home/jonathan/projectsFSAP/project1/project1/output/velocity-z.dat");
    outputtot.open("/home/jonathan/projectsFSAP/project1/project1/output/velocity-vtot.dat");

    int nc = crystal.nc;
    for(int i=0; i< nc ;++i){
        for(int j=0; j< nc; ++j){
            for(int k=0; k<nc ;++k){
                //now at the level of cells
                for(int l=0; l<4; l++){
                    //now at the level of atoms
                    double vtot=0;
                    for(int m=0; m<3; m++){
                        //initializing velocities of the individual atoms
                        vtot+=(crystal.allcells[i][j][k].atoms[l].phasevect(m+3))*(crystal.allcells[i][j][k].atoms[l].phasevect(m+3));
                    }
                    vtot=sqrt(vtot);
                    outputtot << vtot << endl;
                    outputx << crystal.allcells[i][j][k].atoms[l].phasevect(3) << endl;
                    outputy << crystal.allcells[i][j][k].atoms[l].phasevect(4) << endl;
                    outputz << crystal.allcells[i][j][k].atoms[l].phasevect(5) << endl;
                }
            }
        }
    }
    outputx.close();
    outputy.close();
    outputz.close();
    outputtot.close();
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

    Crystal crystal(nc, b);
    printing(crystal);
    printvelocities(crystal);
    return 0;
}



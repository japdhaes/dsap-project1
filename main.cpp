#include <iostream>
#include <fstream>
#include "crystal.h"
#include <armadillo>
#include <ctime>
#include <vector>
#include "verletalgo.h"
#include <sstream>

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

    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        vec3 velocity= crystal.allatoms[i]->getVelocity();
        outputtot << norm(crystal.allatoms[i]->getVelocity(),3) << endl;
        outputx << velocity(0) << endl;
        outputy << velocity(1) << endl;
        outputz << velocity(2)<< endl;
    }
    outputx.close();
    outputy.close();
    outputz.close();
    outputtot.close();
}

string createname( int number){
    stringstream oss;
    oss<<"/home/jonathan/projectsFSAP/project1/project1/output/";
    oss<<"locationatoms.";
    if(number<10){
        oss<<"000"<<number;
    }
    else if(number <100){
        oss<<"00"<<number;
    }
    else if(number <1000){
        oss<<"0"<<number;
    }

    return oss.str();
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

    int seed = time(0);

    cout << seed<< endl;
    int nc=8;
    //latice parameter in unit Angstrom
    double b=5.28;


    cout << createname(100) <<endl;
    Crystal crystal(nc, b, seed);
    printing(crystal);
    VerletAlgo integrator(crystal);
    for(int j=1; j<2000; j++){
        ofstream output;
        output.open(createname(j).c_str());
        for(unsigned int i=0; i<crystal.allatoms.size(); i++){
            integrator.integrateAtom(crystal.allatoms[i]);
        }
        output << crystal << endl;
    }

    printvelocities(crystal);
    vector<Atom*> atoms;
    return 0;
}
/* vec3 r
 *r<<x<<y<<z;
 *
 *mat h = zeros(2,2); //2x2 matrix
 *same as mat h= zeros<mat>(2,2);
 *h << 1 << 2 << endr
 * << 3 << 4 ;
 *
 *h(0,1) = 1e rij 2e kolom element
 *
 *vec r; r<<1<<2<<3;
 *vec p; p<<0<<1<<2;
 *
 *dot(r,p) <- dotproduct
 *
 *r.t()*p <-- same result!
 *
 *norm(r,2) <= norm
 *r.n_elem <= number of elements in the vector
 *
 *we can do matrix A * vector p
 *We can do inv(A) to get inverse of matrix
 *We can do det(a) to get determinant
 *
 *A.col(0) returns the first column as a vector
 *
 */


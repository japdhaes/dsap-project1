#include <iostream>
#include <fstream>
#include "crystal.h"
#include <armadillo>
#include <ctime>
#include <vector>
#include "verletalgo.h"
#include "verletalgo2.h"
#include <sstream>
#include "printing.h"

using namespace arma;
using namespace std;



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

    int seed = -1;

    int nc=8;
    //latice parameter in unit Angstrom
    double b=5.28;
    double temperature=50;

    system("rm /home/jonathan/projectsFSAP/project1/project1/output/*.xyz");
    Printing p;
    Crystal crystal(nc, b, seed, temperature);
    p.printing(crystal);
    //VerletAlgo integrator(crystal);
    VerletAlgo2 integrator(crystal);
    for(int j=1; j<1500; j++){
        integrator.integrate();
        if(j%5==0){
            cout << "now in step " << j << " in the simulation" << endl;
            //cout << "nrofatomsfound "<<crystal.countAtoms()<< endl;
        }

        ofstream output;
        output.open(p.createname(j).c_str());

        //integrator.integrateWithCell();


        output << crystal << endl;

    }
    cout << "the boundary vector is " << crystal.boundary << endl;
    cout << "BC vector is " << crystal.vectorBC << endl;
    p.printvelocities(crystal);
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


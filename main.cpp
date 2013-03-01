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
#include <mpi.h>


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
    time_t tinit = time(0);

    int argc; char** argv; int numprocs, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int seed = -1;

    int nc=9;
    double h=0.0025;
    //latice parameter in unit Angstrom
    double b=5.28;
    double temperature=100;


    //200
    int nrofthermalizingsteps=200;

    //int nrofstepstotakemeasurements=2000;
    int nrofstepstotakemeasurements=1800;

    system("rm /home/jonathan/projectsFSAP/project1/project1/output/*.xyz");
    system("rm /home/jonathan/projectsFSAP/project1/project1/output/temperatures.txt");
    Printing p;
    Crystal crystal(nc, b, seed, temperature);
    p.printing(crystal);
    //VerletAlgo integrator(crystal);
    VerletAlgo2 integrator(crystal, h);
    ofstream measurements;
    measurements.open("/home/jonathan/projectsFSAP/project1/project1/output/measurementsT0.0025.txt");

    double temperatures[1800];

    measurements << "#1:time #2:temperature  #3:total energy #4:kinetic energy #5:potential energy #6:relative energy error"<<endl;
    for(int j=0; j<nrofthermalizingsteps; j++){
        //TURNED OFF THERMOSTAT!!!
        integrator.integrate(false);
        measurements << j << " " <<crystal.temperature()<<" " << integrator.crystall.energy<<" " << integrator.crystall.ke<<" " << integrator.crystall.pe <<endl;
        if(j%100==0){
            cout << "now in step " << j << " in the thermilisation phase" << endl;
            //cout << "nrofatomsfound "<<crystal.countAtoms()<< endl;
        }
        //ofstream output;
        //output.open(p.createname(j).c_str());
        //integrator.integrateWithCell();
        //output << crystal << endl;
    }
    for(int j=0; j<nrofstepstotakemeasurements; j++){

        integrator.integrate(false);
        temperatures[j]=crystal.temperature();
        measurements << j*integrator.h << " " <<crystal.temperature()<<" " << integrator.crystall.energy<<" " << integrator.crystall.ke<<" " << integrator.crystall.pe << " " << abs((integrator.crystall.energy - integrator.crystall.beginenergy)/integrator.crystall.beginenergy) << endl;
        //p.printvelocities(crystal, j);
        //integrator.integrate_noapprox();
        //tempoutput << crystal.temperature()<<" " << integrator.crystall.energy<<" " << integrator.crystall.ke<<" " << integrator.crystall.pe <<endl;
        if(j%50==0){
            cout << "now in step " << j << " doing measurements, i am processor  " << myrank << endl;
            //cout << "nrofatomsfound "<<crystal.countAtoms()<< endl;
        }
        ofstream output;
        output.open(p.createname(j).c_str());
        output << crystal << endl;
    }

    double avgtemp=0;
    for(int j=0; j<nrofstepstotakemeasurements; j++){
        avgtemp+=temperatures[j]/nrofstepstotakemeasurements;
    }

    cout << "average final temperature " << avgtemp << endl;

    double stddev=0;
    for(int j=0; j<nrofstepstotakemeasurements; j++){
        stddev+=(1.0/(nrofstepstotakemeasurements-1))*(temperatures[j]-avgtemp)*(temperatures[j]-avgtemp);
    }
    stddev=sqrt(stddev);

    cout << "stddev final temperature " << stddev << endl;

    //cout << "the boundary vector is " << crystal.boundary << endl;
    //cout << "BC vector is " << crystal.vectorBC << endl;
    //cout << "crystal has " << crystal.numberofatoms << " atoms "<< endl;
    time_t tdone = time(0);
    cout << "I AM PROCESSOR "<< myrank << " and i worked for " << tdone-tinit << " seconds " << endl;

    MPI_Finalize();
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


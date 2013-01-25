#ifndef ATOM_H
#define ATOM_H

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;
class Atom
{
    public:
        Atom(){}
        Atom(double x, double y, double z, double vx, double vy, double vz);
        rowvec phasevect;

        string chemelement;
        friend ostream& operator<<( ostream&, const Atom&);
};

#endif // ATOM_H

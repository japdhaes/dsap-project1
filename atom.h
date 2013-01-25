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
        Atom(rowvec);
        rowvec phasevect;

        string chemelement;
        double m;
        friend ostream& operator<<( ostream&, const Atom&);
};

#endif // ATOM_H

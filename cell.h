#ifndef CELL_H
#define CELL_H

#include <armadillo>
#include <iostream>
#include "atom.h"

using namespace std;
using namespace arma;

class Cell
{
    public:
        Cell(){};
        Cell(double x, double y, double z, double vx, double vy, double vz);
        rowvec phasevect;
        Atom* atoms;

        friend ostream& operator<<( ostream&, const Cell&);

};

#endif // CELL_H

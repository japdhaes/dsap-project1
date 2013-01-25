#include "cell.h"

Cell::Cell(double x, double y, double z, double vx, double vy, double vz)
{
    phasevect = rowvec(6);
    //initializing phasevector
    phasevect(0)=x;
    phasevect(1)=y;
    phasevect(2)=z;
    phasevect(3)=vx;
    phasevect(4)=vy;
    phasevect(5)=vz;

    //there are 4 atoms in a FCC lattice unit cell
    atoms=new Atom[4];
}

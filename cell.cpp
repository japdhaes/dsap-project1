#include "cell.h"

Cell::Cell(double x, double y, double z, double vx, double vy, double vz)
{
    double b = 5.260;
    phasevect = rowvec(6);
    //initializing phasevector
    phasevect(0)=x;
    phasevect(1)=y;
    phasevect(2)=z;
    phasevect(3)=vx;
    phasevect(4)=vy;
    phasevect(5)=vz;

    rowvec r1=rowvec(6), r2=rowvec(6), r3=rowvec(6), r4=rowvec(6);
    r1.fill(0); r2.fill(0); r3.fill(0); r4.fill(0);
    r2(1)=b/2; r2(2)=b/2;
    r3(0)=b/2; r3(1)=b/2;
    r4(0)=b/2; r4(2)=b/2;

    r1=phasevect+r1;
    r2=phasevect+r2;
    r3=phasevect+r3;
    r4=phasevect+r4;

    //there are 4 atoms in a FCC lattice unit cell
    atoms=new Atom[4];
    atoms[0]=Atom(r1);
    atoms[1]=Atom(r2);
    atoms[2]=Atom(r3);
    atoms[3]=Atom(r4);
}

ostream& operator<< (ostream& os , const Cell& cell){
    for (int i=0; i<4; i++){
        os << cell.atoms[i];
    }
    return os;
}

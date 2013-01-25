#include "atom.h"

Atom::Atom(double x, double y, double z, double vx, double vy, double vz)
{
    phasevect = rowvec(6);
    //initializing phasevector
    phasevect(0)=x;
    phasevect(1)=y;
    phasevect(2)=z;
    phasevect(3)=vx;
    phasevect(4)=vy;
    phasevect(5)=vz;

    this->chemelement="Ar";
}

ostream& operator<< (ostream& os , const Atom& atom){
    os << atom.chemelement;
    for(int i=0; i<6;i++){
        os<< " " << atom.phasevect(i);
    }
    return os;
}

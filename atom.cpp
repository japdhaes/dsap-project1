#include "atom.h"

Atom::Atom(rowvec phasearg)
{
    this->phasevect=phasearg;
    this->chemelement="Ar";
}

ostream& operator<< (ostream& os , const Atom& atom){
    os << atom.chemelement;
    for(int i=0; i<6;i++){
        os<< " " << atom.phasevect(i);
    }
    os<<endl;
    return os;
}

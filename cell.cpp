#include "cell.h"


Cell::Cell()
{
    this->first=NULL;
}

void Cell::insertElement(Atom *atom)
{
    ofstream debugging;
    debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglogcell.txt", ios::app);
    debugging << "entering insertElement with atom on position " << atom->getPosition()<< endl;
    //no atoms in this cell so far
    atom->currentcell=this;
    if(this->first==NULL){
        debugging << "no atom in the cell so far"<<endl;

        this->first=atom;
        atom->previousAtom=NULL;
        atom->nextAtom=NULL;
    }

    //already at least 1 atom in the cell
    else{
        debugging << "already atoms in the cell"<<endl;
        atom->nextAtom=this->first;
        this->first->previousAtom=atom;
        this->first=atom;
        atom->previousAtom=NULL;

    }
    if(this->isAtomInCell(atom)){
        debugging << "joepie atom in cell!"<< endl;
    }
    debugging.close();
}

void Cell::removeelement(Atom *atom){
    double x, y, z;
    findCell(atom, x, y, z);
}

void Cell::findCell(Atom* atom, double &x, double &y, double &z){

}

bool Cell::isAtomInCell(Atom *atom){
    vec3 r = atom->getPosition();
    int intpos[3];
    for(int j=0; j<3; j++){
        intpos[j]=int(r(j)/this->vectorBC(j));
    }
    if(x==intpos[0] && y==intpos[1] && z==intpos[2]){
        return true;
    }
    return false;
}

#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <iostream>
#include <vector>
#include "atom.h"
#include "cell.h"
#include "zignor.h"
#include "zigrandom.h"

using namespace std;

const double xunit  =3.405;
const double xunitSI=3.405e-10;
const double tunitSI=2.1569e-12;
const double vunitSI=xunitSI/tunitSI;
const double tempunit =119.74;


class Crystal
{
    public:
        //ofstream* debugging2;

        Crystal()
        {
        }
        Crystal(unsigned int nc, double _b, int& seed, double _temperature);

        vector<Atom*> allatoms;
        vector<vector<vector<Cell> > > allcells;

        int numberofatoms;
        int nc;
        double beginenergy;
        int counter;
        vec3 boundary;

        int countAtoms();
        friend ostream& operator<<( ostream&, const Crystal&);

        vec3 vectorBC;
        double b;
        double energy;
        double pe;
        double ke;
        double inittemp;
        void setvectorBC(double desiredwidth);
        void initializeAtoms(double _temperature);
        void addAllAtomsToCells();
        void initializeCells();
        void findCellOfAtom(Atom *atom, int &x, int &y, int &z);

        double temperature();
        void removeCrystalMomentum();
};

#endif // CRYSTAL_H

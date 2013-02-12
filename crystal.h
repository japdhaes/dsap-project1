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
const double tempunit =119.74;


class Crystal
{
    public:
        //ofstream debugging2;
        Crystal()
        {
        }
        Crystal(unsigned int nc, double _b, int& seed, double _temperature);

        vector<Atom*> allatoms;
        vector<vector<vector<Cell> > > allcells;

        int numberofatoms;
        int nc;
        vec3 boundary;

        int countAtoms();
        friend ostream& operator<<( ostream&, const Crystal&);

        vec3 vectorBC;
        double b;
        double inittemp;
        void setvectorBC(double desiredwidth);
        void initializeAtoms(double _temperature);
        void addAllAtomsToCells();
        void initializeCells();
};

#endif // CRYSTAL_H

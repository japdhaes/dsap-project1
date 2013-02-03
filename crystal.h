#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <iostream>
#include "zignor.h"
#include "zigrandom.h"
#include <vector>
#include "atom.h"
#include "lib.h"

using namespace std;

const double xunit  =3.405;

class Crystal
{
    public:
        Crystal(){}
        Crystal(unsigned int nc, double b, long& seed);

        vector<Atom*> allatoms;

        int numberofcells;
        int nc;
        vec3 boundary;

        friend ostream& operator<<( ostream&, const Crystal&);
};

#endif // CRYSTAL_H

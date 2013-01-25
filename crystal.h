#ifndef CRYSTAL_H
#define CRYSTAL_H

#include <iostream>
#include "cell.h"

using namespace std;

class Crystal
{
    public:
        Crystal(unsigned int nc, double b);
        Cell*** allcells;

        int numberofcells;
        int nc;

        friend ostream& operator<<( ostream&, const Crystal&);
};

#endif // CRYSTAL_H

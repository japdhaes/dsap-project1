#ifndef CRYSTAL_H
#define CRYSTAL_H
#include "cell.h"
class Crystal
{
    public:
        Crystal(unsigned int nc, double b);
        Cell*** allcells;
};

#endif // CRYSTAL_H

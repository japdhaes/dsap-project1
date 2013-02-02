#ifndef VERLETALGO_H
#define VERLETALGO_H

#include <armadillo>
#include "crystal.h"

//time unit of h =2.1569 *  10^3 fs

using namespace arma;
class VerletAlgo
{
public:
    VerletAlgo(Crystal &crystal);
    void integrate(Crystal &crystal);
    void integrateAtom(Atom *atom, vec3 boundvec);
    void calcAcceler(vec3 &position, vec3 &answer);

    double h;
    Crystal crystall;
    void boundCheck(vec3 &position, vec3 &boundvec);
};

#endif // VERLETALGO_H

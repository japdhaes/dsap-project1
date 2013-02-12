#ifndef VERLETALGO_H
#define VERLETALGO_H

#include <armadillo>
#include "crystal.h"

//time unit of h =2.1569 *  10^3 fs

using namespace arma;
class VerletAlgo
{
public:
    //ofstream debugging;
    VerletAlgo(Crystal &crystal);
    void integrate();
    void integrateAtom(Atom *atom, vec3 boundvec);

    double h;
    Crystal crystall;
    void boundCheck(vec3 &position, vec3 &boundvec);
    void calcAcceler(vec3 &position, vec3 &relvec, vec3 &answer);
    void updateVelocity(Atom *atom);
    void updateAcceler(Atom *atom);
    void updatePosition(Atom *atom, vec3 boundvec);
    vec3 findClosestPosition(vec3 &position, vec3 otherposition);
};

#endif // VERLETALGO_H

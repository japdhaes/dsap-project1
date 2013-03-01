#ifndef VERLETALGO2_H
#define VERLETALGO2_H

#include <armadillo>
#include "crystal.h"

const double cutoffacceleration=200;
//time unit of h =2.1569 *  10^3 fs

using namespace arma;
class VerletAlgo2
{
public:
    ofstream debugging;
    VerletAlgo2(Crystal &crystal, double _h);
    void integrate(bool);
    void integrate_noapprox();
    void integrateAtom(Atom *atom, vec3 boundvec);

    double h;
    Crystal crystall;

    void updatePosition(Atom *atom, vec3 &boundvec);
    void updateVelocity(Atom *atom);
    void updateAcceler(Atom *atom);
    void updateAccelerNoApprox(Atom *atom);

    void calcForce(Atom *atom, Atom *otheratom);
    void findXYZCellIndices(int *nrXYZ, int *nrX, int *nrY, int *nrZ);
    vec3 findClosestPosition(vec3 &position, vec3 &otherposition);
    void boundCheck(vec3 &position, vec3 &boundvec);
    double LJpotential(vec3 &relvec);
};

#endif // VERLETALGO2_H


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
    //ofstream debugging;
    VerletAlgo2(Crystal &crystal);
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
    void integrateWithCell();
    void integrateCell(int i, int j, int k, int imax, int jmax, int kmax);
    void integrateAtomToCell(Atom *integratingatom, int lfin, int mfin, int nfin);
    void findXYZCellIndices(int *nrXYZ, int *nrX, int *nrY, int *nrZ);
    void calcForce(Atom *atom, Atom *otheratom);
};

#endif // VERLETALGO2_H


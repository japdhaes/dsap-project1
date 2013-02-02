#include "verletalgo.h"
#define epsilon 119.8
#define sigma 3.405

VerletAlgo::VerletAlgo(Crystal &crystal)
{
    this->crystall=crystal;
    this->h=0.02;
}

void VerletAlgo::integrate(Crystal &crystal){
    int nc=crystal.nc;
    for(int i=0; i<crystal.allatoms.size(); i++){
        integrateAtom(crystal.allatoms[i]);
    }
}

void VerletAlgo::calcAcceler(vec3 &position, vec3 &answer){
    double r=norm(position,3);
    if(r==0){
        return;
    }
    answer-=24*(2*pow(r, -12)-pow(r, -6))*position/(r*r);

}

void VerletAlgo::integrateAtom(Atom *atom){
    vec3 position=atom->getPosition();
    vec3 velocity=atom->getVelocity();
    vec3 acceler=atom->getAcceler();
    vec3 acceleration; acceleration.fill(0);
    double withAcceleration=false;

    velocity+=0.5*acceler*this->h;
    position+=velocity*this->h;

    if(withAcceleration){
        for(int i=0; i<crystall.allatoms.size(); i++){
            vec3 relvec = atom->getPosition()-crystall.allatoms[i]->getPosition();
            calcAcceler(relvec, acceleration);
        }
    }
    else{
        acceleration.fill(0);
    }

    velocity+=0.5*acceleration*this->h;
}

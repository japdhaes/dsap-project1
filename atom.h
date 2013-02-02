#ifndef ATOM_H
#define ATOM_H

#define Xunit 3.405

#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;
class Atom
{
    public:
        Atom(){}
        Atom(vec r, vec v);
        vec3 getPosition();
        vec3 getVelocity();
        void setPosition(const vec3 &newPosition);
        void setVelocity(const vec3 &newVelocity);
        vec3 getAcceler();
        void setAcceler(const vec3 &newAcceler);


        string chemelement;
        double m;
    protected:
        vec3 position;
        vec3 velocity;
        vec3 acceler;

};

#endif // ATOM_H

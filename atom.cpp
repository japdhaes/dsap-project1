#include "atom.h"

Atom::Atom(vec position, vec v): position(position), velocity(v){
    this->chemelement="Ar";
}

vec3 Atom::getPosition(){
    return position;
}

vec3 Atom::getVelocity(){
    return velocity;
}

void Atom::setPosition(const vec3 &newPosition){
    this->position=newPosition;
}

void Atom::setVelocity(const vec3 &newVelocity){
    this->velocity=newVelocity;
}

vec3 Atom::getAcceler()
{
    return this->acceler;
}

void Atom::setAcceler(const vec3 &newAcceler)
{
    this->acceler=newAcceler;
}

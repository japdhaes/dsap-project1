#include "verletalgo2.h"

VerletAlgo2::VerletAlgo2(Crystal &crystal, double _h)
{
    this->debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog.txt");
    this->crystall=crystal;
    this->h=0.0025;
}

void VerletAlgo2::integrate(bool thermalize){
    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    crystall.energy=0;

    for(unsigned int i=0; i<crystall.allcells.size(); i++){
        for(unsigned int j=0; j<crystall.allcells.at(i).size();j++){
            for(unsigned int k=0; k<crystall.allcells.at(i).at(j).size(); k++){
                crystall.allcells.at(i).at(j).at(k).visited=false;
            }
        }
    }

    for(unsigned int i=0; i<crystall.allatoms.size(); i++){
        updateVelocity(crystall.allatoms[i]);
        updatePosition(crystall.allatoms[i], crystall.boundary);
    }

    for(unsigned int i=0; i<crystall.allcells.size(); i++){
        for(unsigned int j=0; j<crystall.allcells.at(i).size();j++){
            for(unsigned int k=0; k<crystall.allcells.at(i).at(j).size(); k++){
                crystall.allcells.at(i).at(j).at(k).visited=true;
                Atom *atom = crystall.allcells.at(i).at(j).at(k).first;
                while(atom!=NULL){
                    updateAcceler(atom);
                    atom=atom->nextAtom;
                }

            }
        }
    }

    for(unsigned int i=0; i<crystall.allatoms.size(); i++){
        updateAcceler(crystall.allatoms[i]);
    }
    crystall.pe=crystall.energy;
    //cout << "potential energy "<< potentialenergy<<endl;

    for(unsigned int i=0; i<crystall.allatoms.size(); i++){
        updateVelocity(crystall.allatoms[i]);

        //kinetic energy of the atom in the crystal
        crystall.energy+=0.5*dot(crystall.allatoms[i]->getVelocity(),crystall.allatoms[i]->getVelocity());
    }
    crystall.ke= crystall.energy-crystall.pe;
    if(thermalize){
        double tem=crystall.temperature();
        double ratio = crystall.inittemp/tem;
        for(int i=0; i<crystall.allatoms.size();i++){
            Atom *atom = crystall.allatoms[i];
            vec3 velocity = atom->getVelocity();
            atom->setVelocity(velocity*ratio);
        }
    }
    if(crystall.beginenergy==0){
        crystall.beginenergy=crystall.energy;
    }
}

void VerletAlgo2::integrate_noapprox(){
    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    for(unsigned int i=0; i<crystall.allatoms.size(); i++){
        updateVelocity(crystall.allatoms[i]);
        updatePosition(crystall.allatoms[i], crystall.boundary);
    }
    for(unsigned int i=0; i<crystall.allatoms.size(); i++){
        updateAccelerNoApprox(crystall.allatoms[i]);

    }
    for(unsigned int i=0; i<crystall.allatoms.size(); i++){
        vec3 acceler;
        acceler.fill(0);
        for(int j=0; j<crystall.allatoms.size();j++){
            for(int k=0;k<3;k++){
//                acceler(k)+=crystall.forces(i,j,k);
            }
        }
        crystall.allatoms[i]->setAcceler(acceler);
        updateVelocity(crystall.allatoms[i]);
    }
}

void VerletAlgo2::updateAcceler(Atom *atom){
    bool debugg=false;

    Atom *oneatom=atom;
    Atom* otheratom;
    while(oneatom!=NULL){
        otheratom=oneatom->nextAtom;
        while(otheratom!=NULL){
            calcForce(oneatom, otheratom);
            otheratom=otheratom->nextAtom;
        }
        oneatom=oneatom->nextAtom;
    }

    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    //indices of cell atom is in nrXYZ:
    int nrXYZ[3];
    vec3 r= atom->getPosition();
    for(int i=0; i<3; i++){
        nrXYZ[i]=int(r(i)/crystall.vectorBC(i));
    }
    int nrX[3], nrY[3], nrZ[3];
    findXYZCellIndices(nrXYZ, nrX, nrY, nrZ);
    //indices of all neighbouring cells are now in nrX, nrY and nrZ

    for(int i=0; i<3; i++){
        for(int j=0; j<3;j++){
            for(int k=0;k<3;k++){
                if(crystall.allcells.at(i).at(j).at(k).visited!=true){
                    otheratom = crystall.allcells.at(nrX[i]).at(nrY[j]).at(nrZ[k]).first;
                    while(otheratom!=NULL){
                        calcForce(atom, otheratom);
                        otheratom=otheratom->nextAtom;
                    }
                }
            }
        }

    }
}

void VerletAlgo2::updateAccelerNoApprox(Atom *atom){
    bool debugg=false;
    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    //indices of cell atom is in nrXYZ:
    int nrXYZ[3];
    vec3 r= atom->getPosition();
    for(int i=0; i<3; i++){
        nrXYZ[i]=int(r(i)/crystall.vectorBC(i));
    }
    int nrX[3], nrY[3], nrZ[3];
    findXYZCellIndices(nrXYZ, nrX, nrY, nrZ);
    //indices of all neighbouring cells are now in nrX, nrY and nrZ

    for(int i=0; i<crystall.allatoms.size();i++){
        Atom* otheratom = crystall.allatoms[i];
        calcForce(atom, otheratom);
    }

}

void VerletAlgo2::calcForce(Atom* atom, Atom* otheratom){
    int i = atom->number;
    int j = otheratom->number;

    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    //stop if integrating atom with itself, should normally not happen
    if(i==j){
        cout << "PROBLEM" <<endl;
        return;
    }

    vec3 position=atom->getPosition();
    vec3 othervec=otheratom->getPosition();
    vec3 closestvector = findClosestPosition(position, othervec);
    vec3 relvec = position - closestvector;
    vec3 relvec2 = position - othervec;
    double r2=dot(relvec,relvec);
    double r6=r2*r2*r2;
    double r12=r6*r6;

    vec3 oneacceler=atom->getAcceler();
    vec3 otheracceler=otheratom->getAcceler();
    for(int k=0; k<3; k++){
        double temp = 24*(2.0/r12-1.0/r6)*relvec(k)/r2;
        crystall.energy+=2.0*LJpotential( relvec2);
        if(temp>cutoffacceleration){
            temp=cutoffacceleration;
            //crystall.forces(i,j,k)=cutoffacceleration;
            //crystall.forces(j,i,k)=-1*cutoffacceleration;
        }
        else if(temp<-1*cutoffacceleration){
            temp=-cutoffacceleration;
//            crystall.forces(i,j,k)=-1*cutoffacceleration;
//            crystall.forces(j,i,k)=cutoffacceleration;
        }
//        else{
//            crystall.forces(i,j,k)=temp;
//            crystall.forces(j,i,k)=-temp;
//        }

        oneacceler(k)+=temp;
        otheracceler(k)-=temp;
    }
    atom->setAcceler(oneacceler);
    otheratom->setAcceler(otheracceler);
}

//this function finds the cellindices of the neighbouring cells
//it takes the minimal image convention into count
//it puts all X, Y and Z indices of the neighbouring cells in the last 3
//function arguments
//!!!!!!!!!!!!!!!!!function is debugged!!!!!!!!!!!!
void VerletAlgo2::findXYZCellIndices(int* nrXYZ, int* nrX, int* nrY, int* nrZ){
    //the maximum indices of cells in x, y and z direction
    int imax = crystall.allcells.size();
    int jmax=crystall.allcells.at(0).size();
    int kmax=crystall.allcells.at(0).at(0).size();

    int i=0;
    for(int l=nrXYZ[0]-1; l<nrXYZ[0]+2; l++){
        if(l<0){
            nrX[i]=l+imax;
        }
        else if(l>=imax){
            nrX[i]=l-imax;
        }
        else{
            nrX[i]=l;
        }
        i++;
    }
    i=0;
    for(int l=nrXYZ[1]-1; l<nrXYZ[1]+2; l++){
        if(l<0){
            nrY[i]=l+jmax;
        }
        else if(l>=jmax){
            nrY[i]=l-jmax;
        }
        else{
            nrY[i]=l;
        }
        i++;
    }
    i=0;
    for(int l=nrXYZ[2]-1; l<nrXYZ[2]+2; l++){
        if(l<0){
            nrZ[i]=l+kmax;
        }
        else if(l>=kmax){
            nrZ[i]=l-kmax;
        }
        else{
            nrZ[i]=l;
        }
        i++;
    }
}

vec3 VerletAlgo2::findClosestPosition(vec3 &position, vec3 &otherposition){
    vec3 answer; answer.fill(0);
    bool debugg=false;
    for(int i=0; i<3; i++){
        double projectionother = otherposition(i);
        double projectionpos = position(i);
        double l=this->crystall.boundary(i);
        double distance=l;
        for(int j=-1; j<2; j++){
            distance=abs(projectionpos-(projectionother+j*l));
            if(distance<=l/2){
                answer(i)=projectionother+j*l;
            }
        }
    }
    return answer;
}

void VerletAlgo2::updateVelocity(Atom *atom){
    vec3 velocity=atom->getVelocity();
    vec3 acceler=atom->getAcceler();
    velocity+=0.5*acceler*this->h;
    atom->setVelocity(velocity);
    /*if(norm(velocity,2)>10){
        cout << "atom number " << atom->number << " has velocity of " << norm(velocity,2)<<endl;
    }*/
}

void VerletAlgo2::updatePosition(Atom *atom, vec3& boundvec){
    vec3 position=atom->getPosition();
    vec3 velocity=atom->getVelocity();

    int nrXYZ[3];
    for(int i=0; i<3; i++){
        nrXYZ[i]=int(position(i)/crystall.vectorBC(i));
    }

    position+=velocity*this->h;
    //realposition updated for the diffusion
    atom->realposition+=velocity*this->h;
    boundCheck(position, boundvec);
    atom->setPosition(position);

    //if atom is not anymore in its cell
    if(!(crystall.allcells.at(nrXYZ[0]).at(nrXYZ[1]).at(nrXYZ[2]).isAtomInCell(atom))){
        crystall.allcells.at(nrXYZ[0]).at(nrXYZ[1]).at(nrXYZ[2]).removeelement(atom);
        int x,y,z;
        crystall.findCellOfAtom(atom, x, y, z);
        crystall.allcells.at(x).at(y).at(z).insertElement(atom);
    }
}


void VerletAlgo2::boundCheck(vec3 &position, vec3 &boundvec){
    for(int i=0; i<3; i++){
        while(position(i)<0){
            //debugging << "summing boundvec at i="<<i<< " position(i)=" << position(i);
            position(i)+=boundvec(i);
            //debugging << "new position(i) " << position(i) << endl;
        }
        while(position(i)>boundvec(i)){
            position(i)-=boundvec(i);
        }
    }
}

double VerletAlgo2::LJpotential(vec3 &relvec){
    double r=norm(relvec,2);
    double r2=r*r;
    double r6=r2*r2*r2;
    double r12=r6*r6;

    return 4*(1.0/r12-1.0/r6);
}

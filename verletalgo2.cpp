#include "verletalgo2.h"

VerletAlgo2::VerletAlgo2(Crystal &crystal)
{
    //this->debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog.txt");
    this->crystall=crystal;
    this->h=0.006;
}

void VerletAlgo2::integrate(){
    bool debugg=false;
    Crystal crystal = this->crystall;
    crystal.forces.fill(0);
    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        updateVelocity(crystal.allatoms[i]);
        updatePosition(crystal.allatoms[i], crystal.boundary);
    }
    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        if(debugg){
        }
        updateAcceler(crystal.allatoms[i]);

    }
    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        updateVelocity(crystal.allatoms[i]);

    }

}
void VerletAlgo2::updateAcceler(Atom *atom){
    bool debugg=false;
    ofstream debugging;
    debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    //indices of cell atom is in:
    int nrXYZ[3];
    vec3 r= atom->getPosition();
    for(int i=0; i<3; i++){
        //nrXYZ[i]=int(crystall.boundary(i)/crystall.vectorBC(i));
        nrXYZ[i]=int(r(0)/crystall.vectorBC(0));
    }
    debugging << "indices of current cell"<< endl;
    for(int i=0; i<3; i++){
        debugging << nrXYZ[i];
    }
    debugging << endl<<endl<< endl;

    int nrX[3], nrY[3], nrZ[3];
    findXYZCellIndices(nrXYZ, nrX, nrY, nrZ);



    debugging << "indices of neighbouring cells"<< endl;
    for(int i=0; i<3; i++){
        for(int j=0; j<3;j++){
            for(int k=0;k<3;k++){
                debugging << nrX[i]<< nrY[j]<< nrZ[k]<<endl;
            }
        }

    }
    debugging << endl<<endl<< endl;


    //atom->setAcceler(acceler);
}

//this function finds the cellindices of the neighbouring cells
//it takes the minimal image convention into count
//it puts all X, Y and Z indices of the neighbouring cells in the last 3
//function arguments
void VerletAlgo2::findXYZCellIndices(int* nrXYZ, int* nrX, int* nrY, int* nrZ){
    //the maximum indices of cells in x, y and z direction
    int imax = crystall.allcells.size();
    int jmax=crystall.allcells.at(0).size();
    int kmax=crystall.allcells.at(0).at(0).size();


    int lfin, mfin, nfin;
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

vec3 VerletAlgo2::findClosestPosition(vec3 &position, vec3 otherposition){
    vec3 answer; answer.fill(0);
    bool debugg=false;
    if(debugg){
        /*debugging << __FUNCTION__ << endl;
        debugging << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        debugging << "position " << position << endl;
        debugging << "otherposition " << otherposition << endl;*/
    }
    for(int i=0; i<3; i++){
        if(debugg){
            //debugging << endl << endl << "integer " << i << endl;
            //debugging << "===================" << endl;
        }
        double projectionother = otherposition(i);
        double projectionpos = position(i);
        double l=this->crystall.boundary(i);
        double distance=l;
        for(int j=-1; j<2; j++){
            //debugging << "integer j" << j << endl;
            distance=abs(projectionpos-(projectionother+j*l));
            /*debugging << "distance " << distance << endl;
            debugging << "projectionpos " << projectionpos << endl;
            debugging << "projectionother " << projectionother << endl;
            debugging << "suggested new position "<< projectionother+j*l << endl;
            debugging << "l " << l << endl;*/
            if(distance<=l/2){
                //debugging << "UPDATING ANSWER with " << projectionother+j*l<< endl;
                answer(i)=projectionother+j*l;
            }
        }
    }
    /*debugging << "answer " << answer << endl;
    debugging << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;*/
    return answer;
}

void VerletAlgo2::updateVelocity(Atom *atom){
    vec3 velocity=atom->getVelocity();
    vec3 acceler=atom->getAcceler();
    velocity+=0.5*acceler*this->h;
    //if(norm(velocity,2)>5){
        /*debugging << "after updating velocity " << endl;
        debugging << "velocity" << velocity << endl;
        debugging << "acceler" <<acceler << endl;*/
    //}
    atom->setVelocity(velocity);
}

void VerletAlgo2::updatePosition(Atom *atom, vec3 boundvec){
    vec3 position=atom->getPosition();
    vec3 velocity=atom->getVelocity();
    position+=velocity*this->h;
    boundCheck(position, boundvec);
    atom->setPosition(position);
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

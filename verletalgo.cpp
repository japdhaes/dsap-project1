#include "verletalgo.h"
#define epsilon 119.8
#define sigma 3.405

VerletAlgo::VerletAlgo(Crystal &crystal)
{
    //this->debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog.txt");
    this->crystall=crystal;
    this->h=0.01;
}

void VerletAlgo::integrate(){
    bool debugg=false;
    Crystal crystal = this->crystall;
    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        updateVelocity(crystal.allatoms[i]);
        updatePosition(crystal.allatoms[i], crystal.boundary);
    }
    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        if(debugg){
            /*debugging << "INTEGRATING NEW ATOM"<<endl;
            debugging << "----------------------------------------------------------------------" << endl;
            debugging << "----------------------------------------------------------------------" << endl;
            debugging << "----------------------------------------------------------------------" << endl;
            debugging << "----------------------------------------------------------------------" << endl;
            debugging << "----------------------------------------------------------------------" << endl;
            debugging << "----------------------------------------------------------------------" << endl;
            debugging << "----------------------------------------------------------------------" << endl;*/
        }
        updateAcceler(crystal.allatoms[i]);

    }
    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        updateVelocity(crystal.allatoms[i]);
    }
}

void VerletAlgo::calcAcceler(vec3 &position, vec3 &othervec, vec3 &answer){
    vec3 closestvector = findClosestPosition(position, othervec);
    vec3 relvec = closestvector-position;
    double r=norm(relvec,2);
    /*debugging << "before answer" << endl;
    debugging << "answer:" <<answer << endl;
    debugging << "r:" << r << endl;
    debugging << "position:" << position << endl;
    debugging << "otherposition" << othervec << endl;
    debugging << "closestvector: "<<closestvector<< endl;*/
    if(r<0.001){
        return;
    }

    if(r<0.01){
        /*debugging << "too close encounter" << endl;
        debugging << "position" << position << endl;
        debugging << "other position " << othervec << endl;
        debugging << "closest position" << closestvector << endl;
        debugging << "relative vector norm " << r << endl;*/
    }

    /*cout << "relvecnorm " << r << endl;
    cout << "r-12 " << pow(r, -12) << endl;
    cout << "r-6 " << pow(r, -6) << endl;
    cout << "1 acceleration calculated " << 24*(2*pow(r, -12)-pow(r, -6))*position/(r*r) << endl;*/


    if(norm(answer,2)>100){
        //debugging << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    }
    answer-=24*(2*pow(r, -12)-pow(r, -6))*relvec/(r*r);
    //debugging << "after answer is updated" << answer << endl;

}

void VerletAlgo::updateAcceler(Atom *atom){
    bool debugg=false;
    vec3 position=atom->getPosition();
    vec3 acceler; acceler.fill(0);
    for(unsigned int i=0; i<crystall.allatoms.size(); i++){
        vec3 othervec = crystall.allatoms[i]->getPosition();
        //cout << "relvec " << relvec << endl;
        if(debugg){
            //debugging << "INTEGRATING ATOM "<< *atom << " with other atom " << *(crystall.allatoms[i])<<endl;
        }
        calcAcceler(position, othervec, acceler);
    }
    //debugging << norm(acceler,2)<< endl;
    if(norm(acceler,2)>500){
        for(int i=0; i<3; i++){
            //cout << "changed acceleration "<< endl;
            //acceler(i)/=10;
        }
    }
    atom->setAcceler(acceler);
}

vec3 VerletAlgo::findClosestPosition(vec3 &position, vec3 otherposition){
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

void VerletAlgo::updateVelocity(Atom *atom){
    vec3 velocity=atom->getVelocity();
    vec3 acceler=atom->getAcceler();
    /*debugging << "before updating velocity " << endl;
    debugging << "velocity" << velocity << endl;*/
    velocity+=0.5*acceler*this->h;
    //if(norm(velocity,2)>5){
        /*debugging << "after updating velocity " << endl;
        debugging << "velocity" << velocity << endl;
        debugging << "acceler" <<acceler << endl;*/
    //}
    atom->setVelocity(velocity);
}

void VerletAlgo::updatePosition(Atom *atom, vec3 boundvec){
    vec3 position=atom->getPosition();
    vec3 velocity=atom->getVelocity();
    position+=velocity*this->h;
    boundCheck(position, boundvec);
    atom->setPosition(position);
}


void VerletAlgo::boundCheck(vec3 &position, vec3 &boundvec){
    for(int i=0; i<3; i++){
        if(position(i)<0){
            //debugging << "summing boundvec at i="<<i<< " position(i)=" << position(i);
            position(i)+=boundvec(i);
            //debugging << "new position(i) " << position(i) << endl;
        }
        else if(position(i)>boundvec(i)){
            position(i)-=boundvec(i);
        }
    }
}

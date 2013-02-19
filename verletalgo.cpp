#include "verletalgo.h"
#define epsilon 119.8
#define sigma 3.405

VerletAlgo::VerletAlgo(Crystal &crystal)
{
    //this->debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog.txt");
    this->crystall=crystal;
    this->h=0.006;
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

void VerletAlgo::integrateWithCell(){
    bool debugg=false;
    ofstream debugging;
    debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt", ios::app);
    Crystal crystal = this->crystall;
    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        updateVelocity(crystal.allatoms[i]);
        updatePosition(crystal.allatoms[i], crystal.boundary);
    }
    int imax = crystal.allcells.size();
    int jmax=crystal.allcells.at(0).size();
    int kmax=crystal.allcells.at(0).at(0).size();
    for(int i=0; i<imax; i++){
        for(int j=0; j<jmax; j++){
            for(int k=0; k<kmax; k++){
                /*debugging << "integrating atoms in cell (i,j,k)=" << i << " ";
                debugging << j << " " << k << " ";
                debugging << " to all atoms in cells ";*/

                integrateCell(i,j,k,imax,jmax,kmax);


                //cout << "i "<<i << " j "<< j << " k " << k << endl;
            }
        }
    }

    for(unsigned int i=0; i<crystal.allatoms.size(); i++){
        updateVelocity(crystal.allatoms[i]);
        vec3 acceler = crystal.allatoms[i]->getAcceler();
        //cout << "do i get here" << endl;
        if(norm(acceler,2)<3){
            //cout << " does this never happen? "<< endl;
            debugging << "particle i "<<i << " not accelerated"<< endl;
            debugging << acceler << endl;
        }
    }
}

void VerletAlgo::integrateCell(int i, int j, int k, int imax, int jmax, int kmax){
    int lfin, mfin, nfin;
    for(int l=i-1; l<i+2; l++){
        if(l<0){
            lfin=l+imax;
        }
        else if(l>=imax){
            lfin=l-imax;
        }
        else{
            lfin=l;
        }
        for(int m=j-1; m<j+2;m++){
            if(m<0){
                mfin=m+jmax;
            }
            else if(m>=jmax){
                mfin=m-jmax;
            }
            else{
                mfin=m;
            }
            for(int n=k-1; n<k+2; n++){
                //cout << "n "<< n << " k "<< k << " kmax "<< kmax << endl;
                if(n<0){
                    nfin=n+kmax;
                }
                else if(n>=kmax){
                    nfin=n-kmax;
                }
                else{
                    nfin=n;
                }

                // lfin, mfin and nfin are the final (l,m,n) indices
                // of the cells we calculate the acceleration of the
                // atom at

                Atom *integratingatom = (crystall.allcells.at(i).at(j).at(k)).first;
                while(integratingatom!=NULL){
                    //cout << "do i get here? before"<< endl;
                    //cout << "lfin "<< lfin << " mfin "<< mfin << " nfin " << nfin << endl;
                    integrateAtomToCell(integratingatom, lfin, mfin, nfin);
                    //cout << "do i get here? after"<< endl;

                    /*if(! crystall.allatoms.at(i).at(j).at(k).isAtomInCell(integratingatom)){
                        crystall.allatoms.at(i).at(j).at(k).removeelement(integratingatom);
                        crystall.findCellOfAtom(integratingatom, );

                    }*/
                    integratingatom=integratingatom->nextAtom;
                }
            }
        }
    }

}

void VerletAlgo::integrateAtomToCell(Atom *atom, int lfin, int mfin, int nfin){
    vec3 position=atom->getPosition();
    vec3 acceler; acceler.fill(0);
    Atom *otheratom = crystall.allcells.at(lfin).at(mfin).at(nfin).first;
    while(otheratom!=NULL){
        vec3 othervec= otheratom->getPosition();
        calcAcceler(position, othervec, acceler);
        otheratom=otheratom->nextAtom;
    }
    atom->setAcceler(acceler);
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

    for(int i=0; i<3; i++){
        if((24*(2*pow(r, -12)-pow(r, -6))*relvec(i)/(r*r))>400){
            answer(i)-=400;
        }
        else if((24*(2*pow(r, -12)-pow(r, -6))*relvec(i)/(r*r))<-400){
            answer(i)+=400;
        }
        else{
            answer(i)-=24*(2*pow(r, -12)-pow(r, -6))*relvec(i)/(r*r);
        }
    }
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

#include "crystal.h"
#include "zignor.c"
#include "zigrandom.c"
#define Xunit 3.405


//Comment for future purposes:
//Mean is indeed approximate 0
//Normal velocities should be around 1.44 A/ps=144 m/s
Crystal::Crystal(unsigned int _nc, double _b, int &seed, double _temperature)
{
    //this->debugging2.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt");
    this->nc=_nc;
    this->b=_b;
    this->inittemp=_temperature/tempunit;
    this->boundary << nc*b/xunit << nc*b/xunit << nc*b/xunit;
    this->numberofatoms=4*nc*nc*nc;
    this->energy=0;
    this->beginenergy=0;
    this->counter=0;
    this->pressure=0;

    this->volume=this->boundary(0)*this->boundary(1)*this->boundary(2);
    this->density=numberofatoms/volume;
    this->msqdplm=0;


    RanNormalSetSeedZigVec(&seed, 100);

    setvectorBC(3);

    this->initializeAtoms(_temperature);
    this->removeCrystalMomentum();
    this->initializeCells();
    this->addAllAtomsToCells();

}

void Crystal::removeCrystalMomentum(){
    vec3 velocitylattice;
    velocitylattice.zeros();
    for(int i=0; i<this->allatoms.size(); i++){
        Atom *atom = this->allatoms[i];
        vec3 vel = atom->getVelocity();
        velocitylattice+=vel;
    }

    velocitylattice/=this->allatoms.size();
    for(int i=0; i<this->allatoms.size(); i++){
        Atom *atom = this->allatoms[i];
        vec3 vel = atom->getVelocity();
        vel-=velocitylattice;
        atom->setVelocity(vel);
    }
}

void Crystal::findCellOfAtom(Atom *atom, int &x, int &y, int &z){
    vec3 r= atom->getPosition();
    x=int(r(0)/this->vectorBC(0));
    y=int(r(1)/this->vectorBC(1));
    z=int(r(2)/this->vectorBC(2));
}

double Crystal::temperature()
{
    double temp=0;
    for(int i=0; i<allatoms.size();i++){
        Atom *atom=allatoms[i];
        //temp+=0.5*dot(atom->getVelocity(), atom->getVelocity());
        vec3 vel = atom->getVelocity();
        double v = norm(vel,2);
        temp+=0.5*v*v;
    }
    return temp*2.0/3/this->numberofatoms;
}

void Crystal::addAllAtomsToCells(){
    //ofstream debugging;
    //debugging.open("/home/jonathan/projectsFSAP/project1/project1/debuglog2.txt");
    //debugging << "boundary" <<this->boundary << endl << "vectorBC" << this->vectorBC << endl;
    for(unsigned int i=0; i<this->allatoms.size(); i++){
        vec3 r = (this->allatoms[i])->getPosition();
        Atom* atom = this->allatoms[i];
        int x, y, z;
        findCellOfAtom(atom, x, y, z);
        this->allcells.at(x).at(y).at(z).insertElement(atom);

    }
}

void Crystal::initializeCells(){
    vec3 vectorBC=this->vectorBC;

    int nrXYZ[3];
    for(int i=0; i<3; i++){
        nrXYZ[i]=round(boundary(i)/vectorBC(i));
    }

    this->allcells.resize(nrXYZ[0]);
    for(int i=0; i<nrXYZ[0]; i++){
        this->allcells[i].resize(nrXYZ[1]);
        for(int j=0; j<nrXYZ[1]; j++){
            this->allcells[i][j].resize(nrXYZ[2]);
            for(int k=0; k<nrXYZ[2]; k++){
                allcells.at(i).at(j).at(k).x=i;
                allcells.at(i).at(j).at(k).y=j;
                allcells.at(i).at(j).at(k).z=k;
                allcells.at(i).at(j).at(k).vectorBC=this->vectorBC;
            }
        }
    }
}

void Crystal::initializeAtoms(double _temperature){
    double tem = _temperature/tempunit;
    int l=0;
    //constructing a nc x nc x nc structure of cells
    for(int i=0; i< nc ;++i){
        for(int j=0; j< nc; ++j){
            for(int k=0; k<nc ;++k){
                vec3 positioncell, r1, r2, r3, r4, v1, v2, v3, v4;
                positioncell.zeros(); r1.zeros(); r2.zeros(); r3.zeros(); r4.zeros();
                positioncell << i*b/xunit << j*b/xunit << k*b/xunit;
                r2 << b/xunit/2 << b/xunit/2 << 0;
                r3 << 0  << b/xunit/2<< b/xunit/2;
                r4 << b/xunit/2 << 0 << b/xunit/2;

                r1=positioncell+r1;
                r2=positioncell+r2;
                r3=positioncell+r3;
                r4=positioncell+r4;


                bool gaussiandistribution=true;
                if(gaussiandistribution){
                    v1.fill(0); v2.fill(0); v3.fill(0); v4.fill(0);
                    v1 << DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem);
                    v2 << DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem);
                    v3 << DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem);
                    v4 << DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem)<< DRanNormalZigVec()*sqrt(tem);
                }
                else{

                    v1=randu(3); v1-=0.5; v1*=2*tem;
                    v2=randu(3); v2-=0.5; v2*=2*tem;
                    v3=randu(3); v3-=0.5; v3*=2*tem;
                    v4=randu(3); v4-=0.5; v4*=2*tem;
                }

                Atom* atom1=new Atom(r1, v1);
                Atom* atom2=new Atom(r2, v2);
                Atom* atom3=new Atom(r3, v3);
                Atom* atom4=new Atom(r4, v4);

                atom1->number=l; l++;
                atom2->number=l; l++;
                atom3->number=l; l++;
                atom4->number=l; l++;

                this->allatoms.push_back(atom1);
                this->allatoms.push_back(atom2);
                this->allatoms.push_back(atom3);
                this->allatoms.push_back(atom4);
            }
        }
    }
}

ostream& operator<< (ostream& os , const Crystal& crystal){
    //first write down the total number of atoms in the simulated crystal
    os << crystal.numberofatoms << endl;
    os << "Some comments here" << endl;
    vector<Atom*> myvector=crystal.allatoms;
    vector<Atom*>::iterator it=myvector.begin();
    for(it=myvector.begin(); it!=myvector.end(); ++it){
        vec3 position=(*it)->getPosition();
        vec3 velocity=(*it)->getVelocity();
        os << (*it)->chemelement << " "<< position(0)*xunit<< " " << position(1)*xunit << " " << position(2)*xunit << " " << velocity(0)<< " " << velocity(1)<< " "<< velocity(2)<< endl;
        //os << (*it)->chemelement << " "<< position(0)<< " " << position(1) << " " << position(2) << " " << velocity(0)<< " " << velocity(1)<< " "<< velocity(2)<< endl;
    }
    return os;
}

int Crystal::countAtoms(){
    int nr=0;
    for(int i=0; i<this->allcells.size();i++){
        for(int j=0; j<this->allcells.size();j++){
            for(int k=0; k <this->allcells.size();k++){
                Atom *atom = this->allcells.at(i).at(j).at(k).first;
                while(atom!=NULL){
                    nr++;
                    atom=atom->nextAtom;
                }
            }
        }
    }
    return nr;
}

//This function will find a width for cells in the crystal
//with a width as close to the desiredwidth as possible
//so that we still have an integer number of equally sized cells
//in the crystal
//if you want to have boxes with width of 3 sigma, set desiredwidth to 3
//this function works in computer units of sigma!
void Crystal::setvectorBC(double desiredwidth)
{
    this->vectorBC.fill(0);
    for(int i=0; i<3; i++){
        int j = int(boundary(i) / desiredwidth);
        this->vectorBC(i)= this->boundary(i)/j;
    }
}

void Crystal::radialDistFunction(){

    double maxlength=this->boundary(0)*0.5*sqrt(3);
    double unit=maxlength*0.01;
    double distribution[100];
    for(int i=0; i<100; i++){
        distribution[i]=0.0;
    }
    for(int i=0; i<this->allatoms.size();i++){
        for(int j=i+1; j<this->allatoms.size();j++){
            Atom *atom = this->allatoms[i];
            Atom *otheratom = this->allatoms[j];
            vec3 position=atom->getPosition();
            vec3 otherposition = this->findClosestPosition(position, otheratom->getPosition());
            double norm=0.0;
            for(int k=0; k<3; k++){
                norm+=(position(k)-otherposition(k))*(position(k)-otherposition(k));
            }
            norm=sqrt(norm);
            distribution[int(norm/unit)]+=2;
        }
    }
    double average=0.0;
    for(int i=0; i<100; i++){
        average+=distribution[i];
    }
    average/=100;
    ofstream measurements;
    measurements.open("/home/jonathan/projectsFSAP/project1/project1/output/radialdistributionsolid.txt");
    for(int i=0; i<100; i++){
        measurements<< i*unit <<" " <<distribution[i]/average<<endl;
    }
    measurements.close();

}

vec3 Crystal::findClosestPosition(vec3 position, vec3 otherposition){
    vec3 answer; answer.fill(0);
    bool debugg=false;
    for(int i=0; i<3; i++){
        double projectionother = otherposition(i);
        double projectionpos = position(i);
        double l=this->boundary(i);
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

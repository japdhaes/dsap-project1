#include "crystal.h"
#include "zignor.c"
#include "zigrandom.c"
#define Tunit 119.74
#define Xunit 3.405


//Comment for future purposes:
//Mean is indeed approximate 0
//Normal velocities should be around 1.44 A/ps=144 m/s
Crystal::Crystal(unsigned int _nc, double b, int &seed)
{
    this->nc=_nc;
    this->numberofcells=nc*nc*nc;
    this->boundary << nc*b/xunit << nc*b/xunit << nc*b/xunit;

    double T=100, tem = T/Tunit;
    RanNormalSetSeedZigVec(&seed, 100);

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

                v1.fill(0); v2.fill(0); v3.fill(0); v4.fill(0);
                v1 << DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem);
                v2 << DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem);
                v3 << DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem);
                v4 << DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem)<< DRanNormalZigVec()*sqrt(3*tem);

                Atom* atom1=new Atom(r1, v1);
                Atom* atom2=new Atom(r2, v2);
                Atom* atom3=new Atom(r3, v3);
                Atom* atom4=new Atom(r4, v4);

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
    os << 4*crystal.nc*crystal.nc*crystal.nc << endl;
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
    for(unsigned int i=0; i<this->allatoms.size(); i++){
        bool atomisok=true;
        for(int j=0; j<3; j++){
            vec3 bounds=this->boundary;
            vec3 position=this->allatoms[i]->getPosition();
            if(position(j)<0 ||position(j)>bounds(j)){
                atomisok=false;
            }
        }
        if(atomisok){
            nr++;
        }
    }
    return nr;
}

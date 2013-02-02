#include "crystal.h"
#define Tunit 119.74
#define Xunit 3.405


//Comment for future purposes:
//Mean is indeed approximate 0
//Normal velocities should be around 1.44 A/ps=144 m/s
Crystal::Crystal(unsigned int _nc, double b, long &seed)
{
    this->nc=_nc;
    this->numberofcells=nc*nc*nc;

    //RanNormalSetSeedZigVec(&seed, 100);
    double T=100, tem = T/Tunit;
    //constructing a nc x nc x nc structure of cells
    for(int i=0; i< nc ;++i){
        for(int j=0; j< nc; ++j){
            for(int k=0; k<nc ;++k){
                vec3 positioncell, r1, r2, r3, r4, v1, v2, v3, v4;
                positioncell.zeros(); r1.zeros(); r2.zeros(); r3.zeros(); r4.zeros();
                positioncell << i*b/Xunit << j*b/Xunit << k*b/Xunit;
                r2 << b/Xunit/2 << b/Xunit/2 << 0;
                r3 << 0  << b/Xunit/2<< b/Xunit/2;
                r4 << b/Xunit/2 << 0 << b/Xunit/2;

                r1=positioncell+r1;
                r2=positioncell+r2;
                r3=positioncell+r3;
                r4=positioncell+r4;

                v1.fill(0); v2.fill(0); v3.fill(0); v4.fill(0);
                double vx[2], vy[2], vz[2];
                gausran(vx, vx+1, &seed);
                gausran(vy, vy+1, &seed);
                gausran(vz, vz+1, &seed);
                v1 << vx[0]*sqrt(3*tem)<< vy[0]*sqrt(3*tem)<< vz[0]*sqrt(3*tem);
                v2 << vx[1]*sqrt(3*tem)<< vy[1]*sqrt(3*tem)<< vz[1]*sqrt(3*tem);
                gausran(vx, vx+1, &seed);
                gausran(vy, vy+1, &seed);
                gausran(vz, vz+1, &seed);
                v3 << vx[0]*sqrt(3*tem)<< vy[0]*sqrt(3*tem)<< vz[0]*sqrt(3*tem);
                v4 << vx[1]*sqrt(3*tem)<< vy[1]*sqrt(3*tem)<< vz[1]*sqrt(3*tem);

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
    /*for(int i=0; i<crystal.nc; i++){
        for(int j=0; j<crystal.nc; j++){
            for(int k=0; k<crystal.nc; k++){
                os << crystal.allcells[i][j][k] << endl;
            }
        }
    }*/
    vector<Atom*> myvector=crystal.allatoms;
    vector<Atom*>::iterator it=myvector.begin();
    int i=0;
    for(it=myvector.begin(); it!=myvector.end(); ++it){
        vec3 position=(*it)->getPosition();
        vec3 velocity=(*it)->getVelocity();
        os << (*it)->chemelement << " "<< position(0)*3.405 << " " << position(1)*3.405 << " " << position(2)*3.405 << " " << velocity(0)<< " " << velocity(1)<< " "<< velocity(2)<< endl;
        i++;
    }
    return os;
}

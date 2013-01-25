#include "crystal.h"

Crystal::Crystal(unsigned int _nc, double b)
{
    this->nc=_nc;
    this->numberofcells=nc*nc*nc;
    vec a=randn<vec>(3*4*8*8*8);
    unsigned int index=0;

    //constructing a nc x nc x nc structure of cells
    this->allcells=new Cell**[nc];
    for(int i=0; i< nc ;++i){
        this->allcells[i]=new Cell*[nc];
        for(int j=0; j< nc; ++j){
            this->allcells[i][j]=new Cell[nc];
            for(int k=0; k<nc ;++k){
                //no velocity for the position of the cells
                this->allcells[i][j][k]=Cell(i*b, j*b, k*b, 0, 0, 0);
                for(int l=0; l<4; l++){
                    for(int m=0; m<3; m++){
                        this->allcells[i][j][k].atoms[l].phasevect(m+3)=a(index);
                        index++;
                    }
                }
            }
        }
    }
}

ostream& operator<< (ostream& os , const Crystal& crystal){
    //first write down the total number of atoms in the simulated crystal
    os << 4*crystal.nc*crystal.nc*crystal.nc << endl;
    os << "Some comments here" << endl;
    for(int i=0; i<crystal.nc; i++){
        for(int j=0; j<crystal.nc; j++){
            for(int k=0; k<crystal.nc; k++){
                os << crystal.allcells[i][j][k] /*<< endl*/;
            }
        }
    }
    return os;
}

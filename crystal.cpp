#include "crystal.h"

Crystal::Crystal(unsigned int nc, double b)
{
    //constructing a nc x nc x nc structure of cells
    this->allcells=new Cell**[nc];
    for(int i=0; i< nc ;++i){
        this->allcells[i]=new Cell*[nc];
        for(int j=0; j< nc; ++j){
            this->allcells[i][j]=new Cell[nc];
            for(int k=0; k<nc ;++k){
                //no velocity for the position of the cells
                this->allcells[i][j][k]=Cell(i*b, j*b, k*b, 0, 0, 0);
            }
        }
    }
}

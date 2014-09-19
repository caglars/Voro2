//
//  CSCalculation.h
//  Voro2
//
//  Created by Caglar Sinayuc on 12/09/14.
//  Copyright (c) 2014 caglars. All rights reserved.
//

#ifndef __Voro2__CSCalculation__
#define __Voro2__CSCalculation__

#include <iostream>
#include "voro++.hh"
using namespace voro;

class Calculation {
    public:
        double *pressure;
        double *rightHandSide;
        double **coefficient;
        double *positionX;
        double *positionY;
        double *positionZ;
        double *cellVolume;
        void simRun(container& aCon);
        void calculate();
        double distance(double x1, double y1, double z1, double x2, double y2, double z2);
    
};




#endif /* defined(__Voro2__CSCalculation__) */

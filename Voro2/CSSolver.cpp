//
//  CSSolver.cpp
//  Voro2
//
//  Created by Caglar Sinayuc on 14/09/14.
//  Copyright (c) 2014 caglars. All rights reserved.
//

#include "CSSolver.h"
#include <math.h>

double * Solution::simpleSolver(int arraySize, double **coefficientArray, double *rightHandSideArray, double *initialValueArray) {
    
    int firstCounter = 0, secondCounter = 0, counter = 0;
	int LHSIndex = 0, parameterCounter = 0;
    
    double sumError = 0.0;
    double leftHandSideValue = 0.0;
    double solution = 0.0;
    double *leftHandSideArray;
    double *oldValueArray;
    double *valueArray;
    valueArray = new double[arraySize];
    oldValueArray = new double[arraySize];
    leftHandSideArray = new double[arraySize];
    
    for (parameterCounter=0; parameterCounter<arraySize; parameterCounter++)
	{
        valueArray[parameterCounter] = initialValueArray[parameterCounter];
        oldValueArray[parameterCounter] = valueArray[parameterCounter];

	}
    
    counter = 0;

    do {
        printf("Counter: %d\n", counter);
        for (firstCounter = 0; firstCounter < arraySize; firstCounter++) {
            LHSIndex = 0;
            sumError = 0.0;
            
            leftHandSideArray[firstCounter] = 0.0;
            for (secondCounter = 0; secondCounter < arraySize-1; secondCounter++) {
                LHSIndex = (firstCounter+secondCounter+1)%arraySize;
                leftHandSideValue = leftHandSideArray[firstCounter]
                                    + coefficientArray[firstCounter][LHSIndex]*valueArray[LHSIndex];
                //printf("Left Hand Side Value: %f\n", leftHandSideValue);
                leftHandSideArray[firstCounter] = leftHandSideValue;
            }
            solution = (rightHandSideArray[firstCounter]-leftHandSideArray[firstCounter])/coefficientArray[firstCounter][firstCounter];
            valueArray[firstCounter] = solution;
            sumError = sumError + fabs(oldValueArray[firstCounter]-valueArray[firstCounter]);
            oldValueArray[firstCounter] = valueArray[firstCounter];
        }
        counter++;
    } while (!(sumError<0.0001 || counter>50));
    
    
    return valueArray;
    
}



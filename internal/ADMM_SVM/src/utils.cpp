#include "utils.h"
#include <algorithm>
#include<iostream>

// Function for Machine Epsilon with an 
// initial value provided as EPS. 
double MachineEpsilon()
{   //calculate machine epsilon 
    // with initial value provided as 0.05 
    double EPS = 0.05;
    // taking a floating type variable 
    double prev_epsilon;

    // run until condition satisfy 
    while ((1 + EPS) != 1)
    {
        // copying value of epsilon into previous epsilon 
        prev_epsilon = EPS;

        // dividing epsilon by 2 
        EPS /= 2;
    }
    return prev_epsilon;
}
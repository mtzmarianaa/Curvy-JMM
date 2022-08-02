#include "SoSFunction.h"
#include <math.h>
#include <stdio.h>

double s_function(double x[])
{
    return 1.0;
}

double s_function_twoSections(double x[], int region){
    double SpoS;
    if(region == 1){
        SpoS = 1.0;
    }
    if( region == 2 ){
        SpoS = 1.2;
    }
    else{
        SpoS = 10;
    }
    return SpoS;
}

double s_function_threeSections(double x[], int region){
    double SpoS;
    if(region == 1){
        SpoS = 1.0;
    }
    else if( region == 2 ){
        SpoS = 1.348*2.5; // for the micro wasp: 1.348
    }
    else if(region == 3){
        SpoS = 1.452*2.5; // for the micro wasp: 1.452
    }
    else{
        SpoS = 100;
    }
    return SpoS;
}
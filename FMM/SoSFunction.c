#include "SoSFunction.h"
#include <math.h>

double s_function(double x[])
{
    return 1.0;
}

double s_function_twoSections(double x[], int region){
    double SpoS;
    if(region == 0){
        SpoS = 1.0;
    }
    if( region == 1 ){
        SpoS = 1.2;
    }
    else{
        SpoS = 10;
    }
    return SpoS;
}

double s_function_threeSections(double x[], int region){
    double SpoS;
    if(region == 0){
        SpoS = 1.0;
    }
    if( region == 1 ){
        SpoS = 1.2;
    }
    if(region == 2){
        SpoS = 1.4;
    }
    else{
        SpoS = 10;
    }
    return SpoS;
}
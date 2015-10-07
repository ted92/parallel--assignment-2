/*
 * RaceTrap.java
 *
 * Created on 22. juni 2000, 13:48
 *
 * Brian Vinter
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// include openmp
#include <omp.h>


// max x and max y are the coordinates max of the map
int max_x = 3000;
int max_y = 3000;
int n_cores;

typedef unsigned char byte;

typedef struct {
    double dLength;
    int iSet;
    byte baPath[128];
} cRouteDefinition;


// size of the path
int iSize;

int *iaaXCoors;
int *iaaYCoors;

// table which contains the distances of the set
double **daaDistanceTable;
double dMinRouteLen=10E18;
double dGlobalBest=10E18;


void for_cycle_zRoute(cRouteDefinition * oRoute, byte c, cRouteDefinition ** reply){
    int i;
    double nl;
    cRouteDefinition Perm, Res, *oResult=&Res, *oPermutation=&Perm;
    cRouteDefinition *oBest=*reply;

    #pragma omp for
    for(i=oRoute->iSet;i<iSize;i++){
        // nl is new_length
        nl=oRoute->dLength+daaDistanceTable[oRoute->baPath[oRoute->iSet-1]][oRoute->baPath[i]];
        if(nl<dGlobalBest) {
            // copy the values of the number of bytes from the location pointed to by source directly to the memory block pointed to by destination
            memcpy(oPermutation->baPath, oRoute->baPath, iSize * sizeof(byte));
            oPermutation->baPath[oRoute->iSet] = oPermutation->baPath[i];
            oPermutation->baPath[i] = c;
            oPermutation->dLength = nl;
            oPermutation->iSet = oRoute->iSet + 1;
            zRoute(oPermutation, &oResult);
            if (oResult->dLength < oBest->dLength) { //Best route so far?
                oBest = oResult;
                if (oBest->dLength < dGlobalBest) {
                    dGlobalBest = oBest->dLength;
                }
            }
        }
    }
    return;

}


// counter for parallelize
int count = 0;
void zRoute(cRouteDefinition * oRoute, cRouteDefinition ** reply){
    byte c;
    int i;
    double nl;
    cRouteDefinition Perm, Res, *oResult=&Res, *oPermutation=&Perm;
    cRouteDefinition *oBest=*reply;

    // if it calculated all the possible routes then exit and return the best value
    if(oRoute->iSet==iSize){
        oRoute->dLength+=daaDistanceTable[oRoute->baPath[oRoute->iSet-1]][oRoute->baPath[0]];
        *reply=oRoute;
        return;
    }

    oBest->dLength=dMinRouteLen;

    c=oRoute->baPath[oRoute->iSet];

    if (count == 0) {
        count++;
        #pragma omp parallel num_threads(n_cores)
        {
        // call the parallelised function
            for_cycle_zRoute(oRoute, c, reply);
        }
    }
    // sequential function
    else{
        for(i=oRoute->iSet;i<iSize;i++){
            // nl is new_length
            nl=oRoute->dLength+daaDistanceTable[oRoute->baPath[oRoute->iSet-1]][oRoute->baPath[i]];
            if(nl<dGlobalBest) {
                // copy the values of the number of bytes from the location pointed to by source directly to the memory block pointed to by destination
                memcpy(oPermutation->baPath, oRoute->baPath, iSize * sizeof(byte));
                oPermutation->baPath[oRoute->iSet] = oPermutation->baPath[i];
                oPermutation->baPath[i] = c;
                oPermutation->dLength = nl;
                oPermutation->iSet = oRoute->iSet + 1;
                zRoute(oPermutation, &oResult);
                if (oResult->dLength < oBest->dLength) { //Best route so far?
                    oBest = oResult;
                    if (oBest->dLength < dGlobalBest) {
                        dGlobalBest = oBest->dLength;
                    }
                }
            }
        }
    }
    return;
}

double zEuclidDist(int from, int to){
    double dx=abs(iaaXCoors[from]-iaaXCoors[to]);
    double dy=abs(iaaYCoors[from]-iaaYCoors[to]);

    return(sqrt(dx*dx+dy*dy));
}


// generate random number given a range from 0 to max
// Assumes 0 <= max <= RAND_MAX
// Returns in the half-open interval [0, max]
int random_at_most(int max) {
    unsigned int
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
            num_bins = (unsigned int) max + 1,
            num_rand = (unsigned int) RAND_MAX + 1,
            bin_size = num_rand / num_bins,
            defect   = num_rand % num_bins;

    int x;
    do {
        x = random();
    }
        // This is carefully written not to overflow
    while (num_rand - defect <= (unsigned int)x);

    // Truncated division is intentional
    return x/bin_size;
}


void zReadRoute(/*int size, int max_x, int max_y*/){
    int i,j;

    daaDistanceTable=(double **)malloc(sizeof(double*)*iSize);

    for(i=0;i<iSize;i++)
        daaDistanceTable[i]=(double *)malloc(sizeof(double)*iSize);

    iaaXCoors=(int *)malloc(sizeof(int)*iSize);
    iaaYCoors=(int *)malloc(sizeof(int)*iSize);

    for(i=0;i<iSize;i++){
        // generate random numbers from 0 to x_max and y_max

        iaaXCoors[i] = random_at_most(max_x);
        iaaYCoors[i] = random_at_most(max_y);

        // printf("\nx: %d\ny: %d", iaaXCoors[i], iaaYCoors[i]);

    }

    for(i=0; i<iSize; i++)
        for(j=0; j<iSize; j++)
            daaDistanceTable[i][j]=zEuclidDist(i,j);
}

double second()
{

#include <sys/time.h>

    struct timeval tv;
    struct timezone tz;
    double t;

    gettimeofday(&tv,&tz);

    t= (double)(tv.tv_sec)+(double)(tv.tv_usec/1.0e6);

    return t;
}

main(int argc, char *argv[]) {
    if (argc == 3){
        // convert the argv to integer with stdlib
        iSize = atoi(argv[1]);
        n_cores = atoi(argv[2]);
    }
    else{
        // define a default size, x and y
        iSize = 14;
        n_cores = 4;
    }

    printf("\nthe number of cities is: %d\n the map is: %d x %d\n", iSize, max_x, max_y);


    int i;
    cRouteDefinition *oOriginalRoute;
    cRouteDefinition res, *r=&res;
    double start, stop;

    zReadRoute();

    oOriginalRoute=(cRouteDefinition*)malloc(sizeof(cRouteDefinition));

    for(i=0; i<iSize; i++)
        oOriginalRoute->baPath[i]=(byte)i;

    oOriginalRoute->dLength=0.0;
    oOriginalRoute->iSet=1;

    start=second();
    zRoute(oOriginalRoute,&r);  //Find the best route:)
    stop=second();
    printf("%d cores used.\nRoute length is %lf found in %lf seconds\n", n_cores, dGlobalBest,stop-start);
}







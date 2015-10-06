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

// size is the number of the cities that the traveller has to cross, max x and max y are the coordinates max of the map
int size;
int max_x;
int max_y;

typedef unsigned char byte;

typedef struct {
  double dLength; 
  int iSet;       
  byte baPath[128];
} cRouteDefinition;


int iSize;
int *iaaXCoors; 
int *iaaYCoors; 
double **daaDistanceTable;
double dMinRouteLen=10E18;
double dGlobalBest=10E18;

void zRoute(cRouteDefinition * oRoute, cRouteDefinition ** reply){
  byte c;
  int i;
  double nl;
  cRouteDefinition Perm, Res, *oResult=&Res, *oPermutation=&Perm;
  cRouteDefinition *oBest=*reply;

  if(oRoute->iSet==iSize){
    oRoute->dLength+=daaDistanceTable[oRoute->baPath[oRoute->iSet-1]][oRoute->baPath[0]];
    *reply=oRoute;
    return;
  }

  oBest->dLength=dMinRouteLen;
  
  c=oRoute->baPath[oRoute->iSet];
  for(i=oRoute->iSet;i<iSize;i++){
    nl=oRoute->dLength+daaDistanceTable[oRoute->baPath[oRoute->iSet-1]][oRoute->baPath[i]];
    if(nl<dGlobalBest){
      memcpy(oPermutation->baPath, oRoute->baPath, iSize*sizeof(byte));
      oPermutation->baPath[oRoute->iSet]=oPermutation->baPath[i];
      oPermutation->baPath[i]=c;
      oPermutation->dLength=nl;
      oPermutation->iSet=oRoute->iSet+1;
      zRoute(oPermutation, &oResult);
      if(oResult->dLength<oBest->dLength){ //Best route so far?
        oBest=oResult;
        if(oBest->dLength<dGlobalBest){
          dGlobalBest=oBest->dLength;
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
long random_at_most(long max) {
  unsigned long
  // max <= RAND_MAX < ULONG_MAX, so this is okay.
          num_bins = (unsigned long) max + 1,
          num_rand = (unsigned long) RAND_MAX + 1,
          bin_size = num_rand / num_bins,
          defect   = num_rand % num_bins;

  long x;
  do {
    x = random();
  }
    // This is carefully written not to overflow
  while (num_rand - defect <= (unsigned long)x);

  // Truncated division is intentional
  return x/bin_size;
}


void zReadRoute(/*int size, int max_x, int max_y*/){
  int i,j;

  iSize = size;
  daaDistanceTable=(double **)malloc(sizeof(double*)*iSize);

  for(i=0;i<iSize;i++)
    daaDistanceTable[i]=(double *)malloc(sizeof(double)*iSize);
  
  iaaXCoors=(int *)malloc(sizeof(int)*iSize);
  iaaYCoors=(int *)malloc(sizeof(int)*iSize);
  
  for(i=0;i<iSize;i++){
    // generate random numbers from 0 to x_max and y_max

    scanf("%d %d",&(iaaXCoors[i]), &(iaaYCoors[i]));
    //printf("%d %d\n", iaaXCoors[i], iaaYCoors[i]);
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

  if (argc == 4){
    // convert the argv to integer with stdlib
    size = atoi(argv[1]);
    max_x = atoi(argv[2]);
    max_y = atoi(argv[3]);
    printf("\n---3 arguments given---\n the size is: %d\n the max y value is: %d\n the max x value is: %d\n", size, max_x, max_y);

    printf("\n---random number x---\n%ld", random_at_most(max_x));
    printf("\n---random number y---\n%ld", random_at_most(max_y));

  }
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
  printf("Route length is %lf found in %lf seconds\n",dGlobalBest,stop-start);
}







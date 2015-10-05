/*
 * RaceTrap.java
 *
 * Created on 22. juni 2000, 13:48
 * 
 * Brian Vinter
 */

#include <math.h>
#include <stdlib.h>

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
    
void zReadRoute(){
  int i,j;
    
  scanf("%d",&iSize);

  daaDistanceTable=(double **)malloc(sizeof(double*)*iSize);
  for(i=0;i<iSize;i++)
    daaDistanceTable[i]=(double *)malloc(sizeof(double)*iSize);
  
  iaaXCoors=(int *)malloc(sizeof(int)*iSize);
  iaaYCoors=(int *)malloc(sizeof(int)*iSize);
  
  for(i=0;i<iSize;i++){
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

main() {
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







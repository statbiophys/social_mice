#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Variables */
  int n, s,s1,s2,n1,n2,e, yInd, Tind, t, n3,
    nNodes, nEdges, nTriplets,
    maxState, 
   *y,
    *edgeEnds, *nStates,
    *tripletEnds;

  //sizeEdgeBel[3], sizeLogZ[2],

  mwSize sizeEdgeBel[3], sizeLogZ[2];
   
   double pot,Z,
     *nodePot, *edgePot, *nodeBel, *edgeBel, *logZ,
     *tripletPot, *tripletBel;
   
   /* Input */
   
   nodePot = mxGetPr(prhs[0]);
   edgePot = mxGetPr(prhs[1]);
   edgeEnds = (int*)mxGetPr(prhs[2]);
   tripletPot = mxGetPr(prhs[3]);
   tripletEnds = (int*)mxGetPr(prhs[4]);
   nStates = (int*)mxGetPr(prhs[5]);
   
   
   if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32"))
        mexErrMsgTxt("edgeEnds, tripletEnds, and nStates must be int32");
   
   /* Compute Sizes */
   
   nNodes = mxGetDimensions(prhs[0])[0];
   maxState = mxGetDimensions(prhs[0])[1];
   nEdges = mxGetDimensions(prhs[2])[0];
   nTriplets = mxGetDimensions(prhs[4])[0];
   
   /* Output */
   sizeEdgeBel[0] = maxState;
   sizeEdgeBel[1] = maxState;
   sizeEdgeBel[2] = nEdges;
   sizeLogZ[0] = 1;
   sizeLogZ[1] = 1;
   plhs[0] = mxCreateNumericArray(2,mxGetDimensions(prhs[0]),mxDOUBLE_CLASS,mxREAL);
   plhs[1] = mxCreateNumericArray(3,sizeEdgeBel,mxDOUBLE_CLASS,mxREAL);
   plhs[2] = mxCreateNumericArray(2,sizeLogZ,mxDOUBLE_CLASS,mxREAL);
   nodeBel = mxGetPr(plhs[0]);
   edgeBel = mxGetPr(plhs[1]);
   logZ = mxGetPr(plhs[2]);
   
   /* Initialize */
   y = mxCalloc(nNodes,sizeof(int));
   Z = 0;
   
   while(1)
   {
      pot = 1;
      
   /* Node */
      for(n = 0; n < nNodes; n++)
      {
         pot *= nodePot[n + nNodes*y[n]];
      }
      
   /* Edges */
      for(e = 0; e < nEdges; e++)
      {
         n1 = edgeEnds[e]-1;
         n2 = edgeEnds[e+nEdges]-1;
         pot *= edgePot[y[n1] + maxState*(y[n2] + maxState*e)];
      }

   /* Triplets */
      for(t = 0; t < nTriplets; t++)
      {
	n1 = tripletEnds[t]-1;
	n2 = tripletEnds[t+nTriplets]-1;
	n3 = tripletEnds[t+nTriplets*2]-1;
	pot *= tripletPot[y[n1] + maxState*(y[n2] + maxState*(y[n3] + maxState*t))];
      }
      

      
   /* Update Z */
      Z += pot;
      
   /* Go to next y */
      
      for(yInd = 0; yInd < nNodes; yInd++)
      {
         y[yInd] += 1;
         if(y[yInd] < nStates[yInd])
         {
          break;  
         }
         else
         {
            y[yInd] = 0;
         }
      }

   /* Stop when we are done all y combinations */
      if(yInd == nNodes)
      {
         break;
      }
   }
   

   *logZ = log(Z);
   
   /* Free memory */
   mxFree(y);
}

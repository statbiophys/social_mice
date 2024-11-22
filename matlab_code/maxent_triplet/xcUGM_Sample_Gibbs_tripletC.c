#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variables */
    int i,n,s,e,n1,n2,Vind,Tind,maxIter,burnIn,
      t,n3,
      nNodes, nEdges, nTriplets, maxState,
    *edgeEnds, *nStates, *V, *E,*y,*S,
    *tripletEnds, *tripleV, *tripleE;

    mwSize dims[2];
    
    double *pot,z,U,u,
    *nodePot, *edgePot, *tripletPot;
    
   /* Input */
    
    nodePot = mxGetPr(prhs[0]);
    edgePot = mxGetPr(prhs[1]);
    tripletPot = mxGetPr(prhs[2]);
    edgeEnds = (int*)mxGetPr(prhs[3]);
    nStates = (int*)mxGetPr(prhs[4]);
    V = (int*)mxGetPr(prhs[5]);
    E = (int*)mxGetPr(prhs[6]);
    tripletEnds = (int*)mxGetPr(prhs[7]);
    tripleV = (int*)mxGetPr(prhs[8]);
    tripleE = (int*)mxGetPr(prhs[9]);
    maxIter = ((int*)mxGetPr(prhs[10]))[0];
    burnIn = ((int*)mxGetPr(prhs[11]))[0];
    y = (int*)mxGetPr(prhs[12]);
		
	if (!mxIsClass(prhs[3],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32")||!mxIsClass(prhs[7],"int32")||!mxIsClass(prhs[8],"int32")||!mxIsClass(prhs[9],"int32")||!mxIsClass(prhs[10],"int32")||!mxIsClass(prhs[11],"int32")||!mxIsClass(prhs[12],"int32"))
        mexErrMsgTxt("edgeEnds, nStates, V, E, tripletEnds, tripleV, tripleE,  maxIter, burnIn, y must be int32");
    
   /* Compute Sizes */
    
    nNodes = mxGetDimensions(prhs[0])[0];
    maxState = mxGetDimensions(prhs[0])[1];
    nEdges = mxGetDimensions(prhs[2])[0];
    nTriplets = mxGetDimensions(prhs[3])[0];

   /* Output */
    pot = mxCalloc(maxState,sizeof(double));
    
    dims[0] = nNodes;
    dims[1] = maxIter;

    // /* Check dims size */
    //size_t temp = sizeof(dims)/sizeof(dims[0]);
    //printf("size of dims = %d\n", temp);
   
    
    plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    
    S = mxGetData(plhs[0]);

 
    
   /* Initialize to States with highest node potentials*/
   /* for(n = 0; n < nNodes; n++)
    {
        u = -1;
        U = 0;
        for(s = 0; s < nStates[n]; s++)
        {
            if(nodePot[n+nNodes*s] > u)
            {
                u = nodePot[n+nNodes*s];
                U = s;
            }
        }
        y[n] = U;
    }
    */
	
	for(n = 0; n < nNodes; n++)
		y[n]--;
    
    for(i = 0; i < burnIn+maxIter; i++)
    {

        for(n = 0; n < nNodes; n++)
        {
            /* Compute Node Potential */
            for(s = 0; s < nStates[n]; s++)
            {
                pot[s] = nodePot[n + nNodes*s];
            }
            
           /* Iterate over Neighbors */
            for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
            {
                e = E[Vind]-1;
                n1 = edgeEnds[e]-1;
                n2 = edgeEnds[e+nEdges]-1;
                 
                /* Multiply Edge Potentials */
                if(n == n1)
                {
                   for(s = 0; s < nStates[n]; s++)
                   {
                        pot[s] *= edgePot[s+maxState*(y[n2] + maxState*e)];
                   }
                    
                }
                else
                {
                    for(s = 0; s < nStates[n]; s++)
                    {
                        pot[s] *= edgePot[y[n1]+maxState*(s + maxState*e)];
                    }
                }  
            }

	    /* Iterate over Triplets */
	    for(Tind = tripleV[n]-1; Tind < tripleV[n+1]-1; Tind++)
	      {
		t = tripleE[Tind]-1;
		n1 = tripletEnds[t]-1;
		n2 = tripletEnds[t+nTriplets]-1;
		n3 = tripletEnds[t+nTriplets*2]-1;

		/* Multiply Triplet Potentials */
		if(t==n1)
		  {
		    for(s = 0; s < nStates[t]; s++)
		      {
			pot[s] *= tripletPot[s+ maxState*(y[n2] + maxState*(y[n3] + maxState*t))];
		      }
		  }
		else if(t==n2)
		  {
		    for(s = 0; s < nStates[n]; s++)
		      {
			pot[s] *= tripletPot[y[n1] + maxState*(s + maxState*(y[n3] + maxState*t))];
		      }
		  }
		else
		  {
		    for(s = 0; s < nStates[n]; s++)
		      {
			pot[s] *= tripletPot[y[n1] + maxState*(y[n2] + maxState*(s + maxState * t))];
		      }
		  }
	      }
            
            /* Normalize */
            z = 0;
            for(s = 0; s < nStates[n]; s++)
                z = z + pot[s];
            for(s = 0; s < nStates[n]; s++)
                pot[s] /= z;
            
            /* Display */
            for(s = 0; s < nStates[n]; s++)
            
            /* Sample Discrete State */
            U = rand()/((double)RAND_MAX + 1);
            u = 0;
            for(s = 0; s < nStates[n]; s++)
            {
                u += pot[s];
                if(u > U)
                {
                    break;
                }
            }
            y[n] = s;
        }
        
        if(i >= burnIn)
        {        /* Record Sample */
            for(n = 0; n < nNodes; n++)
            {
                S[n + nNodes*(i-burnIn)] = y[n]+1;
            }
        }
    }
    
   /* Free memory */
    mxFree(pot);
}

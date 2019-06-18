#include "mex.h"

/*
 * LocalQ.c - LEMming MEX-file
 *
 * Loops through the entire topography and extracts the local fluxes of
 * regolith into each cell
 *
 * This is a MEX-file for MATLAB.
 * 
 */
 
/* $Revision: 0.0.2 for LEMming v.016 $ */

void LocalQ(double Qin[], int x, int y, double Qout[], double Ts[], mwIndex irs[], mwIndex jcs[])
{

  mwIndex xi,yi,from,erg,to,subs[2], stCol, endCol, rowi;
  double w1,w1p,w0,w0p,wsum;
  const double sqrt2 = 1.4142136;
  int diri;

/*  Calculate the linear offsets to each neighbor pixel and store then in the array "neighbors" */
/* 	       		  = { N, S, E, W,  NW  ,  SW  , NE  , SE  } */
  int neighbor[8] = {-1, 1, y,-y,-(y+1),-(y-1),(y-1),(y+1)};


          for (xi = 1; xi < x-1; xi++) {
            for (yi = 1; yi < y-1; yi++) {
            
				
				from = yi+y*(xi-1); /*  Get the linear subscript */
				
				if (Qout[from] > 0)  {  /*  only do this for cells with an outflux */


					  stCol = (jcs[from]); /*  index in Ts and irs of the first "to" pixel in the "from" column */
					  endCol = (jcs[from+1])-1; /*  index of the last "to" pixel in the "from" column */

					if (stCol <= endCol) {	/*  If false, there is a defined outflux at Qout[from] but no corresonding entry in the "from" column of the sparse flow matrix. */

						/* This loop extracts the row indices and values from the sparse matrix,
						determines their proper positions in the output flux matrix, and adds
						them to the appropriate cell there */
						
						wsum = 0.0;		/* zero out weight sum */
							
						/* These first two half-loops sum the distance-weighted flux weights */	

 						for (rowi = stCol; rowi <= endCol; rowi++) {
							
							for (diri = 0; diri <= 3; diri++) {				 /* First do straight N,S,E,W */		
							   
								if (irs[rowi] == (from + neighbor[diri])) {  /* is this the direction? */
								
									wsum += -Ts[rowi];
								
								} /*  if irs[] */
							
							
							} /*  for diri */
																
							for (diri = 4; diri <= 7; diri++) {				 /* Now do diagonal NW,SW,NE,SE */   
							   
								if (irs[rowi] == (from + neighbor[diri])) {  /* is this the direction? */
								
									wsum += -Ts[rowi]/sqrt2;
								
								} /*  if irs[] */
							
							} /*  for diri */

 						} /*  for rowi */
 						
 						
					    /* The second pair of half-loops uses the sum to calculate the final weights and fluxes */

 						for (rowi = stCol; rowi <= endCol; rowi++) {
 						
							for (diri = 0; diri <= 3; diri++) {				 /* First do straight N,S,E,W */		
							   
								if (irs[rowi] == (from + neighbor[diri])) {  /* is this the direction? */
								
									w1 = -Ts[rowi];
									w1p = w1/wsum;			/* Calc distance-weighted flux fraction */
									
									Qin[irs[rowi]] += (w1p * Qout[from]);
								
								} /*  if irs[] */
								
							} /*  for diri */

							for (diri = 4; diri <= 7; diri++) {				 /* Now do diagonal NW,SW,NE,SE */
							   
								if (irs[rowi] == (from + neighbor[diri])) {  /* is this the direction? */
								
									w0 = -Ts[rowi];
									w0p = (w0/sqrt2)/wsum; /* Calc diagonal distance-weighted flux fraction */
								
									Qin[irs[rowi]] += (w0p * Qout[from]);
								
								} /*  if irs[] */
							
							} /*  for diri */

 						} /*  for rowi */
 						
            
					} /*  if StCol */
               } /*  if Qout */
               
            }  /*  for yi */
        } /*  for xi */

}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *Qin,*Qout,*Ts;
  int x,y;
  mwIndex *irs,*jcs;
  
  /* Check for proper number of arguments. */
  if(nrhs!=4) {
    mexErrMsgTxt("Four inputs required.");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments.");
  }


  x = mxGetScalar(prhs[0]);
  y = mxGetScalar(prhs[1]);

  /* Create matrix for the return argument. */
  
   plhs[0] = mxCreateDoubleMatrix(y,x, mxREAL);
 
  
  /* Assign pointers to each input and output. */

  Qout = mxGetPr(prhs[2]); /*  Matrix of fluxes out of each cell */

  Ts = mxGetPr(prhs[3]); /*  Sparse matrix of relative flow from each pixel to each other pixel */
  irs = mxGetIr(prhs[3]); /*  Rows of nonzero elements in Ts */
  jcs = mxGetJc(prhs[3]); /*  Columns of nonzero elements in Ts */
 
  Qin = mxGetPr(plhs[0]);
  
  
  /* Call the LocalQ subroutine. */
	LocalQ(Qin,x,y,Qout,Ts,irs,jcs);

}

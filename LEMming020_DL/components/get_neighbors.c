#include "mex.h"

/*
 * get_neighbors.c - LEMming MEX-file
 *
 *         Find all neighbors (w,z) of (r,c) for which R(w,z) is not NaN, and
 *         for which E(w,z) = E(r,c).  Count how many there are and record
 *         their locations.  Pixel flow will be distributed to all such
 *         neighbors equally. Replaces a slow loop in the plateau_flow_weights()
 *		   component of the flow_matrix() function of the "upslope" package
 *
 * This is a MEX-file for MATLAB.
 * 
 */
 
/* $Revision: 0.0.1 for LEMming v.016 $ */

void get_neighbors(double *x, double *y, double *row, double *col, 
				   double topo[], double R[], 
				   double ww[], double zz[], double *neighbor_count)
{

int w,wstart,wend, z,zstart,zend, nei, src, r,c, nc=0;	/* intialize counters */

/* convert to zero-based indexing */
r = (int)*row-1;
c = (int)*col-1;

/* linear index of r,c */
src = r + ((*y) * (c)); /* Get the linear subscript */

/* establish bounds for loop */
{ if (r-1 > 0)
	wstart = r-1;
else
	wstart = 0;
}
{ if (*y-1 < r+1)
	wend = *y-1;
else
	wend = r+1;
}

{ if (c-1 > 0)
	zstart = c-1;
else
	zstart = 0;
}
{ if (*x-1 < c+1)
	zend = *x-1;
else
	zend = c+1;
}

/* Now get down to business:
 *         Find all neighbors (w,z) of (row,col) for which R(w,z) is not NaN, and
 *         for which topo(w,z) = topo(r,c).  Count how many there are and record
 *         their locations. */
 
    for (w = wstart; w <= wend; w++) {
	    for (z = zstart; z <= zend; z++) {
	    
			nei = (w) + ((*y) * (z)); /* Get the linear subscript of w,z */
	    
            if ( !mxIsNaN(R[nei]) && topo[src] == topo[nei] ) {
            
				ww[nc] = w+1; /* loop indices are zero-based */
				zz[nc] = z+1;
				
				nc = nc++; /* increment neighbor counter */
				
			} /* if mxIsNaN...*/
		} /* for z */
    } /* for w */

*neighbor_count = nc;

} /* end function get_neighbors */


/* Gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *topo,*R, *ww, *zz, *neighbor_count;
  double *row,*col,*x,*y;
  
  /* Check for proper number of arguments. */
  if(nrhs!=6) {
    mexErrMsgTxt("Six inputs required.");
  } else if(nlhs>3) {
    mexErrMsgTxt("Too many output arguments.");
  }


  /* Create matrices for the return arguments. */
  
   plhs[0] = mxCreateDoubleMatrix(8,1, mxREAL); /* ww */
   plhs[1] = mxCreateDoubleMatrix(8,1, mxREAL); /* zz */
   plhs[2] = mxCreateDoubleScalar(0); /* neighbor_count */
   
   
  /* Assign pointers to each input and output. */

  /* inputs */
  x = mxGetPr(prhs[0]); /* number of columns in R */
  y = mxGetPr(prhs[1]); /* number of rows in R */
  row = mxGetPr(prhs[2]); /* row of interest */
  col = mxGetPr(prhs[3]); /* column of interest */
  topo = mxGetPr(prhs[4]); /* Topographic elevation matrix */
  R = mxGetPr(prhs[5]); /* Flow direction matrix */

  /* outputs */
  ww = mxGetPr(plhs[0]);
  zz = mxGetPr(plhs[1]);
  neighbor_count = mxGetPr(plhs[2]);
  
  /* Call the get_neighbors subroutine. */
	get_neighbors(x,y,row,col,topo,R,ww,zz,neighbor_count);

} /* end gateway function */

Matlab-callable C (mex-file) functions to generate predictions for the models described in Smith, P. L., "Linking the Diffusion Model and General Recognition Theory: Circular Diffusion with Bivariate-Normally Distributed Drift Rates." 

I. Resources (Matlab)
=====================

Circular diffusion model with across-trial variability in drift rate

    1. vdcircle600.c and vdcircle600.m: Predictions for the independent-components circular diffusion model.
    2. grtrot600.c and grtrot600.m: Predictions for the rotationally-invariant model.
    3. grtgen600.c and grtgen600.m: Predictions for the general model with arbitrary correlation, rho.

II. Software requirements
=========================
The function were developed and tested using the Gnu C-compiler, gcc, under Linux. They make use of the Gnu scientific library, gsl, and its associated library, gslcblas, to compute things like Bessel functions. These functions have Gnu-specific names like gsl_sf_bessel_J1, etc. All of these functions should be available in other math libraries, but if you wish to use something other than gsl you will need to track these function calls down and change them. Gcc, gsl, and clblas are open-source software and should be freely downloadable. Matlab will complain if you use a version of gcc that is later than the one it currently supports, but it doesn't appear to cause execution problems.

The functions are compiled from the Matlab command window using commands of the form:

    >>  mex vdcircle600.c -lgsl -lgslcblas  -lm

and so on. The discretization of time and phase angle are set internally by constants

    #define nw 50    /* Number of steps on circle */
    #define sz 600   /* Number of time steps */

at the top of the functions. Currently time is discretized on the range [0, tmax] in 600 steps (thus the 600 designation) and phase angle is discretized in 50 steps. You can change these values to anything else you'd like and everything should continue to work.

III. Calling conventions
========================

1. The independent components model
-----------------------------------

The m-file associated with a given C program is a Matlab-compliant header file that echoes the calling convention for the function. So, for example, 

    >> help vdcircle600

will give you:

     ==========================================================================
       Circular diffusion model, indepedendent Gaussian drift rates. 
     
       [T, Gt, Theta, Ptheta, Mt] = vdcircle600(P, tmax, badix);
        P = [v1, v2, eta1, eta2, sigma, a]
 
       Building:
                 mex vdcircle600.c -lgsl -lgslcblas  -lm 
      ===========================================================================

1. Input parameters

The input parameters, P, are the horizontal and vertical components of the mean drift rate (v1 and v2), the horizontal and vertical drift rate standard deviations (eta1, eta2), the infinitesimal standard deviation of the diffusion process (sigma), and the radius of the criterion circle (a). As discussed in the paper, when the stimulus is in canonical orientation, such that v1 = ||v||, v2 = 0, then eta1 and eta2 can be interpreted as the radial and tangential components of drift rate variability.

Predictions are generated on the range [0, tmax] in 600 steps. 

The infinite series representation for the first-passage time density of the Bessel process, dhamana (from Hamana & Matsumoto, 2013) is a "large time" representation, which can be inaccurate at small values of time when the series is truncated. Currently, the series is truncated at 50 steps. This can sometimes lead to a small negative excursion in the density function at small times. In practice this doesn't create problems, because any such excursions are confined to the first few time steps and don't affect the part of the function that's needed for predictions. The parameter "badix" (for bad index) zeros out this part of the density. Setting badix = 5 usually works for most values of the parameters of interest. If there is a problem, it will be obvious graphically. 

2. Output parameters:

These are:

    T: 600 element time vector 
    Gt: 50 x 600 element array of joint density values
    Theta: 50 element vector of phase angles
    Ptheta: 50 element vector of decision outcome densities
    Mt: 50 element vector of mean decision times

A typical function call would look like:

    >> [t,gt,th,pth,mth]=vdcircle600([1.5,1.5,0.5,0.5,1.0,1.0], 2.0, 5);

and would generate predictions on the range [0, 2.0]s

The Matlab command

    >> plot(t, gt)

will plot the joint densities and

    >> plot(th, pth)

and

    >> plot(th, mth)

will plot the distributions of decision outcomes and mean decision times, respectively. 


2. The Rotationally Invariant Model
-----------------------------------

    >> help grtrot600                                               
      ==========================================================================
      Circular diffusion model. Rotational invariant correlated drift rates.
 
       [T, Gt, Theta, Ptheta, Mt] = grtrot600(P, tmax, badix);
        P = [v1, v2, eta1, eta2, sigma, a]
 
       Building:
           mex grtrot600.c -lgsl -lgslcblas  -lm 
      ===========================================================================

The rotationally-invariant model is called in exactly the same way as is the independent components model. No correlation need be given because the program computes it from the mean drift vector as rho = arctan(v2/v2). The values of eta1 and eta2 are interpreted as the radial and tangential components of drift rate standard deviation. The program computes new values of standard deviation in rotated coordinates using Equations (26)-(29). There are some active trace statements in the program that keep track of the transformed values of drift rate mean and standard deviations. If you want to suppress these for fitting purposes you will need to comment out the associated mexPrint calls.

    >> [t1,g1,th1,pth1,mt1]=grtrot600([1.5,0,1.0,0.5,1.0,1.0], 2.0, 5);
    >> [t2,g2,th2,pth2,mt2]=grtrot600([1.0607,1.0607,1.0,0.5,1.0,1.0], 2.0, 5);

generates predictions that are pi/4 out of phase with one another, but are otherwise the same. (There are some small numerical discrepancies in the means at the extreme tails where the probability of a response goes to zero.)

2. The General Model
--------------------

    >> help grtgen600
      ==========================================================================
       Circular diffusion model. Rotational invariant correlated drift rates.
       Correlation in canonical orientation.
       [T, Gt, Theta, Ptheta, Mt] = grtgen600(P, tmax, badix);
        P = [v1, v2, eta1, eta2, sigma, a, rho]
 
       Building:
                 mex grtgen600.c -lgsl -lgslcblas  -lm 
      ===========================================================================

The parameter rho is interpreted as the correlation between the components of drift rate for a stimulus in canonical orientation. Likewise, eta1 and eta2 are interpreted as the drift rate standard deviations in canonical orientation. The rotation angle is computed as phi = arcsin(v2 / v2). Transformed values of rho, eta1, and eta2 are then computed for the drift vector relative to the new ||v|| = (v1, v2) coordinates. These values are currently output as trace statement but may be suppressed by commenting out the associated mexPrint statements. 

All of the usual "as is" caveats apply to the use of this software!



 









/* ==========================================================================
%  Circular diffusion model, rho-correlated Gaussian drift rates. 
%  Correlation in canonical orientation.
%
%  [T, Gt, Theta, Ptheta, Mt] = grtgen600(P, tmax, badix);
%   P = [v1, v2, eta1, eta2, sigma, a]
%
%  Building:
%            mex grtgen600.c -lgsl -lgslcblas  -lm 
% ===========================================================================
*/

#include <mex.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>

#define kmax 50  /* Maximum number of eigenvalues in dhamana */
#define nw 120    /* Number of steps on circle */
#define sz 600   /* Number of time steps */
#define NP 7     /* Number of input parameters */

const double pi = 3.141592653589793;

void dhamana(double *T, double *Gt0, double *P0, double h, int badix) {
    /* 
      ---------------------------------------------------------------
       First-passage-time density for Bessel process.
       Computes roots of J0 using Gnu GSL library.
      ----------------------------------------------------------------
    */
    double J0k[kmax], J0k_squared[kmax], J1k[kmax]; 
    double a, a2, sigma, sigma2, scaler; 
    int i, k;

    a = P0[0];
    sigma = P0[1];
    sigma2 = sigma * sigma;
    a2 = a * a;
    scaler = sigma2 / a2;
    
    /* Roots of J0k */
    for (k = 0; k < kmax; k++) { 
        J0k[k] = gsl_sf_bessel_zero_J0(k+1);
        /* mexPrintf("k = %6d  J0k[k] = %6.4f\n", k, J0k[k]); - OK */
    }
 
    /* Evaluate Bessel function at the roots */
    for (k = 0; k < kmax; k++) {
        J0k_squared[k] = J0k[k] * J0k[k];
        /* J1k[k] = j1(J0k[k]); */  /* besselj in Gnu library */
        J1k[k] = gsl_sf_bessel_J1(J0k[k]); /* GSL library */
    }
    T[0] = 0;
    Gt0[0] = 0;
    for (i = 1; i < sz; i++) {    
        T[i] = i * h;
        Gt0[i] = 0;
        for (k = 0; k < kmax; k++) {
            Gt0[i] += J0k[k] * exp(-J0k_squared[k] * sigma2 * T[i] / (2.0 * a2)) / J1k[k]; 
        }
        Gt0[i] *= scaler;
        if (i <= badix) {
            Gt0[i] = 0;
        }
    }
};

void grtgen600(double *T, double *Gt, double *Theta, double *Ptheta, double *Mt, 
              double *P, double tmax, int badix) {
    /* -----------------------------------------------------------------------------------------------
       Calculate first-passage-time density and response probabilities for circular diffusion process
      ------------------------------------------------------------------------------------------------ */
    
    double Gt0[sz], P0[2];
    double w, two_pi, h, v1, v2, eta1, eta2, sigma, a, sigma2, mt, a1, a2, b,
           rho, kappa, A, B, C, D, E, F, G, H, v1_2, v2_2, eta1_2, eta2_2, GH_minus_rho2_on_kappa,
           exponent_zt, zt, totalmass, covxy, varx_vary, phi, eta1_minus_eta2, eta1_plus_eta2,
           sin_2_phi, cos_2_phi, varx, vary, rhox, old_rho;

    int i, k;

    two_pi = 2.0 * pi;
    w = 2.0 * pi / nw;
  
    /* Parameters */
    h = tmax / sz; 
    v1 = P[0];
    v2 = P[1];
    eta1 = P[2];
    eta2 = P[3];
    if (eta1 <1e-1) {
       eta1 = 0.01;
    }
    if (eta2 <1e-1) {
       eta2 = 0.01;
    }     
    sigma = P[4];
    a = P[5];
    old_rho = P[6];

    /*mexPrintf("w= %6.4f h = %6.4f \n", w, h);    */
    sigma2 = sigma * sigma; 
    eta1_2 = eta1 * eta1;
    eta2_2 = eta2 * eta2;
    v1_2 = v1 * v1;
    v2_2 = v2 * v2;
   /* Correlation depends on rotation of drift vectors. */
    phi = atan(v2 / v1);  /* Rotation angle */
    eta1_minus_eta2 = eta1 - eta2;
    eta1_plus_eta2 = eta1 + eta2;
    sin_2_phi = sin(2.0 * phi);
    cos_2_phi = cos(2.0 * phi);
    covxy = 0.5 * sin_2_phi * eta1_minus_eta2 * eta1_plus_eta2 + cos_2_phi * old_rho * eta1 * eta2;
    varx = cos(phi) * cos(phi) * eta1_2 + sin(phi) * sin(phi) * eta2_2 
           - sin_2_phi * old_rho * eta1 * eta2;
    vary = cos(phi) * cos(phi) * eta2_2 + sin(phi) * sin(phi) * eta1_2
           + sin_2_phi * old_rho * eta1 * eta2;

    varx_vary = varx * vary;
    rho = covxy / sqrt(varx_vary); 
 
    kappa = sqrt(1.0 - rho * rho);
    mexPrintf("v1= %6.4f v2 = %6.4f eta1 = %6.4f eta2 = %6.4f varx = %6.4f vary = %6.4f\n", 
              v1, v2, eta1, eta2, varx, vary); 
    mexPrintf("phi= %6.4f covxy = %6.4f rho = %6.4f varx*vary = %6.3f kappa = %6.4f \n", 
              phi, covxy, rho, varx_vary, kappa); 

    /* Rescale variances */
    eta1 = sqrt(varx);
    eta2 = sqrt(vary);
    eta1_2 = eta1 * eta1;
    eta2_2 = eta2 * eta2;
    mexPrintf("New eta1= %6.4f eta2 = %6.4f\n", eta1, eta2); 

    /* Density for zero-drift process */
    P0[0] = a;
    P0[1] = sigma;
    dhamana(T, Gt0, P0, h, badix); 
    
    /* Response circle (1 x nw) */
    Theta[0] = -pi;
    for (i = 1; i < nw; i++) {
        Theta[i] = Theta[i-1] + w;
    }

   /* Joint RT distribution (nw * sz) - make Matlab conformant */
   for (k = 0; k < sz; k++) {
       for (i = 0; i < nw; i++) {
            a1 = a * cos(Theta[i]) / sigma2;
            a2 = a * sin(Theta[i]) / sigma2;
            b = T[k] / (2.0 * sigma2);
            A = b * eta1_2;
            B = (a1 - 2.0 * b * v1) * eta1;
            C = a1 * v1 - b * v1_2;
            D = b * eta2_2;
            E = (a2 - 2.0 * b * v2) * eta2;
            F = a2 * v2 - b * v2_2;
            G = 1.0 + 2.0 * kappa * A;
            H = 1.0 + 2.0 * kappa * D;
            GH_minus_rho2_on_kappa = (1.0 + 2.0 * (A + D)) + 4.0 * kappa * A * D;
            exponent_zt = (G * E * E + 2.0 * rho * B * E + B * B * H) / (2.0 * GH_minus_rho2_on_kappa);
            zt = exp(C + F) * exp(exponent_zt) / sqrt(GH_minus_rho2_on_kappa);
            Gt[nw * k + i] =  zt * Gt0[k] / two_pi;  
        }
    } 
    /* Total mass */
    totalmass = 0;
    for (i = 0; i < nw; i++) {
       for (k = 1; k < sz; k++) {
           totalmass += (Gt[nw * k + i] + Gt[nw * (k - 1) + i]) / 2.0;
       } 
    }
    totalmass *= w * h;
    mexPrintf("totalmass = %6.4f\n", totalmass);  
    /* Integrate joint densities to get means hitting probabilities */
    for (i = 0; i < nw; i++) {
       Ptheta[i] = 0;
       Mt[i] = 0;
       for (k = 1; k < sz; k++) {
           Ptheta[i] += (Gt[nw * k + i] + Gt[nw * (k - 1) + i]) /2.0;
           Mt[i] += (T[k] * Gt[nw * k + i] + T[k - 1] * Gt[nw * (k - 1)+ i]) / 2.0; 
       }
       Ptheta[i] *= h / totalmass;
       Mt[i] *= h / Ptheta[i] / totalmass;  
   }

} /* vgrtgen600 */
   
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
 /*
     =======================================================================
     Matlab gateway routine.
     =======================================================================
 */
 

int badix; 

double *T, *Gt, *Theta, *Ptheta, *Mt, *P;
 
double tmax, badi;
 
unsigned n, m;

    if (nrhs != 3) {
         mexErrMsgTxt("grtgen600: Requires 3 input args.");
    } else if (nlhs != 5) {
        mexErrMsgTxt("grtgen600: Requires 5 output args."); 
    }

    /*
      -----------------------------------------------------------------------
      Check all input argument dimensions.
      -----------------------------------------------------------------------   
    */

    /* P */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    if (!mxIsDouble(prhs[0]) || !(m * n == NP)) {
        mexPrintf("P is %4d x %4d \n", m, n);
        mexErrMsgTxt("grtgen600: Wrong size P");
    } else {
        P = mxGetPr(prhs[0]);
    }
    /* tmax */
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    if (!mxIsDouble(prhs[1]) || !(m * n == 1)) {
        mexErrMsgTxt("grtgen600: tmax must be a scalar");
    } else { 
        tmax = mxGetScalar(prhs[1]);
    }
    if (tmax <= 0.0) {
        mexPrintf("tmax =  %6.2f \n", tmax);
        mexErrMsgTxt("tmax must be positive");
    } 

    /* badi */
    m = mxGetM(prhs[2]);
    n = mxGetN(prhs[2]);
    if (!mxIsDouble(prhs[2]) || !(m * n == 1)) {
        mexErrMsgTxt("grtgen600: badi must be a scalar");
    } else {
        badi = mxGetScalar(prhs[2]); 
        badix = (int)(badi+0.5); 
    }  
 
    /*
      -----------------------------------------------------------------------
      Create output arrays.
      -----------------------------------------------------------------------
    */
 
    /* T */
    plhs[0] = mxCreateDoubleMatrix(1, sz, mxREAL);
    T = mxGetPr(plhs[0]);
    
    /* Gt */
    plhs[1] = mxCreateDoubleMatrix(nw, sz, mxREAL);
    Gt = mxGetPr(plhs[1]);
    
     /* Theta */
    plhs[2] = mxCreateDoubleMatrix(1, nw, mxREAL);
    Theta = mxGetPr(plhs[2]);


    /* Ptheta */
    plhs[3] = mxCreateDoubleMatrix(1, nw, mxREAL);
    Ptheta = mxGetPr(plhs[3]);

    /* Mt */
    plhs[4] = mxCreateDoubleMatrix(1, nw, mxREAL);
    Mt = mxGetPr(plhs[4]);


    /*
      -----------------------------------------------------------------------
      Run the C-function.
      -----------------------------------------------------------------------
    */

    grtgen600(T, Gt, Theta, Ptheta, Mt, P, tmax, badix);
}



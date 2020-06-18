// HERE we have the THEORY for N choices, alpha, beta >=1
// only theory, so it is very FAST

#include <stdio.h> 
#include <stdlib.h>
//#include <math.h>
#include "randomc.h"                   // define classes for random number generators
#include "userintf.cpp"                // define system specific user interface
#include "mersenne.cpp"                // members of class TRandomMersenne


int main(void)
{
	//Parameters
	const int     N_vec[22] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 100000, 500000};
	const int     num_N = 22;

	//const int     N_vec[18] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000 };
	//const int     num_N = 18;

	//const int     N_vec[14] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200 };
	//const int     num_N = 14;

	const double  alpha = 1; //if different from one, not flat prior
	const double  beta = 1; //only valid for alpha, beta >= 1 (alpha=beta=1, is the uniform case)
	
	//Definitions
	int	     i, j, k, i_max, N, M, L, iN, iM;
	int      num_samples;
	double   Q_max_th3, Q_max_th4;
	double   fact_na, fact_nb, fact_nab;
	double   cum, cum_pre;

	
	FILE *value_actions_th;

	value_actions_th = fopen("value_actions_th2th.m", "w");

	//seed

	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		for (iM = 1; iM <= N; iM++) {

			//that if, if the division gives an integer:
			if ( floor(N/iM)- (1.*N)/(1.*iM) == 0 ) {
					
				M = iM;

				num_samples = floor(N/M);

				/*
				//exact analytical expression
				L = num_samples;
				Q_max_th3 = 0.;
				for (i = 0; i <= L; i++) {
					Q_max_th3 = Q_max_th3 + (pow(i + 1, M) - pow(i, M))*(i+1);
				}
				Q_max_th3 = Q_max_th3 / (pow(L + 1, M)*(L + 2));
				*/

				//Another version of exact analytical expression for alpha=beta=1 (flat prior)
				L = num_samples;
				Q_max_th3 = 0.;
				for (i = 0; i <= L; i++) {
					Q_max_th3 = Q_max_th3 - pow(i + 1, M);
				}
				Q_max_th3 = 1. / (pow(L + 1, M)*(L + 2)) * (pow(L + 1, M+1) + pow(L + 1, M) + Q_max_th3);


				//Analytical expresion for alpha, beta arbitrary, but >1, and integer
				Q_max_th4 = 0; 
				cum_pre = 0;
				for (i = 0; i <= L; i++) {

					//cumulatives
					fact_na = 1;;
					for (k = 1; k < alpha; k++) {
						fact_na = fact_na * (i + alpha - k);
					}
					
					fact_nb = 1;;
					for (k = 1; k < beta; k++) {
						fact_nb = fact_nb * (L - i + beta - k);
					}

					fact_nab = 1;;
					for (k = 1; k < alpha + beta; k++) {
						fact_nab = fact_nab * (L + alpha + beta - k);
					}

					cum = cum_pre + fact_na*fact_nb/fact_nab * tgamma(alpha + beta) / (tgamma(alpha)*tgamma(beta));

					Q_max_th4 = Q_max_th4 + (1.*i + alpha) / (1.*L + alpha + beta) * (pow(cum,M) - pow(cum_pre,M));

					cum_pre = cum;
				
				}
				Q_max_th4 = Q_max_th4;


				printf("%i %i %i %f %f %f  \n", N, M, L, Q_max_th3, Q_max_th4, cum);

				fprintf(value_actions_th, "%i %i %i %f %f \n", N, M, L, Q_max_th3, Q_max_th4);
			}

		}

	}

	fclose(value_actions_th);
	
	return 0;
}

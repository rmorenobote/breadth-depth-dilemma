// HERE we have the THEORY for N choices, alpha, beta >=1
// only theory, so it is very FAST
//here we make a loop over alphas and beta to compute C_critical (the max C for which M/C =1) and exponent for larger values of C_critical

#include <stdio.h> 
#include <stdlib.h>
//#include <math.h>
#include "randomc.h"                   // define classes for random number generators
#include "userintf.cpp"                // define system specific user interface
#include "mersenne.cpp"                // members of class TRandomMersenne


int main(void)
{
	//Parameters
	const int     N_ini = 1;
	const int     N_end = 2000;
	const int     num_N = 2000;

	//loop for beta
	const double  alpha_0 = 1; //if different from one, not flat prior

	const double  beta_ini = 1; //only valid for alpha, beta >= 1 (alpha=beta=1, is the uniform case)
	const double  beta_end = 20;

	//loop for alpha
	const double  beta_0 = 1; //if different from one, not flat prior

	const double  alpha_ini = 1; //only valid for alpha, beta >= 1 (alpha=beta=1, is the uniform case)
	const double  alpha_end = 20;

	
	//Definitions
	int	     i, j, k, i_max, N, M, L, iN, iM, M_opt;
	int      num_samples, i_alpha, i_beta;
	double   Q_max_th3, Q_max_th4, Q_opt;
	double   fact_na, fact_nb, fact_nab;
	double   cum, cum_pre, alpha, beta;

	
	FILE *value_actions_th, *value_actions_th_opt;

	value_actions_th = fopen("value_actions_th2th_loop_alphas.m", "w");
	value_actions_th_opt = fopen("value_actions_th2th_loop_alphas_opt.m", "w");


	//seed

	//loop for beta, below loop for alpha
	for (i_beta = beta_ini; i_beta < beta_end + 1; i_beta++) {

		alpha = alpha_0;
		beta = i_beta;

		for (iN = 0; iN < num_N; iN++) {

			N = iN + 1;

			//innitial conditions
			Q_opt = 0;
			M_opt = 0;

			for (iM = 1; iM <= N; iM++) {

				//that is, if the division gives an integer, except for low values:
				if (N > 100 | floor(N / iM) - (1.*N) / (1.*iM) == 0) {

					M = iM;

					num_samples = floor(N / M); //aproximate to the closet lower integer

					//Exact analytical expression for alpha=beta=1 (flat prior)
					L = num_samples;
					Q_max_th3 = 0.;
					for (i = 0; i <= L; i++) {
						Q_max_th3 = Q_max_th3 - pow(i + 1, M);
					}
					Q_max_th3 = 1. / (pow(L + 1, M)*(L + 2)) * (pow(L + 1, M + 1) + pow(L + 1, M) + Q_max_th3);


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

						cum = cum_pre + fact_na * fact_nb / fact_nab * tgamma(alpha + beta) / (tgamma(alpha)*tgamma(beta));

						Q_max_th4 = Q_max_th4 + (1.*i + alpha) / (1.*L + alpha + beta) * (pow(cum, M) - pow(cum_pre, M));

						cum_pre = cum;

					}
					Q_max_th4 = Q_max_th4;

					if (Q_max_th4 > Q_opt) {
						Q_opt = Q_max_th4; //optimal value of reward across all M
						M_opt = M; //optimal M
					}

					//fprintf(value_actions_th, "%f %i %i %i %f %f \n", alpha/(alpha+beta), N, M, L, Q_max_th3, Q_max_th4);

				}

			}

		printf("%f %i %i %f \n", alpha / (alpha + beta), N, M_opt, Q_opt);

		fprintf(value_actions_th_opt, "%f %i %i %f  \n", alpha / (alpha + beta), N, M_opt, Q_opt);

		}

	}


	//loop for alpha
	for (i_alpha = alpha_ini+1; i_alpha < alpha_end + 1; i_alpha++) { // we start from alpha_ini+1 so that we do not repeat a=b=1

		alpha = i_alpha;
		beta = beta_0;

		for (iN = 0; iN < num_N; iN++) {

			N = iN + 1;

			//innitial conditions
			Q_opt = 0;
			M_opt = 0;

			for (iM = 1; iM <= N; iM++) {

				//that is, if the division gives an integer, except for low values:
				if (N > 100 | floor(N / iM) - (1.*N) / (1.*iM) == 0) {

					M = iM;

					num_samples = floor(N / M); //aproximate to the closet lower integer

												//Exact analytical expression for alpha=beta=1 (flat prior)
					L = num_samples;
					Q_max_th3 = 0.;
					for (i = 0; i <= L; i++) {
						Q_max_th3 = Q_max_th3 - pow(i + 1, M);
					}
					Q_max_th3 = 1. / (pow(L + 1, M)*(L + 2)) * (pow(L + 1, M + 1) + pow(L + 1, M) + Q_max_th3);


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

						cum = cum_pre + fact_na * fact_nb / fact_nab * tgamma(alpha + beta) / (tgamma(alpha)*tgamma(beta));

						Q_max_th4 = Q_max_th4 + (1.*i + alpha) / (1.*L + alpha + beta) * (pow(cum, M) - pow(cum_pre, M));

						cum_pre = cum;

					}
					Q_max_th4 = Q_max_th4;

					if (Q_max_th4 > Q_opt) {
						Q_opt = Q_max_th4; //optimal value of reward across all M
						M_opt = M; //optimal M
					}

					//fprintf(value_actions_th, "%f %i %i %i %f %f \n", alpha/(alpha+beta), N, M, L, Q_max_th3, Q_max_th4);

				}

			}

			printf("%f %i %i %f \n", alpha / (alpha + beta), N, M_opt, Q_opt);

			fprintf(value_actions_th_opt, "%f %i %i %f  \n", alpha / (alpha + beta), N, M_opt, Q_opt);

		}

	}


	fclose(value_actions_th);
	fclose(value_actions_th_opt);

	
	return 0;
}

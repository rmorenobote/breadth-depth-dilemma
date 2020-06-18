// here we have N options and we select M options (sequential?/focused), to compare with sample all of them once (parallel)
// here we introduce different non-uniform priors beta for prob, for arbitrary alpha and beta
//here we also have the analytic for any beta dist, alpha, beta >= 1, integer

#include <stdio.h> 
#include <stdlib.h>
//#include <math.h>
#include "randomc.h"                   // define classes for random number generators
#include "userintf.cpp"                // define system specific user interface
#include "mersenne.cpp"                // members of class TRandomMersenne


int main(void)
{
	//Parameters
	const int     N_vec[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000};
	const int     num_N = 20;

	//const int     N_vec[18] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000 };
	//const int     num_N = 18;

	//const int     N_vec[14] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200 };
	//const int     num_N = 14;

	//const int     N_vec[12] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20 };
	//const int     num_N = 12;

	//const int     N_vec[2] = { 2, 4};
	//const int     num_N = 2;

	const double  alpha = 1; //if different from one, not flat prior
	const double  beta = 1; //only valid for alpha, beta >= 1 (alpha=beta=1, is the uniform case)

	const int     per = 1;

	const int     samples = 10000; //10000 is a good value

	const int     seed =  46021; //random seed
	
	//Definitions
	int	     i, j, k, i_max, N, M, L, iN, iM, n, n_aux, i_aux, flag;
	int      num_samples; //this is L in our math notation
	double   r, sum, aux;
	double   r1_sum, r2_sum, r1, r2;
	double   Q_max, Q_max_aver, p_max, p_select_aver, p_max_aver;
	double   Q_max_th2, Q_max_th3, Q_max_th4;
	double   fact_na, fact_nb, fact_nab;
	double   cum, cum_pre;
	double   R, R_aver; 

	int      *y_vec, *i_vec;
	double   *Q_vec, *p_vec;

	
	FILE *value_actions, *value_actions_per, *value_actions_th;

	value_actions = fopen("value_actions2.m", "w");
	value_actions_per = fopen("value_actions_per2.m", "w");
	value_actions_th = fopen("value_actions_th2.m", "w");

	//seed
    TRandomMersenne rg(seed);

	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		for (iM = 1; iM <= N; iM++) {

			//that if, if the division gives an integer:
			if ( floor(N/iM)- (1.*N)/(1.*iM) == 0 ) {
					
				M = iM;

				num_samples = floor(N/M);

				//Allocation of vectors
				y_vec = (int*) malloc(N*sizeof(double));
				Q_vec = (double*) malloc(N*sizeof(double));
				p_vec = (double*) malloc(N*sizeof(double));

				//initial condiction
				Q_max_aver = 0;
				p_select_aver = 0;
				R_aver = 0;

				
				//SIMULATIONS monte carlo simulations
				for (k = 0; k < samples; k++) {

					for (i = 0; i < N; i++) {
						y_vec[i] = 0; 

						flag = 0;
						while (flag == 0) {
							r1 = rg.Random();
							r2 = rg.Random();
							if (r2 <= pow(1.*r1, alpha - 1.)*pow(1. - r1, beta - 1.)) {
								flag = 1;
							}	
						}
						p_vec[i] = r1;

						//printf("%i %i %f \n", i, flag, p_vec[i]);

						/*
						//drawing uniform p for all Bernouillis from a beta
						r1_sum = 0;
						r2_sum = 0;
						for (j = 0; j < floor(alpha); j++) {
							r = rg.Random();
							r1_sum = r1_sum - log(r);
							r = rg.Random();
							r2_sum = r2_sum - log(r);
						}
						p_vec[i] = r1_sum / (r1_sum + r2_sum);
						*/

						//drawing uniform p for all Bernouillis
						//r = rg.Random();
						//p_vec[i] = r;
					}

					//for each sampled Bernoilli, we take num_samples
					for (i = 0; i < M; i++) {

						Q_max = 0;

						for (j = 0; j < num_samples; j++) {
							r = rg.Random();
							if(r < p_vec[i]) {
								y_vec[i] = y_vec[i] + 1;
							}
						}

						//Q values:
						Q_vec[i] = (alpha + y_vec[i])/(alpha + beta + num_samples); //predictive posterior
					}

					p_max = 0;
					for (i = 0; i < N; i++) {
						if (p_vec[i] > p_max) {
							p_max = p_vec[i];
						}
					}

					//check this 0 instead of 1/2, and also below, changes a lot!!!!!
					Q_max = 0; // 1. / 2.; //maybe here put 1/2, because we can always take one of the not-sampled options, with reward prob 1/2
					for (i = 0; i < M; i++) {
						if (Q_vec[i] > Q_max) {
							Q_max = Q_vec[i];
							i_max = i;
						}
					}

					Q_max_aver = Q_max_aver + Q_max;
					p_select_aver = p_select_aver + p_vec[i_max];
					p_max_aver = p_max_aver + p_max;

					//simulating reward
					r = rg.Random();
					if (r < p_vec[i_max]) {
						R = 1;
					}
					if (r >= p_vec[i_max]) {
						R = 0;
					}
					R_aver = R_aver + R;

					//printf("%i %f %f %f \n", k, Q_max, p_vec[i_max], p_max);	
				}
				//printf("\n");

				Q_max_aver = Q_max_aver/samples;
				p_select_aver = p_select_aver/samples;
				p_max_aver = p_max_aver/samples;
				R_aver = R_aver / samples;

				//note that num_samples = L in our math
				printf("%i %i %i %f %f %f %f \n", N, M, num_samples, p_max_aver, Q_max_aver, p_select_aver, R_aver);

				fprintf(value_actions, "%i %i %i %f %f %f %f \n", N, M, num_samples, p_max_aver, Q_max_aver, p_select_aver, R_aver);

				//approx analytical expression
				Q_max_th2 = (M / (M + 1.) + 1. / (M + 1.)*pow(0.5, M + 1))*((num_samples+alpha)/(num_samples+alpha+beta)); //hack, approx

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
				//alternative
				//Q_max_th3 = 1. + Q_max_th3/(pow(L + 1, M)*(L + 2));


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


				printf("%i %i %i %f %f %f %f %f \n", N, M, L, Q_max_aver, Q_max_th2, Q_max_th3, Q_max_th4, cum);

				fprintf(value_actions_th, "%i %i %i %f %f %f \n", N, M, L, Q_max_th2, Q_max_th3, Q_max_th4);
			}

		}

	}

	//this just tests the perturbation L+1 for the 1 option, and L-1 for the second option, and the others with L samples, as above
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		for (iM = 1; iM <= N; iM++) {

			//that if, if the division gives an integer:
			if (floor(N / iM) - (1.*N) / (1.*iM) == 0) {

				M = iM;

				num_samples = floor(N / M);

				//Allocation of vectors
				y_vec = (int*)malloc(N * sizeof(double));
				Q_vec = (double*)malloc(N * sizeof(double));
				p_vec = (double*)malloc(N * sizeof(double));

				//initial condiction
				Q_max_aver = 0;
				p_select_aver = 0;
				R_aver = 0;


				//SIMULATIONS monte carlo simulations
				for (k = 0; k < samples; k++) {

					for (i = 0; i < N; i++) {
						y_vec[i] = 0;

						flag = 0;
						while (flag == 0) {
							r1 = rg.Random();
							r2 = rg.Random();
							if (r2 <= pow(r1, alpha - 1)*pow(1 - r1, beta - 1)) {
								flag = 1;
							}
						}
						p_vec[i] = r1;

						/*
						//drawing uniform p for all Bernouillis from a beta
						r1_sum = 0;
						r2_sum = 0;
						for (j = 0; j < floor(alpha); j++) {
							r = rg.Random();
							r1_sum = r1_sum - log(r);
							r = rg.Random();
							r2_sum = r2_sum - log(r);
						}
						p_vec[i] = r1_sum / (r1_sum + r2_sum);
						*/

						//drawing uniform p for all Bernouillis
						//r = rg.Random();
						//p_vec[i] = r;
					}


					//for each sampled Bernoilli, we take num_samples
					for (i = 0; i < M; i++) {

						//we only take L+1 and L-1 if there are at least to alternatives
						if (M >= 2) {

							Q_max = 0;

							if (i == 0) {
								for (j = 0; j < num_samples + per; j++) { //we take here L+1 samples
									r = rg.Random();
									if (r < p_vec[i]) {
										y_vec[i] = y_vec[i] + 1;
									}
								}

								//Q values:
								Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + num_samples + per); //predictive posterior
							}
							else if (i == 1) {
								for (j = 0; j < num_samples - per; j++) { //we take here L-1 samples
									r = rg.Random();
									if (r < p_vec[i]) {
										y_vec[i] = y_vec[i] + 1;
									}
								}

								//Q values:
								Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + num_samples - per); //predictive posterior
							}
							else {
								for (j = 0; j < num_samples; j++) {
									r = rg.Random();
									if (r < p_vec[i]) {
										y_vec[i] = y_vec[i] + 1;
									}
								}

								//Q values:
								Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + num_samples); //predictive posterior
							}
						}

						if (M < 2) {

							Q_max = 0;

							for (j = 0; j < num_samples; j++) {
								r = rg.Random();
								if (r < p_vec[i]) {
									y_vec[i] = y_vec[i] + 1;
								}
							}

							//Q values:
							Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + num_samples); //predictive posterior
						}

					}

					p_max = 0;
					for (i = 0; i < N; i++) {
						if (p_vec[i] > p_max) {
							p_max = p_vec[i];
						}
					}

					//check this 0 instead of 1/2, and also below, changes a lot!!!!!
					Q_max = 0; // 1. / 2.; //maybe here put 1/2, because we can always take one of the not-sampled options, with reward prob 1/2
					for (i = 0; i < M; i++) {
						if (Q_vec[i] > Q_max) {
							Q_max = Q_vec[i];
							i_max = i;
						}
					}

					Q_max_aver = Q_max_aver + Q_max;
					p_select_aver = p_select_aver + p_vec[i_max];
					p_max_aver = p_max_aver + p_max;

					//simulating reward
					r = rg.Random();
					if (r < p_vec[i_max]) {
						R = 1;
					}
					if (r >= p_vec[i_max]) {
						R = 0;
					}
					R_aver = R_aver + R;

					//printf("%i %f %f %f \n", k, Q_max, p_vec[i_max], p_max);	
				}
				//printf("\n");

				Q_max_aver = Q_max_aver / samples;
				p_select_aver = p_select_aver / samples;
				p_max_aver = p_max_aver / samples;
				R_aver = R_aver / samples;

				printf("%i %i %i %f %f %f %f \n", N, M, num_samples, p_max_aver, Q_max_aver, p_select_aver, R_aver);

				fprintf(value_actions_per, "%i %i %i %f %f %f %f \n", N, M, num_samples, p_max_aver, Q_max_aver, p_select_aver, R_aver);

			}

		}

	}

	fclose( value_actions );
	fclose(value_actions_per );
	fclose(value_actions_th);
	
	return 0;
}

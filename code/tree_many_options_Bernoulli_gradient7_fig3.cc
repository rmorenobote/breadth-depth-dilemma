// here we have N options and we select M options (sequential?/focused), to compare with sample all of them once (parallel)
// here we introduce different non-uniform priors beta for prob, for arbitrary alpha and beta
//here we USE STOCHASTIC GRADIENT DESCENT to looks for the OPTIMAL distribution of resources
//here we start the search from the uniform distribution, or close to it (identical code to ...gradient6, no difference, just to separate results)
//We only take non-increasing distributions!!!!


#include <stdio.h> 
#include <stdlib.h>
//#include <math.h>
#include "randomc.h"                   // define classes for random number generators
#include "userintf.cpp"                // define system specific user interface
#include "mersenne.cpp"                // members of class TRandomMersenne


int main(void)
{
	//Parameters

	//const int     N_vec[18] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000 };
	//const int     num_N = 18;

	//use squares to compare
	const int     N_vec[18] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 25, 49, 100, 14 * 14, 20 * 20, 30 * 30, 50*50 };
	const int     num_N = 18;

	//use squares to compare
	//const int     N_vec[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 25, 49, 100, 14*14 };
	//const int     num_N = 15;

	//const int     N_vec[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	//const int     num_N = 10;

	//use squares to compare
	//const int     N_vec[20] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
	//const int     num_N = 20;



	const double  alpha = 1; //if different from one, not flat prior
	const double  beta = 1; //only valid for alpha, beta >= 1 (alpha=beta=1, is the uniform case)

	const int     iterations = 3000;

	const int     samples = 50000;

	const int     seed =  40212; //random seed
	
	//Definitions
	int	     i, j, k, i_max, N, M, L, iN, iM, n, n_aux, i_aux, flag, iter;
	int      i1, i2, N_res;
	double   r, sum, aux, r1, r2, Q_max_ver_ini;
	double   Q_max, Q_max_aver, Q_max_aver_old, Q_max_aver_new, Q_max_aver_ini, p_select_aver;
	double   Q_max_th3, Q_max_th32;
	double   R, R_aver; 

	int      *y_vec, *i_vec, *L_vec;
	double   *Q_vec, *p_vec;

	
	FILE *value_actions, *distributions;

	value_actions = fopen("value_actions_gradient7.m", "w");
	distributions = fopen("distributions7.m", "w");

	//seed
    TRandomMersenne rg(seed);

	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		//allocation of vector
		L_vec = (int*)malloc(N * sizeof(double)); //vector of samples used for option 1, 2, ... N
		y_vec = (int*)malloc(N * sizeof(double));
		Q_vec = (double*)malloc(N * sizeof(double));
		p_vec = (double*)malloc(N * sizeof(double));

		/*
		//this is the uniform initial condition
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 1;
		}
		*/
		//sqrt initial Sol, or close to it
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 0; //to start with
		}
		L = floor(sqrt(N));
		for (iM = 0; iM < L; iM++) {
			L_vec[iM] = L;
		}
		N_res = N - pow(L, 2); //residual because sqrt is not exact partitioning
		if (N - pow(L, 2) > 0) {
			L_vec[L] = N - pow(L, 2);
		}
		//uniform L=1 for N<=7
		if (N <= 7) {
			for (iM = 0; iM < N; iM++) {
				L_vec[iM] = 1; //to start with
			}
		}

		//Gradient iterations
		for (iter = 0; iter < iterations; iter++) {

			/////////////////
			//L_vec initial

			//Computing the value of the action L_vec

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
				}

				//for each sampled Bernoilli, we take num_samples L_vec
				for (i = 0; i < N; i++) {

					Q_max = 0;

					for (j = 0; j < L_vec[i]; j++) {
						r = rg.Random();
						if (r < p_vec[i]) {
							y_vec[i] = y_vec[i] + 1;
						}
					}

					//Q values:
					Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + L_vec[i]); //predictive posterior

					//that is, if we do not sample this option, then we cannot take it:
					if (L_vec[i] == 0) {
						Q_vec[i] = 0;
					}
				}

				//check this 0 instead of 1/2, and also below, changes a lot!!!!!
				Q_max = 0; // 1. / 2.; //maybe here put 1/2, because we can always take one of the not-sampled options, with reward prob 1/2
				for (i = 0; i < N; i++) {
					if (Q_vec[i] > Q_max) {
						Q_max = Q_vec[i];
						i_max = i;
					}
				}
				Q_max_aver = Q_max_aver + Q_max;

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

			Q_max_aver_old = Q_max_aver / samples;
			p_select_aver = p_select_aver / samples;
			R_aver = R_aver / samples;
			if (iter == 0) {
				Q_max_aver_ini = Q_max_aver_old; //initial reward for uniform distribution with sqrt rule
			}


			/////////////////
			//L_vec perturbed

			flag = 0;
			while (flag == 0) {
				//chose randomly two options, and increase and decrease between them, if L_vec>1
				r = rg.Random();
				i1 = floor(1.*r*N);
				r = rg.Random();
				i2 = floor(1.*r*N);

				// we only take it if it is positive so we can remove one sample from this option
					//    ...and if L_vec is non-decreasing
				if (L_vec[i1] > 0) {
					if (i1 >= i2) {
						flag = 1;
					}
					if (i1 < i2) {
						if (L_vec[i1] - 1 >= L_vec[i2] + 1 && L_vec[i2-1] >= L_vec[i2] + 1) {
							flag = 1;
						}
					}

				}
			}

			L_vec[i1] = L_vec[i1] - 1;
			L_vec[i2] = L_vec[i2] + 1;

			//Computing the value of the action L_vec

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
				}

				//for each sampled Bernoilli, we take num_samples L_vec
				for (i = 0; i < N; i++) {

					Q_max = 0;

					for (j = 0; j < L_vec[i]; j++) {
						r = rg.Random();
						if (r < p_vec[i]) {
							y_vec[i] = y_vec[i] + 1;
						}
					}

					//Q values:
					Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + L_vec[i]); //predictive posterior

					//that is, if we do not sample this option, then we cannot take it:
					if (L_vec[i] == 0) {
						Q_vec[i] = 0; 
					}
				}


				//check this 0 instead of 1/2, and also below, changes a lot!!!!!
				Q_max = 0; // 1. / 2.; //maybe here put 1/2, because we can always take one of the not-sampled options, with reward prob 1/2
				for (i = 0; i < N; i++) {
					if (Q_vec[i] > Q_max) {
						Q_max = Q_vec[i];
						i_max = i;
					}
				}
				Q_max_aver = Q_max_aver + Q_max;

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

			Q_max_aver_new = Q_max_aver / samples;
			p_select_aver = p_select_aver / samples;
			R_aver = R_aver / samples;

			//we revert change of L_vec if it did not improve Q_max
			if (Q_max_aver_new < Q_max_aver_old) {
				L_vec[i1] = L_vec[i1] + 1;
				L_vec[i2] = L_vec[i2] - 1;
			}
			
			//printf("%i %i %f %f %i %i %i %i \n", N, iter, Q_max_aver_new, Q_max_aver_old, L_vec[0], L_vec[1], L_vec[2], L_vec[3]);
		
		}

		//exact analytical expression for the approximate optimal M*=sqrt(C), to compare with the optimal found by stoch grad descent
		L = floor(sqrt(N));
		Q_max_th3 = 0.;
		for (i = 0; i <= L; i++) {
			Q_max_th3 = Q_max_th3 + (pow(i + 1, L) - pow(i, L))*(i + 1);
		}
		Q_max_th3 = Q_max_th3 / (pow(L + 1, L)*(L + 2));

		//uniform choice
		M = N; //L=1 here
		Q_max_th32 = 0.;
		for (i = 0; i <= 1; i++) {
			Q_max_th32 = Q_max_th32 + (pow(i + 1, M) - pow(i, M))*(i + 1);
		}
		Q_max_th32 = Q_max_th32 / (pow(1 + 1, M)*(1 + 2));

		printf("%i %i %i %f %f %f %f %f \n", N, L, iter, Q_max_aver_new, Q_max_aver_old, Q_max_th3, Q_max_th32, Q_max_aver_ini);

		fprintf(value_actions, "%i %i %i %f %f %f %f %f %f \n", N, L, iter, Q_max_aver_new, R_aver, Q_max_aver_old, Q_max_th3, 
													Q_max_th32, Q_max_aver_ini);

		for (i = 0; i < N; i++) {
			printf("%i %i %i \n", N, i, L_vec[i]);
			fprintf(distributions, "%i %i %i \n", N, i, L_vec[i]);
		}

	}

	fclose( value_actions );
	fclose(distributions);
	
	return 0;
}

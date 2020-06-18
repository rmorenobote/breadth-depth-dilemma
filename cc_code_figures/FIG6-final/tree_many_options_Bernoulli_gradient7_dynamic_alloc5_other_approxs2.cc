// here we have N options and we select M options (sequential?/focused), to compare with sample all of them once (parallel)
// here we introduce different non-uniform priors beta for prob, for arbitrary alpha and beta
//here we USE STOCHASTIC GRADIENT DESCENT to looks for the OPTIMAL distribution of resources
//here we start the search from the uniform distribution, or close to it (identical code to ...gradient6, no difference, just to separate results)
//We only take non-increasing distributions!!!!
//Here we run the DYNAMIC ALLOCATION analaysis, with WAVES
//in this code, we resampled the top M_vec[j] of the previous WAVE; we allow for non-decreasing M_i


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
	//const int     N_vec[18] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 25, 49, 100, 14 * 14, 20 * 20, 30 * 30, 50*50 };
	//const int     num_N = 18;

	//use squares to compare
	//const int     N_vec[15] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 25, 49, 100, 14*14 };
	//const int     num_N = 15;

	//const int     N_vec[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
	//const int     num_N = 10;

	//use squares to compare
	//const int     N_vec[20] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
	//const int     num_N = 20;

	//use squares to compare
	const int     N_vec[27] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 49, 100, 14 * 14, 20 * 20, 30 * 30, 50 * 50 };
	const int     num_N = 27;

	const double  alpha = 1; //if different from one, not flat prior
	const double  beta = 1; //only valid for alpha, beta >= 1 (alpha=beta=1, is the uniform case)

	const int     iterations = 3000;

	const int     samples = 50000;

	const int     seed =  4021; //random seed
	
	//Definitions
	int	     i, j, k, i_max, N, M, L, iN, iM, n, n_aux, i_aux, flag, flag_aux, iter, N_aux;
	int      i1, i2, N_res, a, b, a_L, b_L, i_index, sum_L;
	int      L_max, L_max_per, M_max, M_max_per, index;
	double   r, sum, aux, r1, r2, Q_max_ver_ini, p_a, p_b, q_a, q_b;
	double   Q_max, Q_max_aver, Q_max_aver_old, Q_max_aver_new, Q_max_aver_ini, p_select_aver;
	double   Q_max_th3, Q_max_th32;
	double   R, R_aver; 

	int      *y_vec, *i_vec, *M_vec, *L_vec;
	double   *Q_vec, *p_vec;

	
	FILE *value_actions, *value_actions_breadth, *value_actions_depth,
		*value_actions_random, *value_actions_triangular, *value_actions_sq_root, *distributions;

	value_actions = fopen("value_actions_gradient7_dyn_alloc5_other_approxs2.m", "w");
	value_actions_breadth = fopen("value_actions_gradient7_dyn_alloc5_other_approxs_breadth2.m", "w");
	value_actions_depth = fopen("value_actions_gradient7_dyn_alloc5_other_approxs_depth2.m", "w");
	value_actions_random = fopen("value_actions_gradient7_dyn_alloc5_other_approxs_random2.m", "w");
	value_actions_triangular = fopen("value_actions_gradient7_dyn_alloc5_other_approxs_triangular2.m", "w");
	value_actions_sq_root = fopen("value_actions_gradient7_dyn_alloc5_other_approxs_sq_root2.m", "w");
	distributions = fopen("distributions7_dyn_alloc5_other_approxs2.m", "w");

	//seed
    TRandomMersenne rg(seed);


	//pure Breadth search
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		//allocation of vector
		L_vec = (int*)malloc(N * sizeof(double)); //vector of samples used for option 1, 2, ... N
		y_vec = (int*)malloc(N * sizeof(double));
		Q_vec = (double*)malloc(N * sizeof(double));
		p_vec = (double*)malloc(N * sizeof(double));

		//innitial condition
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 0;
		}

		//breadth search
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 1;
		}

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
		}

		Q_max_aver = Q_max_aver / samples;

		printf("%i %f \n", N, Q_max_aver);

		fprintf(value_actions_breadth, "%i %f \n", N, Q_max_aver);

	}


	//pure Depth search (2 options sampled, at most)
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		//allocation of vector
		L_vec = (int*)malloc(N * sizeof(double)); //vector of samples used for option 1, 2, ... N
		y_vec = (int*)malloc(N * sizeof(double));
		Q_vec = (double*)malloc(N * sizeof(double));
		p_vec = (double*)malloc(N * sizeof(double));

		//innitial condition
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 0;
		}

		//depth search (2 options max)
		L = floor(N / 2);
		if (N > 1) {
			L_vec[0] = L;
			L_vec[1] = N - L;
		}
		if (N == 1) {
			L_vec[0] = 1;
		}

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
		}

		Q_max_aver = Q_max_aver / samples;

		printf("%i %f %i \n", N, Q_max_aver, L);

		fprintf(value_actions_depth, "%i %f %i \n", N, Q_max_aver, L);

	}



	//pure Random (with replacement) search
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		//allocation of vector
		L_vec = (int*)malloc(N * sizeof(double)); //vector of samples used for option 1, 2, ... N
		y_vec = (int*)malloc(N * sizeof(double));
		Q_vec = (double*)malloc(N * sizeof(double));
		p_vec = (double*)malloc(N * sizeof(double));

		//Computing the value of the action L_vec

		//initial condiction
		Q_max_aver = 0;
		p_select_aver = 0;
		R_aver = 0;

		//SIMULATIONS monte carlo simulations
		for (k = 0; k < samples; k++) {

			//A different random allocation each time:
			//innitial condition
			for (iM = 0; iM < N; iM++) {
				L_vec[iM] = 0;
			}
			//loading L_vec
			for (i = 0; i < N; i++) {
				r = rg.Random();
				i_index = floor(1.*r*N);
				L_vec[i_index] = L_vec[i_index] + 1;
			}

			//check
			sum_L = 0;
			for (iM = 0; iM < N; iM++) {
				sum_L = sum_L + L_vec[iM];
			}
			if (sum_L != N) {
				printf("%i %i %i \n", -9999, sum_L, N);
			}

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
		}

		Q_max_aver = Q_max_aver / samples;

		printf("%i %f %i \n", N, Q_max_aver, L);

		fprintf(value_actions_random, "%i %f \n", N, Q_max_aver);

	}


	//triangular search
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		//allocation of vector
		L_vec = (int*)malloc(N * sizeof(double)); //vector of samples used for option 1, 2, ... N
		y_vec = (int*)malloc(N * sizeof(double));
		Q_vec = (double*)malloc(N * sizeof(double));
		p_vec = (double*)malloc(N * sizeof(double));

		//innitial condition
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 0;
		}

		//triangular search
		if (N > 5) {
			L = floor(sqrt(2 * N));
			sum_L = 0;
			for (iM = 0; iM < L; iM++) { //asign L - i capacity to each alternative i until capacity is exhausted
				L_vec[iM] = L - iM - 1; //note the minus -1
				sum_L = sum_L + L_vec[iM];
				if (sum_L >= N) {
					L_vec[iM] = L_vec[iM] - (sum_L - N);
					exit;
				}
			}
			if (sum_L < N) {
				L_vec[iM + 1] = N - sum_L;
			}
		}
		if (N <= 5) {
			L = N;
			for (iM = 0; iM < N; iM++) {
				L_vec[iM] = 1;
			}
		}

		//check
		sum_L = 0;
		for (iM = 0; iM < N; iM++) {
			sum_L = sum_L + L_vec[iM];
		}
		if (sum_L != N) {
			printf("%i %i %i \n", -1111, sum_L, N);
		}

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
		}

		Q_max_aver = Q_max_aver / samples;

		printf("%i %f \n", N, Q_max_aver);

		fprintf(value_actions_triangular, "%i %f \n", N, Q_max_aver);

	}


	//square root search
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		//allocation of vector
		L_vec = (int*)malloc(N * sizeof(double)); //vector of samples used for option 1, 2, ... N
		y_vec = (int*)malloc(N * sizeof(double));
		Q_vec = (double*)malloc(N * sizeof(double));
		p_vec = (double*)malloc(N * sizeof(double));

		//innitial condition
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 0;
		}

		//sqrt initial Sol, or close to it
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 0; //all zeros
		}
		L = floor(sqrt(N));
		for (iM = 0; iM < L; iM++) {
			L_vec[iM] = L; //aprox uniform allocations
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
		}

		Q_max_aver = Q_max_aver / samples;

		printf("%i %f \n", N, Q_max_aver);

		fprintf(value_actions_sq_root, "%i %f \n", N, Q_max_aver);

	}


	/// GRADIENT DESCENT METHOD
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		//allocation of vector
		M_vec = (int*)malloc(N * sizeof(double)); //vector of samples used in WAVEs 1, 2, ... N
		L_vec = (int*)malloc(N * sizeof(double)); //accumulated samples for each options
		y_vec = (int*)malloc(N * sizeof(double));
		Q_vec = (double*)malloc(N * sizeof(double));
		p_vec = (double*)malloc(N * sizeof(double));

		//sqrt initial Sol, or close to it
		for (iM = 0; iM < N; iM++) {
			L_vec[iM] = 0; //all zeros
		}
		L = floor(sqrt(N));
		for (iM = 0; iM < L; iM++) {
			L_vec[iM] = L; //aprox uniform allocations
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

		//we now transform L_vec initial into M_vec initial, compatible with the above aprox square root 
		for (i = 0; i < N; i++) {
			M_vec[i] = 0;
		}
		for (i = 0; i < N; i++) { //M_vec is decreasing even if L_vec is not (M_j = {i: L_i>=j}) /note that it is > because i,j start at zero
			for (j = 0; j < N; j++) {
				if (j < L_vec[i]) {
					M_vec[j] = M_vec[j] + 1;
				}
			}
		}
		//for (i = 0; i < N; i++) {
		//	printf("%i %i %i %i %i \n", N, i, L_vec[i], M_vec[i], -9999);
		//}		

		//Gradient iterations
		for (iter = 0; iter < iterations; iter++) {

			/////////////////
			//M_vec initial

			//Computing the value of the action M_vec

			//initial condition
			Q_max_aver = 0;
			p_select_aver = 0;
			R_aver = 0;

			//SIMULATIONS monte carlo simulations
			for (k = 0; k < samples; k++) {

				for (i = 0; i < N; i++){
					y_vec[i] = 0;
					L_vec[i] = 0;
					Q_vec[i] = 0;
				}

				for (i = 0; i < N; i++) { 

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

				//maximum number of waves, that depends on M_vec
				L_max = 0;
				for (i = 0; i < N; i++) { //this is the rule L_j = {i: M_i>=j}
					if (M_vec[i] > 0) {
						L_max = L_max + 1;
					}
				}
				//maximum in the number of options sampled in any wave
				M_max = 0;
				for (i = 0; i < N; i++) { //this is the rule L_j = {i: M_i>=j}
					if (M_vec[i] > M_max) {
						M_max = M_vec[i];
					}
				}

				//loop over WAVES: max number of waves = N
				for (j = 0; j < L_max; j++) { 

					//order the options from highest number of positive outcomes
					flag = 0;
					while (flag == 0) {
						flag = 1; // assume order is correct
						for (i = 0; i < M_max-1; i++) {
							if (y_vec[i + 1] > y_vec[i]) {
								a = y_vec[i];
								b = y_vec[i + 1]; 
								a_L = L_vec[i];
								b_L = L_vec[i + 1];
								p_a = p_vec[i];
								p_b = p_vec[i + 1]; 
								q_a = Q_vec[i];
								q_b = Q_vec[i + 1];

								//we change order
								y_vec[i] = b; //we change values if y was not ordered from high to low
								y_vec[i + 1] = a;
								L_vec[i] = b_L; //we change values if y was not ordered from high to low
								L_vec[i + 1] = a_L;
								p_vec[i] = p_b;
								p_vec[i + 1] = p_a;
								Q_vec[i] = q_b;
								Q_vec[i + 1] = q_a;

								flag = 0; //incorrect order
							}
						}
					}
					/*
					if (j == N - 1 & k == samples -1 & iter == iterations - 1) {
						for (i = 0; i < N; i++) {
							printf("%i %i %i %i %i %i %i %i  \n", N, iter, k, j, i, y_vec[i], L_vec[i], -9999);
						}
					}
					*/
				
					
					//for each wave j , we sampled the TOP/Higuest M_vec[j] Bernoilli var, 
					//	and we take num_samples = 1 per sampled option
					for (i = 0; i < M_vec[j]; i++) {

						Q_max = 0;

						//we sample just 1 per sampled option:
						r = rg.Random();
						if (r < p_vec[i]) {
							y_vec[i] = y_vec[i] + 1;
						}
						L_vec[i] = L_vec[i] + 1;
					}
				}

				for (i = 0; i < N; i++) {
					if (L_vec[i] > 0) { //only sampled ones (note that non-sampled have Q_vec = 0)
						Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + L_vec[i]); //predictive posterior
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
				Q_max_aver_ini = Q_max_aver_old;
			}


			/////////////////
			//M_vec perturbed

			//limit looking for i1 and i2 only in the L_max + 1 interval, not in 1 to N
			if( L_max < N ){
				N_aux = L_max + 1;
			}
			if (L_max == N) {
				N_aux = L_max;
			}
			//printf("%i %i %i \n", N, L_max, N_aux);

			flag = 0;
			while (flag == 0) {
				//chose randomly two WAVES, and increase and decrease between them, if M_vec>1
				r = rg.Random();
				i1 = floor(1.*r*N_aux); 
				r = rg.Random();
				i2 = floor(1.*r*N_aux);

				if (N == 1) {
					M_vec[i1] = M_vec[i1] - 1;
					M_vec[i2] = M_vec[i2] + 1;

					flag = 1;
				}
				if (N > 1) {
					// we only take it if it is positive so we can remove one sample from this option
					if (M_vec[i1] > 0 && i1 != i2) {

						M_vec[i1] = M_vec[i1] - 1;
						M_vec[i2] = M_vec[i2] + 1;
					
						flag_aux = 1;
						//check now that it is non-increasing
						for (i = 0; i < N-1; i++) {
							if (M_vec[i] < M_vec[i + 1]) { //if not non-increasing, then look for another proposal
								flag_aux = 0;
								exit;
							}
						}
						if (flag_aux == 0){
							//undo if it was not succesful above
							M_vec[i1] = M_vec[i1] + 1;
							M_vec[i2] = M_vec[i2] - 1;
						}
						flag = flag_aux;
					}
				}
			}

			//Computing the value of the action M_vec
			//initial condition
			Q_max_aver = 0;
			p_select_aver = 0;
			R_aver = 0;

			//SIMULATIONS monte carlo simulations
			for (k = 0; k < samples; k++) {

				for (i = 0; i < N; i++) {
					y_vec[i] = 0;
					L_vec[i] = 0;
					Q_vec[i] = 0;
				}

				for (i = 0; i < N; i++) { //we only sample M_vec[0] options in the first wave, and we stick on them

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

				//maximum number of waves, that depends on M_vec
				L_max_per = 0;
				for (i = 0; i < N; i++) { //this is the rule L_j = {i: M_i>=j}
					if (M_vec[i] > 0) {
						L_max_per = L_max_per + 1;
					}
				}
				//maximum in the number of options sampled in any wave
				M_max_per = 0;
				for (i = 0; i < N; i++) { //this is the rule L_j = {i: M_i>=j}
					if (M_vec[i] > M_max_per) {
						M_max_per = M_vec[i];
					}
				}


				//loop over WAVES: max number of waves = N
				for (j = 0; j < L_max_per; j++) {

					//order the options from highest number of positive outcomes
					flag = 0;
					while (flag == 0) {
						flag = 1; // assume order is correct
						for (i = 0; i < M_max_per-1; i++) {
							if (y_vec[i + 1] > y_vec[i]) {
								a = y_vec[i];
								b = y_vec[i + 1];
								a_L = L_vec[i];
								b_L = L_vec[i + 1];
								p_a = p_vec[i];
								p_b = p_vec[i + 1];
								q_a = Q_vec[i];
								q_b = Q_vec[i + 1];

								y_vec[i] = b; //we change values if y was not ordered from high to low
								y_vec[i + 1] = a;
								L_vec[i] = b_L; //we change values if y was not ordered from high to low
								L_vec[i + 1] = a_L;
								p_vec[i] = p_b;
								p_vec[i + 1] = p_a;
								Q_vec[i] = q_b;
								Q_vec[i + 1] = q_a;

								flag = 0; //wrong order
							}
						}
					}

					//for each wave j , we sampled the TOP M_vec[j] Bernoilli var, 
					//	and we take num_samples = 1 per sampled option
					for (i = 0; i < M_vec[j]; i++) {

						Q_max = 0;

						//we sample just 1 per sampled option:
						r = rg.Random();
						if (r < p_vec[i]) {
							y_vec[i] = y_vec[i] + 1;
						}
						L_vec[i] = L_vec[i] + 1;
					}
				}

				for (i = 0; i < N; i++) {
					if (L_vec[i] > 0) { //only sampled ones (note that non-sampled have Q_vec = 0)
						Q_vec[i] = (alpha + y_vec[i]) / (alpha + beta + L_vec[i]); //predictive posterior
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
				M_vec[i1] = M_vec[i1] + 1;
				M_vec[i2] = M_vec[i2] - 1;
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

		//converting M_vec into L_vec:
		for (i = 0; i < N; i++) {
			L_vec[i] = 0;
		}
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				if (j < M_vec[i]) {
					L_vec[j] = L_vec[j] + 1;
				}
			}
		}

		for (i = 0; i < N; i++) {
			printf("%i %i %i %i \n", N, i, L_vec[i], M_vec[i]);
			fprintf(distributions, "%i %i %i %i \n", N, i, L_vec[i], M_vec[i]);
		}

	}

	fclose( value_actions );
	fclose(distributions);
	
	return 0;
}

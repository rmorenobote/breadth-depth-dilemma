// here we have N options from a Gaussian
// we compute the optimum over M for different capacities
// here we have an analytical expression for the expeced value, computed numerically through integration
// THIS is a new program where Capacity is controlled by the number of samples, not by sigma!!!!!

#include <stdio.h> 
#include <stdlib.h>
//#include <math.h>
#include "randomc.h"                   // define classes for random number generators
#include "userintf.cpp"                // define system specific user interface
#include "mersenne.cpp"                // members of class TRandomMersenne


int main(void)
{
	//Parameters

	//const int     N_vec[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	//const int     num_N = 10;
	
	const int     N_vec[23] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 50, 100, 150, 200, 300, 400, 500, 1000, 2000, 5000 };
	const int     num_N = 23;

	const double  mu_0 = 1;
	const double  sigma2_0 = 1;
	const double  sigma2_c = 1;

	const double  z_ini = -10;
	const double  z_end = 10;
	const int     num_z = 5000;

	const int     N_max = 500; //not used so far

	const int     iterations = 100000;

	const int     seed =  4011; //random seed
	
	//Definitions
	int	     i, j, k, iN, iM, N, M, n, num_samples;

	double   Q_max, Q_max_aver, Q_max_analytic, J_max_aver_th, J_max_aver;
	double   r1, r2, r3, x1, x2, x3, pi, mu, x, sigma2;
	double   sigma_0, sigma_c, mu_chosen, x_old;
	double   z, hz, sum, function_M, coeff;

	pi = 4.0 * atan(1.0);
	sigma_0 = pow(sigma2_0, 0.5);
	sigma_c = pow(sigma2_c, 0.5);

	
	FILE *value_actions, *value_actions_th;

	value_actions = fopen("value_actions_gauss_capacity_fixed2.m", "w");
	value_actions_th = fopen("value_actions_gauss_capacity_fixed_th2.m", "w");
	
	//seed
    TRandomMersenne rg(seed);

	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		for (iM = 1; iM <= N; iM++) {

			//that if, if the division gives an integer:
			if (floor(N / iM) - (1.*N) / (1.*iM) == 0) {

				M = iM;

				num_samples = floor(N / M);

				//initial condition
				Q_max_aver = 0;

				for (i = 0; i < iterations; i++) {

					//we draw mu for every option, and from this we draw x observed, and chose based on x

					x_old = -9999999999;
					for (j = 0; j < M; j++) { //we sample M options

						r1 = rg.Random();
						r2 = rg.Random();
						x1 = pow(-2.*log(r1), 0.5) * cos(2 * pi * r2);
						x2 = pow(-2.*log(r1), 0.5) * sin(2 * pi * r2);

						//Gaussian prior
						//drawing mu for every option
						//mu = mu_0 + sigma_0 * x1;

						//flat prior
						r3 = rg.Random();
						mu = mu_0 * r3;

						//drawing the sample mean x for that option based on a number of samples = num_samples
						x = mu + sigma_c * x2 / sqrt(num_samples);

						//if x is larger than the other option, then we select it, and then we get the corresponding mu
						if (x > x_old) {
							mu_chosen = mu;
							x_old = x;
						}

					}

					Q_max_aver = Q_max_aver + mu_chosen;

				}

				J_max_aver = Q_max_aver / iterations;


				printf("%i %i %i %f \n", N, M, num_samples, J_max_aver);

				fprintf(value_actions, "%i %i %i %f \n", N, M, num_samples, J_max_aver);

			}
		}
	}

	//Theory, Numerical Integration of Analytical expression
	hz = (z_end - z_ini) / num_z;

	//for (iN = 0; iN < N_max; iN++) {

		//N = iN;
	for (iN = 0; iN < num_N; iN++) {

		N = N_vec[iN];

		for (iM = 1; iM <= N; iM++) {

			//that if, if the division gives an integer:
			if (floor(N / iM) - (1.*N) / (1.*iM) == 0) {

				M = iM;

				num_samples = floor(N / M);

				//initial condition
				sum = 0;

				//numerical integration
				for (i = 0; i < num_z; i++) {
					z = z_ini + hz * i;

					sum = sum + z * exp(-pow(z, 2.) / 2.) * pow((1. + erf(z / sqrt(2.))) / 2., M - 1);
				}
				function_M = M * sum * hz / pow(2.*pi, 0.5);

				coeff = sigma_0 / sqrt(1 + sigma2_c / sigma2_0 / num_samples); 
				
				Q_max_analytic = mu_0 + coeff * function_M;

				printf("%i %i %i %f \n", N, M, num_samples, Q_max_analytic);

				fprintf(value_actions_th, "%i %i %i %f \n", N, M, num_samples, Q_max_analytic);
			
			}
	
		}
	
	}

	fclose(value_actions);
	fclose(value_actions_th);
	
	return 0;
}

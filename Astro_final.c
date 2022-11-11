#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Physical and math constant */
#define E 2.7182818
#define PI 3.1415926
#define C 2.99792458e8
#define PLANCK 6.62607004e-34
#define BOLTZ 1.38064852e-23
#define GOLDEN 0.61803

/* Basic properties of observing */
#define WAVE_LEGNTH_B 435.3e-9 
#define WAVE_LEGNTH_V 547.7e-9
#define WAVE_LEGNTH_R 634.9e-9

/* Other constants */
#define N 59 

void star_initialize(void);
void read_source(void);
double vector_component(double wave_up, double wave_down, double temp_var, double ratio);
double built_norm(double temp_var, double B_V_ratio, double V_R_ratio, double R_B_ratio);
double golden_section(double B_V_ratio, double V_R_ratio, double R_B_ratio);	// Solve for Temp.
double solve_radius(double flux_B, double flux_V, double flux_R, double temp);	// Solve for radius.
void list_radius_temp(void);

typedef struct STAR_SOURCE{
	int index;
	double time;
	
	/* Initial data from LOT's */
	double flux_B;
	double flux_V;
	double flux_R;
	
	/* processing: flux ratio vector */
	double B_V_ratio;
	double V_R_ratio;
	double R_B_ratio;
	
	/* Output */
	double radius;
	double radius_rate;
	double temp;
} STAR_SOURCE;

STAR_SOURCE data[N];

int main()
{
	int i = 0;
	
	star_initialize();
	read_source();
	
	data[0].radius_rate = 100;
	
	for (i = 0; i < N; i++)	
	{
		data[i].temp = golden_section(data[i].B_V_ratio, data[i].V_R_ratio, data[i].R_B_ratio);
		data[i].radius = solve_radius(data[i].flux_B, data[i].flux_V, data[i].flux_R, data[i].temp);
		data[i].radius_rate = 100 * (data[i].radius - data[0].radius) / data[0].radius; 	
	}
	
	list_radius_temp();
	
	return 0;
}

void star_initialize(void)
{
	int i = 0;
	
	for (i = 0; i < N; i++)
	{
		data[i].index = 0;
		data[i].time = 0;
		data[i].flux_B = 0;
		data[i].flux_V = 0;
		data[i].flux_R = 0;
		data[i].B_V_ratio = 0;
		data[i].V_R_ratio = 0;
		data[i].R_B_ratio = 0;
		data[i].radius = 0;
		data[i].radius_rate = 0;
		data[i].temp = 0;
	}
}

void read_source(void)
{
	int i = 0;
	
	while (getchar() != '\n');
	
	for (i = 0; i < N; i++)
	{
		data[i].index = i+1;
		scanf("%lf", &data[i].time);
		scanf("%lf", &data[i].flux_B);
		scanf("%lf", &data[i].flux_V);
		scanf("%lf", &data[i].flux_R);
		scanf("%lf", &data[i].B_V_ratio);
		scanf("%lf", &data[i].V_R_ratio);
		scanf("%lf", &data[i].R_B_ratio);
	}
	
}

double vector_component(double wave_up, double wave_down, double temp_var, double ratio)
{
	double coeff = 0;
	double expo = 0;
	double component = 0;
	
	coeff = pow(wave_down, 5) / pow(wave_up, 5);
	expo = (exp(PLANCK * C / (wave_down*BOLTZ*temp_var)) - 1) / (exp(PLANCK * C / (wave_up*BOLTZ*temp_var)) - 1);
	
	component = coeff * expo - ratio;
	
	return component;	
}

double built_norm(double temp_var, double B_V_ratio, double V_R_ratio, double R_B_ratio)
{
	double x = 0, y = 0, z = 0;
	double norm;
	
	x = vector_component(WAVE_LEGNTH_B, WAVE_LEGNTH_V, temp_var, B_V_ratio);
	y = vector_component(WAVE_LEGNTH_V, WAVE_LEGNTH_R, temp_var, V_R_ratio);
	z = vector_component(WAVE_LEGNTH_R, WAVE_LEGNTH_B, temp_var, R_B_ratio);
	
	norm = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	
	return norm;
}

double golden_section(double B_V_ratio, double V_R_ratio, double R_B_ratio)
{
	double lower_bound = 0, upper_bound = 25000;
	double x1 = 0, x2 = 0;
	double err = 1e-4;
	
	x1 = lower_bound + (upper_bound - lower_bound) * (1 - GOLDEN);
	x2 = lower_bound + (upper_bound - lower_bound) * GOLDEN;
	
	while ((x2 - x1) > err)
	{
		if (built_norm(x1, B_V_ratio, V_R_ratio, R_B_ratio) > built_norm(x2, B_V_ratio, V_R_ratio, R_B_ratio))
		{
			lower_bound = x1;
			x1 = x2;
			x2 = lower_bound + (upper_bound - lower_bound) * GOLDEN;
		}
		else
		{
			upper_bound = x2;
			x2 = x1;
			x1 = lower_bound + (upper_bound - lower_bound) * (1 - GOLDEN);
		}
		
	}

	return (x1+x2) / 2;
}

double solve_radius(double flux_B, double flux_V, double flux_R, double temp)
{
	// [0]:B, [1]:V, [2]:R.
	double lamda[3] = {0, 0, 0};
	double coeff[3] = {0, 0, 0};
	double expo[3] = {0, 0, 0};
	double radius[3] = {0, 0, 0};
	int i = 0;
	
	lamda[0] = pow(WAVE_LEGNTH_B , 5);
	lamda[1] = pow(WAVE_LEGNTH_V , 5);
	lamda[2] = pow(WAVE_LEGNTH_R , 5);
	
	
	for (i = 0; i < 3; i++)	
		coeff[i] = 2 * PLANCK * C * C / lamda[i];

	expo[0] = PLANCK * C / (WAVE_LEGNTH_B * BOLTZ * temp);
	expo[1] = PLANCK * C / (WAVE_LEGNTH_V * BOLTZ * temp);
	expo[2] = PLANCK * C / (WAVE_LEGNTH_R * BOLTZ * temp);
	
	for (i = 0; i < 3; i++)
		expo[i] = 1 / (exp(expo[i]) - 1);
	
	radius[0] = sqrt(flux_B / (4 * PI * coeff[0] * expo[0]) );
	radius[1] = sqrt(flux_V / (4 * PI * coeff[1] * expo[1]) );
	radius[2] = sqrt(flux_R / (4 * PI * coeff[2] * expo[2]) );
	
	return (radius[0] + radius[1] + radius[2]) / 3;
}

void list_radius_temp(void)
{
	int i = 0;
	
	printf("###--time---------flux_B-------flux_V-------flux_R-------B/V-----V/R-----R/G-----radius----radius_rate-----temp\n");
	
	for (i = 0; i < N; i++)
	{
		printf("#%2d  ", i);
		printf("%10.4lf", data[i].time);
		printf("%13.4lf", data[i].flux_B);
		printf("%13.4lf", data[i].flux_V);
		printf("%13.4lf", data[i].flux_R);
		printf("%8.4lf", data[i].B_V_ratio);
		printf("%8.4lf", data[i].V_R_ratio);
		printf("%8.4lf", data[i].R_B_ratio);
		printf("%12.4e", data[i].radius);
		printf("%8.2lf%%", data[i].radius_rate);
		printf("%15.4lf", data[i].temp);
		printf("\n");
	}
	
}

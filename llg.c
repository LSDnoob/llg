/**********************************************************************************************************

Date: 11 March 2021
Name: Raunak Prasant

Purpose: To solve the Landau-Lifshitz-Gilbert equation using RK45 method and the huen method.

Inputs: None
Output: Solution of the LLG (values of theta, alpha)
		Requirements for graph: error vs iteration number, step size vs RMS error, alpha vs Switching time
		
*********************************************************************************************************/

#include <conio.h> 
#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 
#include <math.h>

#define gamma -1.76e11
#define hz 7.95e4                  //after substituting the value of u0

float mxr[75000];
float myr[75000];
float mzr[75000];				 // mx, my, mz are the cartesian coordinates for the magnetisation (m) vector, r stands for rk45 method

float mxh[75000];
float myh[75000];
float mzh[75000];                // same as above, except that m stands for midpoint method

float rk45(float theta, float h, float time, float alpha, float phi, float theta_final)
{
 	int i = 0;
 	float tsw;
 	int flag = 0; //to check for the switching point
	do
	{
		//for theta
		
		float c = gamma*alpha*hz/(alpha*alpha + 1);
		float k1 = c*sin(theta*M_PI/180);
		float k2 = c*sin((theta + (1/5)*k1*h)*M_PI/180);
		float k3 = c*sin((theta + (3/40)*k1*h + (9/40)*k2*h)*M_PI/180);
		float k4 = c*sin((theta + (3/10)*k1*h - (9/10)*k2*h + (6/5)*k3*h)*M_PI/180);
		float k5 = c*sin((theta - (11/54)*k1*h + 2.5*k2*h - (70/27)*k3*h + (35/27)*k4*h)*M_PI/180);
		float k6 = c*sin((theta + (1631/55296)*k1*h + (175/512)*k2*h + (575/13824)*k3*h + (44275/110592)*k4*h + (253/4096)*k5*h)*M_PI/180);
	
		float theta_next = theta + ((2825/27648)*k1 + (18575/48384)*k3 + (13525/55296)*k4 + (277/14336)*k5 + 0.25*k6)*h;   // fifth order formula
	
		//for phi
	
		k1 = -1*gamma*hz/(alpha*alpha + 1);
		k2 = k1;
		k3 = k1;
		k4 = k1;
		k5 = k1;
		k6 = k1;                         //d(phi)/dt is constant
		
		float phi_next = phi + ((2825/27648)*k1 + (18575/48384)*k3 + (13525/55296)*k4 + (277/14336)*k5 + 0.25*k6)*h;
		
		if (phi_next >= 360)
			phi_next = phi_next - 360;
	
		//printf("%f      %f       %f\n", theta_next, phi_next, time*1e16);
		
		if((theta <= 90.0) && (flag == 0))      // switching time is the moment the value of theta falls below 90 degrees
		{
			tsw = time;
			flag = 1;
		}
		
		time = time + h;
		theta = theta_next;
		phi = phi_next;
		
		mxr[i] = sin(theta*M_PI/180)*cos(phi*M_PI/180);
		myr[i] = sin(theta*M_PI/180)*sin(phi*M_PI/180);
		mzr[i] = cos(theta*M_PI/180);
		i++; 
		                     
	}while (theta >= theta_final);
	return tsw;
}

int heun(float theta, float h, float time, float alpha, float phi, float theta_final)
{
	int i = 0;
 	
	do
	{
		//for theta
		
		float c = gamma*alpha*hz/(alpha*alpha + 1);
		float theta0_next = theta + c*sin(theta*M_PI/180)*h;
		float theta_next = theta + (c*sin(theta*M_PI/180) + c*sin(theta0_next*M_PI/180))*h/2;  // the predictor - corrector equation
		
		//for phi
		
		float phi0_next = -1*gamma*hz/(alpha*alpha + 1);
		float phi_next = phi + (-1*gamma*hz/(alpha*alpha + 1))*h;
	
		if (phi_next >= 360)
			phi_next = phi_next - 360;
			
		//printf("%f      %f       %f\n", theta_next, phi_next, time*1e16);
		time = time + h;
		theta = theta_next;
		phi = phi_next;
		
		mxh[i] = sin(theta*M_PI/180)*cos(phi*M_PI/180);
		myh[i] = sin(theta*M_PI/180)*sin(phi*M_PI/180);
		mzh[i] = cos(theta*M_PI/180);
		i++;  
	}while (theta >= theta_final);
	return i;
}

int main(int argc, char **argv)
{
	
	float theta_in = 179.0, phi_in = 1.0, theta_final = 1.0, alpha;
	float theta = theta_in, time = 0.0, phi = phi_in;
	float h, d;
	printf("Enter the value for the step size as X * 10^-16 (enter only X): ");
	scanf("%f", &h);
	//h = 0.5e-16;
	h = h*1e-16;
	//alpha = 0.2;
	printf("Enter the value of alpha: ");
	scanf("%f", &alpha);
	float h_user = h;                                     // storing what the user enters as step size for later use
	float tsw = rk45(theta, h, time, alpha, phi, theta_final);
	int number_of_points_h = heun(theta, h, time, alpha, phi, theta_final);
	
	// writing the coordinates of the points into a file (obtained by using the rk45 method)
	
	FILE* fp = fopen("m_vector.txt", "w");
	int i;
	for(i = 0; i < 4*number_of_points_h; i++)
		fprintf(fp, "%f %f %f\n", mxr[i], myr[i], mzr[i]);
	fclose(fp);
	
	// calculation of error : error is the distance between the point obtained by rk45 and heun method, (x, y, z) coordinates, and writing into text file
	
	float error;
	fp = fopen("error.txt", "w");
	for (i = 0; i < number_of_points_h; i++)
	{
		error = sqrt((mxr[4*i] - mxh[i])*(mxr[4*i] - mxh[i]) + (myr[4*i] - myh[i])*(myr[4*i] - myh[i]) + (mzr[4*i] - mzh[i])*(mzr[4*i] - mzh[i]));
		fprintf(fp, "%d %f\n", i+1, error);
	}
	fclose(fp);
	
	// variation of RMS error with step size
	
	float rms = 0.0;
	fp = fopen("rms.txt", "w");
	for (h = 0.5e-16; h <= 10e-16; h = h + 0.1e-16)
	{
		tsw = rk45(theta, h, time, alpha, phi, theta_final);               // tsw to just accept the returning value
		number_of_points_h = heun(theta, h, time, alpha, phi, theta_final);  
		for (i = 0; i < number_of_points_h; i++)
			rms += ((mxr[4*i] - mxh[i])*(mxr[4*i] - mxh[i]) + (myr[4*i] - myh[i])*(myr[4*i] - myh[i]) + (mzr[4*i] - mzh[i])*(mzr[4*i] - mzh[i]));
		rms = sqrt(rms);
		fprintf(fp, "%f %f\n", rms, h*1e16);
	}
	fclose(fp);
	
	h = h_user;
	
	// variation of  switching time with alpha 
	
	fp = fopen("time.txt", "w");
	for (alpha = 0.1; alpha <= 1.0; alpha += 0.05)
	{
		tsw = rk45(theta, h, time, alpha, phi, theta_final); 
		fprintf(fp, "%f %f\n", tsw*1e13, alpha);
	}
	fclose(fp);
	
	printf("\nThe required data for plotting the required graphs and the 3D coordinates of the magnetisation vector have been written into 4 seperate files.");
	
	return 0;
}

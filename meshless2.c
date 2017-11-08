/* 
*Author: Maike Jung
*Date: 30.10.2017

Test meshless membrane model

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define Dt 0.001 // time step
#define sigma 1.0 // excluded volume diameter

#define L (20.0*sigma) // box size
#define L4 (L*L/4.0)
#define N 300 // number of particles

//constants for updating velocities
#define xi 1.0
#define m 1.0
#define a0 ( (1-xi*Dt/(2*m))/(1+xi*Dt/(2*m)) )
#define a1 ( (Dt/m)/(1+(xi*Dt/(2*m))) )


void createConfiguration(double (*pos)[3], double (*vel)[3]){
	/*place particles randomly into the box - considering excluded volume*/
	int i, j;
	double x, y, z;
	double dx, dy, dz, dr;
	bool overlap;
	srand ( time(NULL) );
    for (i=0; i<N; i++){
		do{
			overlap = false;
			x = (double)rand()/(double)(RAND_MAX/L);
			y = (double)rand()/(double)(RAND_MAX/L);
			z = (double)rand()/(double)(RAND_MAX/L);
			for (j=0; j<i; j++){
				dx = pos[j][0]-x;
				dy = pos[j][1]-y;
				dz = pos[j][2]-z;
				// include periodic boundary conditions
				if (dx*dx>=L4) dx = dx - round(dx/L)*L;
				if (dy*dy>=L4) dy = dy - round(dy/L)*L;
				if (dz*dz>=L4) dz = dz - round(dz/L)*L;
				dr = sqrt(dx*dx + dy*dy + dz*dz);
				if (dr < 1.5) overlap = true;
			}
		} while(overlap);
        pos[i][0] = x;
		pos[i][1] = y;
		pos[i][2] = z;
		vel[i][0]=0.01;
		vel[i][1]=0.01;
		vel[i][2]=0.01;
		printf("placed particle %d \n", i);
    }	
}

void read_configuration(double (*pos)[3], double (*vel)[3]){
	/*read in starting configuration*/
	printf("yea\n");
    FILE *myFile;
    myFile = fopen("final_configuration6.dat", "r");
	int i;
    for (i = 0; i < N; i++) {
        fscanf(myFile, "%lf %lf %lf", &pos[i][0], &pos[i][1], &pos[i][2]);
    }
    fclose(myFile);
}

double f_cut(double dr){
	double r_cut = 2.4*sigma;
	double r_att = 1.9*sigma;
	int n = 6;
	double A = log(2)* ( pow( r_cut/r_att, n) - 1 );

	double f;
	if (dr < r_cut){
		f = exp(A*(1 + 1/( pow( dr/r_cut, n) - 1) ) );
	}
	else{
		f = 0.0;
	}
	return f;
}

double df_cut(double dr){
	double r_cut = 2.4*sigma;
	double r_att = 1.9*sigma;
	int n = 6;
	double A = log(2)* ( pow( r_cut/r_att, n) - 1 );

	double f, df;
	if (dr < r_cut){
		f = exp(A*(1 + 1/( pow( dr/r_cut, n) - 1) ) );
	}
	else{
		f = 0.0;
	}
	df = -A*n/r_cut * pow(dr/r_cut,n-1)/pow((pow(dr/r_cut,n) -1) ,2) * f;
	return df;
}


void calcNewVelocities(double (*pos)[3], double (*vel)[3]){

	int i, j, k;
	double dx, dy, dz;
	double f, dr, r_ik, r_hut, U_rep;
	double epsilon = 2.0;
	double rho = 0.0;
	double rhos = 14.0;

	// for each particle component: calculate new velocity (TODO Gaussian noise)
	for (i=0; i<N; i++){ 
		for (j=0; j<3; j++){
			f = 0;

			// cal rho
			for (k = 0; k<N; k++){
				dx = pos[i][0]-pos[k][0];
				dy = pos[i][1]-pos[k][1];
				dz = pos[i][2]-pos[k][2];
				// include periodic boundary conditions
				if (dx*dx>=L4) dx = dx - round(dx/L)*L;
				if (dy*dy>=L4) dy = dy - round(dy/L)*L;
				if (dz*dz>=L4) dz = dz - round(dz/L)*L;
				dr = sqrt(dx*dx + dy*dy + dz*dz);
				if(i!=k) rho += f_cut(dr);
			}
			for (k=0; k<N; k++){ 
				if (i!=k){
					dx = pos[i][0]-pos[k][0];
					dy = pos[i][1]-pos[k][1];
					dz = pos[i][2]-pos[k][2];
					// include periodic boundary conditions
					if (dx*dx>=L4) dx = dx - round(dx/L)*L;
					if (dy*dy>=L4) dy = dy - round(dy/L)*L;
					if (dz*dz>=L4) dz = dz - round(dz/L)*L;
					dr = sqrt(dx*dx + dy*dy + dz*dz);
					//if (dr < 1.0) printf("dr %.10f", dr);
					if (dr <= 2.4*sigma){
						r_ik = pos[i][j]-pos[k][j];
						r_hut = (r_ik - round(r_ik/L)*L )/dr;
						U_rep = exp(-20.0*( (dr/sigma) -1));
						if (U_rep > 10000000) printf("Urep %f \n", U_rep);
						f += (20.0/sigma)*U_rep*r_hut + epsilon*df_cut(dr)*r_hut*epsilon*exp(4*rhos)/( exp(4*rho) + exp(4*rhos) ); 
					}
				}
			}
			rho = 0.0;
			vel[i][j] = a0*vel[i][j] + a1*f;
		}
    }	
}

void calcNewPositions(double (*pos)[3], double (*vel)[3]){
	/* for each particle component: calculate new position*/
	int i, j;

	for (i=0; i<N; i++){ 
		for (j=0; j<3; j++){
			pos[i][j] = pos[i][j] + vel[i][j]*Dt;
		}
    }
	
}

void printConfiguration(char *name, double (*pos)[3]){
	/*create file with current configuration*/
	int i; 
	double x, y, z;
	FILE *f = fopen(name, "w+");
	for(i=0; i<N; i++){
		/*place particles in small box*/
		x = fabs(pos[i][0]);
		y = fabs(pos[i][1]);
		z = fabs(pos[i][2]);
		if (x > L) x = x - floor(x)/1.0;
		if (y > L) y = y - floor(y)/1.0;
		if (z > L) z = z - floor(z)/1.0;
    	fprintf(f, "%f %f %f \n", x, y, z);
	}
    fclose(f);
}

void printVelocity(char *name, double (*vel)[3]){
	/*create file with current particle velocities*/
	int i; 
	double x, y, z;
	FILE *f = fopen(name, "w+");
	for(i=0; i<N; i++){
		x = vel[i][0];
		y = vel[i][1];
		z = vel[i][2];
    	fprintf(f, "%f %f %f \n", x, y, z);
	}
    fclose(f);
}

int main(void){
	printf("START\n");

	double position[N][3];
	double velocity[N][3];
	int i, j;
	char start[] = "starting_configuration_2.dat";
	char mid[] = "intermediate_configuration_2.dat";
	char mid2[] = "intermediate_configuration2_2.dat";
	char mid3[] = "intermediate_configuration3_2.dat";
	char mid4[] = "intermediate_configuration4_2.dat";
	char end[] = "final_configuration_2.dat";
	char endv[] = "final_velocity_2.dat";
	
	/*initialize configuration*/
	createConfiguration(position, velocity);
	printConfiguration(start, position);

	//read_configuration(position, velocity);

	/*calculate new positions and velocities*/
	for(i=0; i<5000000; i++){
		printf("%d \n",i);
		calcNewVelocities(position, velocity);
		calcNewPositions(position, velocity);
		if(i==1000) printConfiguration(mid, position);
		if(i==10000) printConfiguration(mid2, position);
		if(i==500000) printConfiguration(mid3, position);
		if(i==1000000) printConfiguration(mid4, position);
	}
	
	/*store final configuration*/
	printConfiguration(end, position);
	printVelocity(endv, velocity);	

    printf("DONE\n");
    return 0;
}

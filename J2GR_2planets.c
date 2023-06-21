//Two-body Simulation, Accounting for J2 Gravitational Moment
//by Olivia Maynard
//
//
//

#include "./oliviarebound/rebound/rebound.h"
#include "../reboundx/reboundx/reboundx.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

//Currently testing system: TOI-172
//All length units need to be in AU


	//Global variables
	const double Mplanet = 0.005173869; //Solar masses
        const double Rplanet = 0.00046116819; //AU
	
	const double Mstar = 1.128; //Mass, solar mass
	const double Rstar = 0.0082638803; //Radius, AU

	const double J2 = (7.152543E-7); //J2 of the star

        const double tmax = 5000;
	const double tout = 0.025;
	
        //const double J2 = 649.914; //J2 of the star
        //const double J2 = 0.001115;
        const double ObliquityStar = 0; //Obliquities of the star are unknown?

int main(int argc, char* argv[]) {
	
	struct reb_simulation* r = reb_create_simulation();

	double pi  = 4*atan(1);

	FILE *file;
	file = fopen("input.txt", "r");
	int j = 0;
	float array[30] = {0};
	float number = 0;
	while (!feof(file)) {
		fscanf(file, "%f", &number);
		array[j] = number;
		//printf("%.4f\n", array[j]);
		j++;
	}
		
	//Star variables (Body 1)
	double Mstar = array[0]; //Mass, solar mass
        double Rstar = array[1]; //Radius, AU
	
	//Planet 1 variables (Body 2)	
	double Mplanet = array[2]; //Mass, solar mass
        double Rplanet = array[3]; //Radius, AU
	double a = array[4]; //Semimajor axis, AU
        double e = array[5]; //Eccentricity, unitless
        double i = array[6]; //Inclination, radians
	double omega = 0; //Argument of pericenter
	double Omega = pi/2; //Longitude of ascending node 
	double f = 0; //True anomaly

	//Planet 2 variables (Body 3)
	double Mthird = Mplanet/3; //Mass, solar mass
	double Rthird = Rplanet/3; //Radius, AU
	double a_tri = 2*a; //Semimajor axis, AU
	double e_tri = 2*e; //Eccentricity, unitless
	double i_tri = i; //Inclination, radians
	double f_tri = 0; //True anomaly
	double omega_tri = 0; //Argument of pericenter
	double Omega_tri = pi/2;  //Longitude of ascending node 	

	//Declare functions 
	void heartbeat(struct reb_simulation* r); //Prints to the screen
	void force_J2(struct reb_simulation* r); //Gets J2 external force

	//Declare simulation parameters
	r->integrator = REB_INTEGRATOR_IAS15;
	r->dt = 1e-8;
	r->N_active = 3;
	r->G = 4*(pow(pi,2));

	//Adding simulation elements
	reb_add_fmt(r, "m r", Mstar, Rstar); //Star - particles[0]
	reb_add_fmt(r, "m r a e inc Omega omega f", Mplanet, Rplanet, a, e, i, Omega, omega, f); //First planet - particles [1]
	reb_add_fmt(r, "m r a e inc Omega omega f", Mthird, Rthird, a_tri, e_tri, i_tri, Omega_tri, omega_tri, f_tri); //Second planet - particles[2]
	 
	//Setup functions and integrate
	

	//GENERAL RELATIVITY ADDITION
        struct rebx_extras* rebx = rebx_attach(r);
    	// Could also add "gr" or "gr_full" here.  See documentation for details.
    	struct rebx_force* gr = rebx_load_force(rebx, "gr");
    	rebx_add_force(rebx, gr);
    	// Have to set speed of light in right units (set by G & initial conditions), AU/yr
    	rebx_set_param_double(rebx, &gr->ap, "c", 63240.2);

	reb_move_to_com(r); //Move to center of mass
	system("rm -v a.txt"); //Deletes old file
	r->heartbeat = heartbeat;
	r->additional_forces = force_J2;
	reb_integrate(r, tmax);
	rebx_free(rebx);
	reb_free_simulation(r);
}


void force_J2(struct reb_simulation* r){
        if (J2==0) return;

        const struct reb_particle star = r->particles[0]; //Star            // cache
        const int N = r->N; //Number of particles in the system ?
#pragma omp parallel for //Code segment will be executed by multiple threads
        for (int i=1;i<N;i++){
                const struct reb_particle p = r->particles[i];  // cache
                const double sprx = p.x-star.x;
                const double spry = p.y-star.y;
                const double sprz = p.z-star.z;
                const double prx  = sprx*cos(-ObliquityStar) + sprz*sin(-ObliquityStar);
                const double pry  = spry;
                const double prz  =-sprx*sin(-ObliquityStar) + sprz*cos(-ObliquityStar);
                const double pr2  = prx*prx + pry*pry + prz*prz;                // distance^2 relative to planet
                const double fac  = 3.*r->G*J2*star.m*Rstar*Rstar/2./pow(pr2,3.5);

                const double pax  = fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
                const double pay  = fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
                const double paz  = fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);

                r->particles[i].ax += pax*cos(ObliquityStar) + paz*sin(ObliquityStar);
                r->particles[i].ay += pay;
                r->particles[i].az +=-pax*sin(ObliquityStar) + paz*cos(ObliquityStar);
        }
}

void heartbeat(struct reb_simulation* r){
        if(reb_output_check(r, (tmax/100))){                           //Triggers output at regular time intervals.  output something to screen

                reb_output_timing(r, tmax);

        }
        if(reb_output_check(r, tout)){                           // output semimajor axis to file
                FILE* f = fopen("a.txt","a");
                const struct reb_particle star = r->particles[0];
                const int N = r->N;
                for (int i=1;i<N;i++){
                        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],star);
			float x = r->particles[i].x;
			float y = r->particles[i].y;
			float z = r->particles[i].z;
			float dist = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
                        fprintf(f,"%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",r->t,o.a,o.e,o.inc, o.Omega, o.omega, o.f, dist); 
			//1. Current sim time //2. Semimajor axis //3. Eccentricity //4. Inclination
			//5. Inclination of ascending node //6. Argument of pericenter //7. True anomaly
                }
                fclose(f);
        }
}






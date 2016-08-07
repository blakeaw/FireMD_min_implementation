/*
Implementation of the FIRE MD minimizer 
(Bitzek, Koskinen, Gahler, Moseler, Gumbsch, Phys Rev Lett, 97, 170201 (2006).)
with a periodic Lennard-Jones system.

The original MD code was written by Amir Taghavi Nasrabadi (axt128930@utdallas.edu)
and subsequently modified by Blake Wilson (blake.wilson@utdallas.edu) to
implement the FIRE minimization routine 

email:blake.wilson@utdallas.edu
Please e-mail me if you have comments, questions, and suggestions, or to report errors/bugs.

reccomended compile:
g++ -O2 pFireMDmin_v1.cpp

*/

#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <time.h>
// set the number of atoms
#define n_atoms_n 17
// Values for RNG - 64 bit Mersenne Twister
#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */
using namespace std;

// Define Classes
// class Atom
class Atom {
	public:
		double x, y, z;  // position
		double v_x, v_y, v_z; //velocity (angstrom/fs)
		double f_x, f_y, f_z; //force
		double sigma , epsilon ; // sigma (Angstrom) & epsilon (kcal/mol)
		double mass; // mass (kg/mol)
	Atom (): sigma(1.0), epsilon(1.0), mass (1.0) { 
               }
           };
           
//class Frame which contains an array of atoms!
class Frame {
	public:
		int natoms; //number of atoms
          double boxx, boxy, boxz; //box dimensions!
          double rCutOff, rSkin; //values of cutoff and skin!
		double t; //time
		Atom *atom;
		//constructor
		Frame(): natoms(n_atoms_n), boxx(5.0), boxy(5.0), boxz(5.0),rCutOff (2.5) {
 			atom = new Atom[natoms]; //builds an array of atoms 
		}
        //destructor
		~Frame(){
			delete [] atom;
		}
        //deep copy of another frame object - assumes same number of atoms
          void equate (Frame &qqq) {
                 for ( int i = 0; i < natoms; i++)
                 { 
                      atom[i].x = qqq.atom[i].x;
                      atom[i].y = qqq.atom[i].y;
                      atom[i].z = qqq.atom[i].z; 
                      atom[i].v_x = qqq.atom[i].v_x;
                      atom[i].v_y = qqq.atom[i].v_y;	
	                 atom[i].v_z = qqq.atom[i].v_z;	
                      atom[i].f_x = qqq.atom[i].f_x;	
                      atom[i].f_y = qqq.atom[i].f_y;
                      atom[i].f_z = qqq.atom[i].f_z;	
                 }
               }

        //wraps the coordinates back into the box
		void WrapCoordinates(void){

		    for (unsigned int i = 0; i < natoms; ++i)
		    {
			    long double xc = atom[i].x;
			    long double yc = atom[i].y;
			    long double zc = atom[i].z;

			    long double hboxx = boxx/2.0;
			    while(xc > hboxx || xc < -hboxx) {	
				    if(xc > hboxx){
					    xc = xc - boxx;
				    }
				    else if(xc < -hboxx){
					    xc = xc + boxx;
				    }
			    }	


	
			    long double hboxy = boxy/2.0;
			    while(yc > hboxy || yc < -hboxy){
				    if(yc > hboxy){
					    yc = yc - boxy;
				    }
				    else if(yc < -hboxy){
					    yc = yc + boxy;
				    }
			    }


	
			    long double hboxz = boxz/2.0;
			    while(zc > hboxz || zc < -hboxz){
				    if(zc > hboxz){
					    zc = zc - boxz;
				    }
				    else if(zc < -hboxz){
					    zc = zc + boxz;
				    }
			    }

			    atom[i].x = xc;
			    atom[i].y = yc;
			    atom[i].z = zc;
		    }
		
		    return;

        }	
};

// RNG - 64 bit Mersenne Twister
class MTRandomNum {
      public:
	unsigned long long int mt[NN];
        int mti;
	MTRandomNum() : mti(NN+1){}
	// initializes with a seed
        void initialize(unsigned long long int seed){
        	mt[0] = seed;
   		 for (mti=1; mti<NN; mti++){ 
      			  mt[mti] =  (unsigned long long int)(6364136223846793005ULL * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
		}
	
		return;
        }
	// called to generate random number
	// on [0,1] 
	long double generate(void){
		 int i;
   		 unsigned long long x;
   		 static unsigned long long mag01[2]={0ULL, MATRIX_A};

    		if (mti >= NN) { /* generate NN words at one time */	

        	/* if initialize has not been called, */
        	/* a default initial seed is used     */
        		if (mti == NN+1){
        		    initialize(5489ULL);
			} 

        		for (i=0;i<NN-MM;i++) {
        		    x = (mt[i]&UM)|(mt[i+1]&LM);
        		    mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        		}
        		for (;i<NN-1;i++) {
        		    x = (mt[i]&UM)|(mt[i+1]&LM);
        		    mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
        		}
        		x = (mt[NN-1]&UM)|(mt[0]&LM);
        		mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&1ULL)];
	
        		mti = 0;
    		}
  
   		 x = mt[mti++];

    		x ^= (x >> 29) & 0x5555555555555555ULL;
    		x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    		x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    		x ^= (x >> 43);
		// [0,1] real
//		return (x >> 11) * (1.0/9007199244640991.0);
		// [0,1) real
//		return (x >> 11) * (1.0/9007199254740992.0);
		// (0,1) real 
		return ((x >> 12) + 0.5) * (1.0/4503599627370496.0);
	
	}
};

// Prototype Functions
	double PotentialEnergy (Frame &qqq);
	void ForceCalculation (Frame &qqq);
	void verlet (Frame &qqq, Frame &nnn, double timestep);  
     double kinetic (Frame &qqq);
     void pbc (double &dx, double &dy, double &dz, double boxx, double boxy, double boxz);
	long double CalcP(Frame& f);
    void ScaleVelocityFire(Frame& f, long double alpha);
	void ZeroVelocityFire(Frame& f);
// ** Main Function **
int main(int argc, char *argv[])
{
    // Initialization
	   Frame *frame;
        int nsteps = 30000; //number of time steps!
        frame = new Frame[2]; //builds 2 frames! one as current frame and the other one as next frame!
        double timestep = 0.0001; //timestep (femtosecond)
        frame[0].t = 0.0; // sets initial time to zero!
	   double potentialEnergy; // potential energy (kcal/mol)	
       int screenfreq = 5000; // frequency of timesteps to output to screen
        
		//FIRE parms
		long double alphastart=0.1;
		long double finc = 1.1;
		long double fdec = 0.5;
		long double falpha = 0.99;
		int Nmin = 5;
		long double dtmax = 10.0*timestep;
		long double alpha = alphastart;
		int nz=0;
		
		
		MTRandomNum rng;
		unsigned long long seed = time(0);
		rng.initialize(seed);
		
		//generate random coordinates
		for (unsigned int i = 0; i < frame[0].natoms; ++i)
		{
			frame[0].atom[i].x=rng.generate()*frame[0].boxx;
			frame[0].atom[i].y=rng.generate()*frame[0].boxy;
			frame[0].atom[i].z=rng.generate()*frame[0].boxz;
			frame[0].atom[i].v_x=0.0;
			frame[0].atom[i].v_y=0.0;
			frame[0].atom[i].v_z=0.0;
			frame[0].atom[i].f_x=0.0;
			frame[0].atom[i].f_y=0.0;
			frame[0].atom[i].f_z=0.0;
		}
    	
		frame[0].WrapCoordinates();        
        ForceCalculation (frame[0]); //calculates forces for the first frame!
        
 //  ** MD MAIN LOOP **
        for (int i = 0; i < nsteps; i++)
            {
         
             
            // updates the position, velocity and force for the next frame!
			
               verlet (frame[0], frame[1], timestep);
				frame[1].WrapCoordinates();
				if (i>0){
					
				
                  	//F1
              		long double P = CalcP(frame[1]);
				   
				    //F2
				    ScaleVelocityFire(frame[1], alpha);
				    //F3
				
				    if (P>0.0 && nz>Nmin){
					    timestep*=finc;
					    alpha*=falpha;
					    if (timestep>dtmax){
						    timestep=dtmax;
					    }
				    }
				    //F4
				    if (P<=0.0){
					
					    timestep*=fdec;
					    alpha=alphastart;
					    ZeroVelocityFire(frame[1]);
				    }	
							
				    if (P>=0.0){
					    ++nz;					
				    }
				    else{
					    nz=0;
				    }
                    //screen output
                    if ( (i%screenfreq)==0 ){
                     cout<<"P "<<P<<" alpha "<<alpha<<" timestep "<<timestep<<endl;
                       // calculates potential energy!
                     potentialEnergy = PotentialEnergy (frame[1]);
    			     cout<<"step "<<i<<" potential energy "<<potentialEnergy<<endl;
             }   
				}
         
            
          
             //check timestep and break if too small
             if (timestep<1.0e-10){
                    break;
                }

              // turns next frame[1] to new current frame[0] for the next timestep!
              frame[0].equate(frame[1]);
              

            } // **End of MD LOOP**
        
      // output the final minized structure
        ofstream outfile ("coordinates_min.xyz");
        outfile << frame[0].natoms << endl;
        outfile << "MD simulation xyz coordinates" << endl;
          // prints out the position of each atom in the output file!
              for (int l = 0; l< frame[0].natoms; l++)
              { 
                outfile << "LJ " << frame[0].atom[l].x << "  " << frame[0].atom[l].y << "  " << frame[0].atom[l].z << endl; 
              }
   
     
     
return EXIT_SUCCESS;

}

// **Define Functions**

// Potential energy function: returns potential energy!

double PotentialEnergy (Frame &qqq){
   int n;
   n = qqq.natoms;
   double r2, dx, dy, dz, sigma, epsilon, sigma2,rCutOff2, u;

   u = 0.0;
   for (int i = 0; i < n - 1 ; i++)
  	  {
 
             for (int j = i + 1; j < n; j++)
               {  
	            dx = qqq.atom[i].x - qqq.atom[j].x;
                 dy = qqq.atom[i].y - qqq.atom[j].y;
                 dz = qqq.atom[i].z - qqq.atom[j].z;
                 pbc(dx, dy, dz, qqq.boxx, qqq.boxy, qqq.boxz);
    	            r2 = pow (dx , 2.0) + pow (dy, 2.0) + pow (dz , 2.0);
                 rCutOff2 = pow (qqq.rCutOff, 2.0);
              // test cutoff
                 if (r2 < rCutOff2)
                    {
  		     sigma = (qqq.atom[i].sigma + qqq.atom[j].sigma) / 2.0;
    		     sigma2 = sigma * sigma;
   		     epsilon = sqrt (qqq.atom[i].epsilon * qqq.atom[j].epsilon);
   		     u += 4.0 * epsilon * (pow ((sigma2 / r2), 6.0) - pow ((sigma2 / r2), 3.0) );
                    }
              }	
  	}		
return u;
 
}


//Force Calculation function

void ForceCalculation (Frame &qqq) {
	int n;
	n = qqq.natoms;
	double r2, dx, dy, dz, sigma, epsilon, sigma2,rCutOff2, u, Fx, Fy, Fz;


 //Dump old force values and reset to zero

   for(int i = 0; i < n; ++i)
         {
	qqq.atom[i].f_x = 0.0;
	qqq.atom[i].f_y = 0.0;
	qqq.atom[i].f_z = 0.0;
         }

//Calculating the new force values

       for (int i = 0; i < n - 1 ; i++)
            {
 
             for (int j = i + 1; j < n; j++)
                 {  
	            dx = qqq.atom[i].x- qqq.atom[j].x;
                 dy = qqq.atom[i].y- qqq.atom[j].y;
                 dz = qqq.atom[i].z- qqq.atom[j].z;
                pbc(dx, dy, dz, qqq.boxx, qqq.boxy, qqq.boxz);
    	           r2 = pow (dx , 2.0) + pow (dy, 2.0) + pow (dz , 2.0);
                 rCutOff2 = pow (qqq.rCutOff, 2.0);
          // test cutoff 
                 if (r2 < rCutOff2)
                    {
    		      sigma = (qqq.atom[i].sigma + qqq.atom[j].sigma) / 2.0;
   		      sigma2 = sigma * sigma;
   		      epsilon = sqrt (qqq.atom[i].epsilon * qqq.atom[j].epsilon);

  		      u = 4.0 * epsilon * (12.0 * pow ((sigma2 / r2), 6.0) - 6.0 * pow ((sigma2 / r2), 3.0));
  		      Fx = u * dx * (1.0/r2);// * 4.18400e-7;
  		      Fy = u * dy * (1.0/r2);// * 4.18400e-7;
  		      Fz = u * dz * (1.0/r2);// * 4.18400e-7;

  		        qqq.atom[i].f_x +=  Fx;
 	             qqq.atom[j].f_x += -Fx;
 	             qqq.atom[i].f_y +=  Fy;
 	             qqq.atom[j].f_y += -Fy;
 	             qqq.atom[i].f_z +=  Fz;
 	             qqq.atom[j].f_z += -Fz; 
                   } 
                 }	
              }		
}

// Verlet function: integrates equations of motion!

void verlet (Frame &qqq, Frame &nnn, double timestep) {
       
      nnn.t = qqq.t + timestep; // time evolution
   // 1st Step
     for ( int i = 0; i < qqq.natoms; i++)
       {
    	      nnn.atom[i].x = qqq.atom[i].x + qqq.atom[i].v_x * timestep + (pow (timestep, 2.0)/(2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_x;
           nnn.atom[i].y = qqq.atom[i].y + qqq.atom[i].v_y * timestep + (pow (timestep, 2.0)/(2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_y;
           nnn.atom[i].z = qqq.atom[i].z + qqq.atom[i].v_z * timestep + (pow (timestep, 2.0)/(2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_z;
   // Half-step
           nnn.atom[i].v_x = qqq.atom[i].v_x + (timestep / (2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_x;
           nnn.atom[i].v_y = qqq.atom[i].v_y + (timestep / (2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_y; 
           nnn.atom[i].v_z = qqq.atom[i].v_z + (timestep / (2.0 * qqq.atom[i].mass)) * qqq.atom[i].f_z;
     
         }
  // Updates force
     
       ForceCalculation (nnn);

 // Second half-step
     for ( int i = 0; i < qqq.natoms; i++)
         {
  	      nnn.atom[i].v_x = nnn.atom[i].v_x + (timestep/(2.0 * qqq.atom[i].mass)) * nnn.atom[i].f_x;
           nnn.atom[i].v_y = nnn.atom[i].v_y + (timestep/(2.0 * qqq.atom[i].mass)) * nnn.atom[i].f_y;
           nnn.atom[i].v_z = nnn.atom[i].v_z + (timestep/(2.0 * qqq.atom[i].mass)) * nnn.atom[i].f_z;
         } 
}
  
    
// Kinetic energy function: returns the total kinetic energy for each frame

double kinetic (Frame &qqq) {
  double k, v2;
  
  k = 0.0;
  for (int i = 0; i<qqq.natoms; i++)
     { 
       v2 = pow (qqq.atom[i].v_x, 2.0) + pow (qqq.atom[i].v_y, 2.0) + pow (qqq.atom[i].v_z, 2.0);
       k += 0.5 * qqq.atom[i].mass * v2;
     }
  return k * 0.239005736e7; // conversion to (kcal/mol) 
}

void pbc (double &dx, double &dy, double &dz, double boxx, double boxy, double boxz){
	int pbc_int (double);
	int xInt = pbc_int(dx/boxx);
	int yInt = pbc_int(dy/boxy);
	int zInt = pbc_int(dz/boxz);

	dx = dx - boxx * (double)(xInt);
	dy = dy - boxy * (double)(yInt);
	dz = dz - boxz * (double)(zInt);
}

// pbc_int function : rounds a real number to the nearest integer!

int pbc_int (double x) {
  int i;
  double y;
  i = (int)(x);
  y = x-(double)(i);
     if (y > 0.5)  {i++;}
     if (y < -0.5) {i--;}      
     return i;
}

//----- FIRE specific Functions---
// Caluclate the P term
long double CalcP(Frame& f){

	long double P=0.0;
	for (unsigned int i = 0; i < f.natoms; ++i)
	{
		P+=f.atom[i].v_x*f.atom[i].f_x;
		P+=f.atom[i].v_y*f.atom[i].f_y;
		P+=f.atom[i].v_z*f.atom[i].f_z;
	}
	return P;

}
// Scale Velocities according to the FIRE protocol
void ScaleVelocityFire(Frame& f, long double alpha){

		for (unsigned int i = 0; i < f.natoms; ++i)
		{
			long double vx = abs(f.atom[i].v_x);
			long double vy = abs(f.atom[i].v_y);
			long double vz = abs(f.atom[i].v_z);
			long double fx = f.atom[i].f_x;
			long double fy = f.atom[i].f_y;
			long double fz = f.atom[i].f_z;
			long double fm = sqrt(fx*fx+fy*fy+fz*fz);
			long double fhx = fx/fm;
			long double fhy = fy/fm;
			long double fhz = fz/fm;
			f.atom[i].v_x*=(1.0-alpha);
			f.atom[i].v_x+=alpha*fhx*vx;
			f.atom[i].v_y*=(1.0-alpha);
			f.atom[i].v_y+=alpha*fhy*vy;
			f.atom[i].v_z*=(1.0-alpha);
			f.atom[i].v_z+=alpha*fhz*vz;
			
		}

}
//Set all velocities to zero
void ZeroVelocityFire(Frame& f){
	for (unsigned int i = 0; i < f.natoms; ++i)
	{
		f.atom[i].v_x=0.0;
		f.atom[i].v_y=0.0;	
		f.atom[i].v_z=0.0;
	}

}

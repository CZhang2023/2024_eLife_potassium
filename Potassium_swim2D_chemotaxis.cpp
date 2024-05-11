//******* surface chemotaxis*********

/* model description 
The model of swim chemotaxis based on Tu et al. [1] 

*/

#include "MersenneTwister.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

const double v0 = 25.0;           // swimming speed um/s
const double pi = 3.14159;
const double Dr = 0.062; // rotation diffusion coefficient. (rad^2/s)
const double d = 1500.0; // size of space 300*300 um^2  x-axis: solid wall ; y-axis: periodic boundary.
const int N = 10000; // cell number

const double D = 1333.3; // in um^2/s. Diffusion constant of potassium ions in water. 19.7×10^(-6) cm^2 s^(-1).  Ref.49 in Prindle15-nature.   80000 um^2/min Fell & Hutchison 1971  from Ref. 2017 cell
//const double T = 7200.0; // in seconds. peroid of fluctuation of potassium concentration. 1-4 hrs fig. 1E in Ref. Humphries2017-cell.
//const double L0 = 50.0; // amplitude of potassium concentration in x = d.  (mM)
//const double w = 2*pi/T;

const int NofS = 1;

// Chemotaxis related parameters.       SPEC model, Tu, Plos conputational biology,2010.
const double kr = 0.005;  // 0.005/s for MeAsp
const double kTumble = 5.0; // switch rate of Run to Tumble[average tumble period=0.2 s]

int main() {
	
	MTRand randtheta,randRunOrTumble,randDr,randTumbleAngle,randpos;
	
	double a, yp,bias,kRun,lig;  // switch rate of Run to Tumble 1.3636
	double dt,t,tsave;       // time step
	
    double x[N],y[N],theta[N],m[N],Run[N];
	double Tp[2] = {80.0,1200.0};
	double L0p[1] = {1.0};
	double T = 0.0;
	double L0 = 0.0;
	
	for (int TI=0;TI<2;TI++){
		T = Tp[TI];
		double w = 2*pi/T;
		
		for (int LI=0;LI<1;LI++){
			L0 = L0p[LI];
			char s[200]; //sprintf(s, "temp_index%04d.txt",Itime+1);
			
			for(int Itime=0;Itime<NofS;Itime++) {
			
				cout<<"Simulating L0 = "<<L0<<"/ D = "<<D<<"/ T = "<<T<<"/ "<<Itime+1<<'\\'<<NofS<<" ..."<<endl; 
				sprintf(s, "temp_D%.0f_L0%.04f_T%.0f_%03d.txt",D,L0,T,Itime+1);
	
	
				 // ************ initialized cell states ****************
			   for(int i=0;i<N;i++) {
			   		x[i]=d*randpos.rand();
			   		y[i]=d*randpos.rand();
			   		theta[i] = pi*(2*randtheta.randExc()-1);
			   		lig = L0*(1-cos(-sqrt(w/2/D)*x[i])*exp(-sqrt(w/2/D)*x[i])); // t=0;
			   		m[i] = 1+(log((1+lig/0.05)/(1+lig/10000.0)))/1.7-(log(2))/1.7/0.35;
					Run[i] = 1;   
			   } 
			   
			
				ofstream outfile("Simulation.txt");
				dt = 0.01; 
				t = 0.0;
				tsave = 0.0;
				while(t<=1000+4*T) {
					
					
					if(tsave>=1.0 or t==0){
					
						double mx = 0.0;
						double my = 0.0;
						double mm = 0.0;
						double mcellnum = 0.0;
						for(int i=0;i<N;i++) {
							mx = mx+x[i];
							my = my+y[i];
							mm = mm+m[i];
							if (x[i]<=100){ //100 um以内的计数。 
								mcellnum = mcellnum+1;
							}
						}
						mx = mx/N;
						my = my/N;
						mm = mm/N;
						
						outfile<<t<<' '<<mx<<' '<<my<<' '<<mm<<' '<<mcellnum<<endl;
						tsave = 0;
					}
					
					
					tsave+=dt;
				 	t+=dt;
					
					for(int i=0;i<N;i++) {
						// 确定此时第i个细菌的运动状态是run还是tumble （Run = 1 or 0）;
						lig=L0*(1-cos(w*t-sqrt(w/2/D)*x[i])*exp(-sqrt(w/2/D)*x[i]));
						//lig = 50.0*exp(-x[i]/d);
						a=1.0/(1.0+exp(0.35*(1.7*(1.0-m[i])+log((1+lig/0.05)/(1+lig/10000.0)))));  //参数来自于MWCfit 
						yp=7.86*a; //yp_nostimuli=2.62uM -> bias=0.15.
						bias=pow(7.86*a,10.3)/(pow(7.86*a,10.3)+pow(3.1,10.3));
						
						kRun = bias/0.11;
						if(Run[i]==1){
							if(kRun*dt>=randRunOrTumble.randDblExc()){
								Run[i] = 0;
							}
						}	
						else{
							if(kTumble*dt>=randRunOrTumble.randDblExc()){
							  	Run[i] = 1;
							}
						}
						
						
						m[i] += (kr*(1-a)-2*kr*a)*dt;
						
						
						if(m[i]>4.0){
							m[i] = 4.0;
						}
						
						if(m[i]<0.0){
							m[i] = 0.0;
						}
						
						
						// 通过以上的所有影响确定当前细菌的位置与速度方向状态 
						//均匀分布的tumble角度。
						//theta[i] += Run[i]*(sqrt(2*Dr*dt)*randDr.randNorm(0,1))+(1.0-Run[i])*pi*(1.0-2.0*randTumbleAngle.randExc());
						
						//细菌的Tumble角度满足分布 0.5*(1+cos(theta))*sin(theta)   Ref. [2-4]
						theta[i] += Run[i]*(sqrt(2*Dr*dt)*randDr.randNorm(0,1))+(1.0-Run[i])*acos(2*sqrt(1-randTumbleAngle.rand())-1);
						
						x[i] += Run[i]*v0*cos(theta[i])*dt;
						y[i] += Run[i]*v0*sin(theta[i])*dt;
				  	  	
				 		
						if(x[i]<0){
							x[i]=0;
						}
						
						if(x[i]>d){
							x[i]=d;
						}
						
						if(y[i]>d){
							y[i]=y[i]-d;
						}	
						
						if(y[i]<0){
							y[i]=y[i]+d;
						}
					}
				}	
				outfile.close();
				rename("Simulation.txt",s);
			}
		}
	}
	return 0;
}

/* References:
[1] Jiang, L., Ouyang, Q. & Tu, Y. Quantitative modeling of Escherichia coli chemotactic motion in environments varying in space and time. PLoS Comput Biol 6, e1000735 (2010).
[2] Berg, H. C. & Brown, D. A. Chemotaxis in Escherichia coli analysed by three-dimensional tracking. Nature 239, 500-504 (1972).
[3] Chen KC, Ford RM, Cummings PT (1998) The global turning probability density function for motile bacteria and its applications. J Theor Biol 195: 139–155.
[4] Vladimirov N, Lebiedz D, Sourjik V. Predicted auxiliary navigation mechanism of peritrichously flagellated chemotactic bacteria[J]. PLoS computational biology, 2010, 6(3): e1000717.

*/

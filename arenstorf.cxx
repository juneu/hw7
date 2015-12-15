// Homework 7
// RK 4 and RK 5 with step size control
//--------------------
#include <cmath>
#include <iostream>
#include <fstream>
//--------------------

//Vorlage fürs Programm von der Lösung von hw6
void f(double* y0, const double t);
void RKstep(double* yn, const double x, const double dx, const int dim, const double* const b);
void stepsize( double* const, double* const, double&, const int);
//--------------------
using namespace std;
//--------------------
int main(void)
{
        
	ofstream out("solution.dat");
        const int dim = 4; 
        const double L=17.065216560157;
	double dt = 0.001,t=0;
	
       double y4[dim]; // values x1, x2, y1, y2 calculated with RK 4
       double y5[dim]; // values x1, x2, y1, y2 calculated with RK 5
       
    
       int number=7; // how many b variables in butcher tableau?
       double b4[number], b5[number]; // b values from butcher tableau for RK 4 and RK5, respectively
       
       b4[0]=5179.0/57600;
       b4[1]=0;
       b4[2]=7571.0/16695;
       b4[3]=393.0/640;
       b4[4]=-92097.0/339200;
       b4[5]=187.0/2100;
       b4[6]=1.0/40;
       
       b5[0]=35.0/384;
       b5[1]=0;
       b5[2]=500.0/1113;
       b5[3]=125.0/192;
       b5[4]=-2187.0/6784;
       b5[5]=11.0/84;
       b5[6]=0;

	
	y4[0]=0.994; y4[1]=0; y4[2]=0; y4[3]=-20015851063790.8/10000000000000; //Initial conditions for RK4
	
                        
                        //Wieso kann ich nicht einfach die Kommazahl mit mehr als vier Nachkommastellen hier eingeben???
	
        out << t << "\t" << y4[0] << "\t" << y4[2] << "\t" << endl; //print the initial values in file
        
	while(t<=L)
	{
                    for(int i=0; i<dim; i++) y5[i]=y4[i]; //further calculation with results from RK4
                    t += dt;
                
		RKstep(y4, t, dt, dim, b4); // calculate value (first ks and then function f) with RK 4
                RKstep(y5, t, dt, dim, b5); // same for RK5
                
                stepsize(y4, y5, dt, dim);
            
                for(int i=0; i<dim; i++) y5[i]=y4[i]; //further calculation with results from RK4
                
		out << t << "\t" << y4[0] << "\t" << y4[2] << "\t" <<dt<< endl;
	}
	
	out.close();
	return(0);
}

//-------------------

void stepsize( double* const y4,  double* const y5, double& dt, const int dim){
    
    const double tol=0.00001; //tolerance local error
    const double q=0.5; //safety factor
    double error,errormax=0;
    for(int i=0; i<dim; i++){
    error = abs(y4[i]-y5[i]);	//calculate local error
    if (error>errormax){        //maximum local error
      errormax = error;}
       }
        dt*=pow(tol/errormax,0.2);// 1/5 oder 1/6 ???? safety factor?
}

//-------------------

void RKstep(double* yn, const double x, const double dx, const int dim, const double* const b)
{
	double k1[dim], k2[dim], k3[dim], k4[dim], k5[dim], k6[dim], k7[dim];
	double y0[dim];
        
        for(int i=0; i<dim; i++) y0[i]=yn[i]; //save old values in y0 to calculate k correctly

  for(int i=0;i<dim; i++) k1[i] = y0[i];
	f(k1, x);

	for(int i=0;i<dim; i++) k2[i] = y0[i] + dx* 1.0/5 * k1[i];
  f(k2, x+ (1.0/5 *dx));

	for(int i=0;i<dim; i++) k3[i] = y0[i] + dx * ( (3.0/40 *k1[i]) + (9.0/40 *k2[i]) );
	f(k3, x+ (3.0/10 *dx));

  for(int i=0;i<dim; i++) k4[i] = y0[i] + dx * ( (44.0/45 * k1[i]) + (-56.0/15 *k2[i]) + (32.0/9 *k3[i]));
	f(k4,  x+ (4.0/5 *dx));
  
  for(int i=0;i<dim; i++) k5[i] = y0[i] + dx * ( (19372.0/6561 * k1[i]) + (-25360.0/2187 *k2[i]) + (64448.0/6561 *k3[i]) + (-212.0/729 *k4[i]));
	f(k5,  x+ (8.0/9 *dx));

  for(int i=0;i<dim; i++) k6[i] = y0[i] + dx * ( (9017.0/3168 * k1[i]) + (-355.0/33 *k2[i]) + (46732.0/5247 *k3[i]) + (49.0/176 *k4[i]) + (-5103.0/18656 *k5[i]));
	f(k6,  x+ dx);

  for(int i=0;i<dim; i++) k7[i] = y0[i] + dx * ( (35.0/384 * k1[i]) + (500.0/1113 *k3[i])  + (125.0/192 *k4[i]) + (-2187.0/6784 *k5[i]) + + (11.0/84 *k5[i]));
	f(k7,  x+ dx);


	for(int i=0;i<dim; i++)
	 yn[i] = y0[i] + dx*(k1[i]*b[0] + k2[i]*b[1] + k3[i]*b[2] + k4[i]*b[3] + k5[i]*b[4] + k6[i]*b[5] + k7[i]*b[6]);
	  
}

//-------------------
// 
void f(double* y0, const double t) //hier unabhängig von t also Übergabe nicht notwendig
{
	const int dim=4;
	double y[dim] = { y0[0], y0[1], y0[2], y0[3]}; //save y-values that were calculated with k in y (sonst nimmt man bei der Berechnung im nächsten Schritt die falschen Werte)
        double r, s;
        const double m=0.012277471; 
        r=sqrt((y[0]+m)*(y[0]+m)+(y[2]*y[2]));
        s=sqrt((y[0]-1+m)*(y[0]-1+m)+(y[2]*y[2]));
        
  y0[0] = y[1];
  y0[1] = y[0]+2*y[3]-(((1-m)*(y[0]+m)/(r*r*r))) - ((m*(y[0]-1+m))/(s*s*s));
  y0[2] = y[3];
  y0[3] = y[2]-2*y[1]-(((1-m)*y[2])/r*r*r) - ((m*y[2])/(s*s*s));
	
}
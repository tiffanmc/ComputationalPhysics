#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include "time.h"
#include <armadillo>
using namespace arma;

using namespace std;
int main(int argc, char *argv[])
// PROJECT 1
// Exersize b)
{
    ofstream outputFile;
    outputFile.open("v_1000.txt");
    outputFile << setiosflags(ios::showpoint | ios::uppercase);

    clock_t start, finish;        //calculating algorithm time for 1c

//DEFINE all variables

    int N = 10;
    double *a = new double[N+2];    // Last element is N+1
    double *b = new double[N+2];
    double *c = new double[N+2];
    double *f = new double[N+2];
    double *u = new double[N+2];
    double *x = new double[N+2];
    double h = 1.0/(N+1);
    double *ux= new double[N+2];
    double time= time;
    double *E= new double[N+2];


  for (int i=0; i<N+2; i++){
      a[i]=2;                   //define variables again for all i
      b[i]=-1;
      c[i]=-1;
      f[i]=h*h*100*exp(-10*i*h);
      x[i]=i*h;
      ux[i]=1-(1-exp(-10))*i*h-exp(-10*i*h);
  }
  //general form of algorithm
  for (int i=2; i<N+2; i++){
      a[i] = a[i] - b[i-1]*c[i-1]/a[i-1];       //3 flops
      f[i] = f[i] - (f[i-1]*c[i-1])/a[i-1];         //3 flops

  }

  //definition of u from equation Ax=b
  u[N]= f[N]/a[N]; //1 flop

  //finding U (solution)
  for (int i=N-1; i>=1; i--){                       //backwards substitution
      u[i]=(f[i]-b[i]*u[i+1])/a[i];//3 flops



      E[i]=log10(abs(u[i]-ux[i])/ux[i]);  //fine the error of the exact value 1d

      //cout <<E[i] << endl;

  }
  for(int i = 0; i < N+2; i++){
     outputFile << setprecision(10) << setw(20) << u[i] << endl;
  }


  //QUESTION 1D


  //Printf ("CPU time to computer algorithm: %d\n",finish-start);
  //143 clicks (0.143000 seconds) is the conversion rate
  //1154615 CPU for N=10^6
  //132 CPU for N=1000
  //40 CPU for N=100
  //19 CPU for N=10

  //QUESTION 1E

//error

 //N=10=-1.1797
 //N=100=-3.08804
 //N=1000=-5.08005
 //N=10000=-7.07929
 //N=100000=(-8.88828-9.0049)
 //N=10000000=(-6.77136-6.54196)
 //N=10^7=-5.1708 - 5.96231

  // LU decomposition

  start=clock();

  mat A = zeros<mat>(N,N);
  vec b1 = randu<vec>(N);


          A.diag().fill(2);
          A.diag(1).fill(-1);
          A.diag(-1).fill(-1);
  start=clock();
   // A.print("A =");
   // b1.print("b=");
  // solve Ax = b
    vec x1 = solve(A,b1);
    // print x
   // x1.print("x=");
         // find LU decomp of A, if needed, P is the permutation matrix
   mat L, U;
   lu(L,U,A);
   //print l
   //L.print(" L= ");
  // print U
   //U.print(" U= ");
  //Check that A = LU
   //(A-L*U).print("Test of LU decomposition");

 finish=clock();
 time= ((finish-start)/CLOCKS_PER_SEC) ;

 cout << finish-start << endl;
  return 0;

  //n=10 (157 CPU)
  //n=100 (1759)
  //n=1000 (161293)
  //n=10^6 (it doesnt work)
}

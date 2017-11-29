#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

int main (int argc, char *argv[])
{
    FILE    *fout;
    double  r, r_m, r_cut, x, B;
    double diffR;
    double T_r, Tprime_r;

    if (argc != 5){
        printf("Usage: %s B r_cut r_m ofile\n", argv[0]);
	printf("where B has 1/nm units in equation: V(r)=A*exp(-B*r)-C/r^6\n");
	printf("r_cut, r_m, are in nm (r_cut > r_m), used in Tapering function.\n");
	printf("Tapering function f(r)=1.0                   for r < r_m\n");
	printf("                  f(r)=(1.0-x)^3*(1+3x+6x^2) for r_m < r < r_cut\n");
	printf("                  f(r)=0.0                   for r >= r_cut\n");
	printf("               where x=(r-r_m)/(r_cut-r_m)\n");	
	printf("ofile name of the outputfile. E.g: table.xvg\n");
	return(0);
    }
    else
    {
      B = atof(argv[1]);
      r_cut = atof(argv[2]);
      r_m = atof(argv[3]);
      printf("B=%g, r_cut=%g, r_m=%g\n",B,r_cut,r_m);
      diffR = r_cut-r_m;
    }

    fout = fopen(argv[4], "w");
    fprintf(fout, "#\n# Tapered Buckingham Potential\n#\n");
    fprintf(fout, "#\n# B = %g,  r_cut = %g,  r_m = %g \n#\n", B, r_cut, r_m);

    for (r=0; r<=(r_cut+2.0); r+=0.002) {

        double f = 1/r;
        double fprime = 1/(pow(r,2));

        double g = -1/(pow(r,6));
        double gprime = -6/(pow(r,7));

        double h = exp(-B*r);
        double hprime = B*h;
	
	if (r > r_m && r < r_cut){
	    double x = (r-r_m)/diffR;
	    T_r = pow((1-x),3)*(1+3*x+6*x*x);
	    Tprime_r = 30*x*x*(1-x)*(1-x)*r/diffR;

	    //fprintf(stdout,"%g %g %g %g\n",diffR, x,T_r,Tprime_r);
	    gprime = gprime*T_r + g*Tprime_r;
	    hprime = hprime*T_r + h*Tprime_r;
	    // important: modify g,h only after gprime,hprime are modified
	    g = g * T_r;
	    h = h * T_r;
	}

	
	/* print output */
        if (r<0.04 || r >=r_cut) {
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,0.0,0.0,0.0,0.0,0.0,0.0);
        } else {
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,f,fprime,g,gprime,h,hprime);
        }
    }

    fclose(fout);
    return(0);
}

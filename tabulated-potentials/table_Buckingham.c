#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

int main (int argc, char *argv[])
{
    FILE    *fout;
    double  r, r_cut, B;

    if (argc != 4){
        printf("Usage: %s B ofile\n", argv[0]);
	printf("where B has 1/nm units in equation: V(r)=A*exp(-B*r)-C/r^6\n");
	printf("ofile name of the outputfile. E.g: table.xvg\n");
	return(0);
    }
    else
    {
      B = atof(argv[1]);
      r_cut = atof(argv[2]);
      printf("B=%g, r_cut=%g\n",B,r_cut);
    }

    fout = fopen(argv[3], "w");
    fprintf(fout, "#\n# Buckingham Potential\n#\n");
    fprintf(fout, "#\n# B = %g,  r_cut = %g\n#\n", B, r_cut);

    for (r=0; r<=(r_cut+2.0); r+=0.002) {

        double f = 1/r;
        double fprime = 1/(pow(r,2));

        double g = -1/(pow(r,6));
        double gprime = -6/(pow(r,7));

        double h = exp(-B*r);
        double hprime = B*h;

        /* print output */
        if (r<0.04) {
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,0.0,0.0,0.0,0.0,0.0,0.0);
        } else {
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,f,fprime,g,gprime,h,hprime);
        }
    }

    fclose(fout);
    return(0);
}

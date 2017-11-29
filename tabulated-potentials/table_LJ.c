#include <stdio.h>
#include <stdlib.h> 
#include <math.h>

int main (int argc, char *argv[])
{
    FILE    *fout;
    double  r, r_cut;

    if (argc != 3){
        printf("Usage: %s r_cut ofile\n", argv[0]);
	printf("r_cut in nm.\n");
	printf("ofile name of the outputfile. E.g: table.xvg\n");
	return(0);
    }
    else
    {
      r_cut = atof(argv[1]);
      printf("r_cut=%g \n",r_cut);
    }

    fout = fopen(argv[2], "w");
    //fout = fopen("table_LJ6_12.xvg", "w");
    fprintf(fout, "#\n# LJ 6-12 Potential\n#\n");
    fprintf(fout, "#\n# r_cut = %g \n#\n", r_cut);

    for (r=0; r<=(r_cut+2.0); r+=0.002) {

        double f = 1/r;
        double fprime = 1/(pow(r,2));

        double g = -1/(pow(r,6));
        double gprime = -6/(pow(r,7));

        double h = 1/(pow(r,12));
        double hprime = 12/(pow(r,13));

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

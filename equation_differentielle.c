#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "control_fonction.h"


void scanFonction(float y[], float t[], int n, int m)
{
    int i;
    for(i=0; i<(n+1); i++)
    {
        printf("Entrer le %de coefficient de y: ", (i+1));
        scanf("%f", &y[i]);

        printf("\n");
    }

    for(i=0; i<(m+1); i++)
    {
        printf("Entrer le %de coefficient de t: ", (i+1));
        scanf("%f", &t[i]);

        printf("\n");
    }

    /*Affichage*/
    printf("La fonction entrée est y'(t) = ");
    for(i=0; i<(n+1); i++)
    {
        printf("%fy^%d", y[i], i);
        if(i == (n-1) && t[0] > 0)
            printf(" + ");
        else if(i == (n-1) && t[0] < 0)
            printf("  ");
        else if(i != (n-1))
            printf(" + ");
    }
    for(i=0; i<(m+1); i++)
    {
        printf("%ft^%d", t[i], i);
        if(i != m)
            printf(" + ");
    }



}

float calcul_fonction(float t0, float y0)
{
    /*int i, k;
    float result = 0;

    for(i=0, k=0; i<(n+1); i++, k++)
    {
        result += y[i]*pow(y0,k);
        k--;

    }

    for(i=0, k=0; i<(m+1); i++, k++)
    {
        result += t[i]*pow(t0,k);
        k--;

    }*/

    return -y0 + t0 + 1;



}

void euler(float t0, float y0, double h, int iteration)
{
    float yn, tn;
    int i;
    printf("A t%d = %f, on a y%d = %f\n",0, t0, 0, y0);
    for(i=0; i<iteration; i++)
    {
        yn = y0 + h * calcul_fonction(t0, y0);
        tn = t0 + h;
        printf("A t%d = %f, on a y%d = %f\n",(i+1), tn, (i+1), yn);
        t0 = tn;
        y0 = yn;

    }
}

void runge_kutta_milieu(float t0, float y0, double h, int iteration)
{
    float yn, tn, yint;
    int i;
    printf("A t%d = %f, on a y%d = %f\n",0, t0, 0, y0);
    for(i=0; i<iteration; i++)
    {
        yint = y0 + (h/2) * calcul_fonction(t0, y0);
        yn = y0 + h * calcul_fonction((t0 + h/2), yint);
        tn = t0 + h;
        printf("A t%d = %f, on a y%d = %f\n",(i+1), tn, (i+1), yn);
        t0 = tn;
        y0 = yn;



    }

}

void runge_kutta_euler(float t0, float y0, double h, int iteration)
{
    float yn, tn, yint;
    int i;
    printf("A t%d = %f, on a y%d = %f\n",0, t0, 0, y0);
    for(i=0; i<iteration; i++)
    {
        yint = y0 + h * calcul_fonction(t0, y0);
        yn = y0 + (h/2) * (calcul_fonction(t0, y0) + calcul_fonction(t0 + h, yint));
        tn = t0 + h;
        printf("A t%d = %f, on a y%d = %f\n",(i+1), tn, (i+1), yn);
        t0 = tn;
        y0 = yn;
    }

}

void recopie(double hfinal, double t0final, double y0final, int iterationfinal, double h, double t0, double y0, int iteration)
{
    hfinal = h;

    t0final = t0;

    y0final = y0;

    iterationfinal = iteration;
}

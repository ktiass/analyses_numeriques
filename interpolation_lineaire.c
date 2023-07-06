#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "equation_lineaire.h"

void decomposition_choleski_interpo(float ptr[][100], float ptr2[][100], size_t N, size_t M, size_t P, size_t Q)
{
    float L[100][100];
    float U[100][100];
    int i, k ,j;
    float w, z;

    L[0][0] = sqrt(ptr[0][0]);

    for(i=1; i<N; i++)
    {
        L[i][0] = ptr[i][0] / L[0][0];
    }

    for(k=1; k<N; k++)
    {
        for(j=0, w=0; j<k; j++)
        {
            w += pow(L[k][j], 2);
        }
        L[k][k] = sqrt(ptr[k][k] - w);

        for(i=(k+1); i<N; i++)
        {
            for(j=0, z=0; j<k; j++)
            {
                z += L[i][j] * L[k][j];
            }
            L[i][k] = (ptr[i][k] - z)/(L[k][k]);
        }
    }


    /*Tansposée de L, constitution de la matrice U*/
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            U[j][i] = L[i][j];
        }
    }


    /*Descente triangulaire*/
    float t, Tab[100];
    Tab[0]=(ptr2[0][0])/(L[0][0]);
    for(i=1; i<P; i++)
    {
        for(k=0, t=0; k<i; k++)
        {
           t += (L[i][k]) * Tab[k];
        }
        Tab[i] = ((ptr2[i][0]) - t) / (L[i][i]);
    }

    /*printf("\nLes solutions sont: \n");
    for(i=0; i<P; i++)
    {
        printf("y%d= %f\n", (i+1), Tab[i]);
    }*/


    /*Remontée triangulaire*/
    float Tab1[100];
    Tab1[(P-1)] = (Tab[(P-1)])/(U[(N-1)][(M-1)]);
    for(i=(P-2); i>-1; i--)
    {
        for(k=(i+1), t=0; k<(P); k++)
        {
           t += (U[i][k]) * Tab1[k];
        }
        Tab1[i] = (Tab[i] - t) / (U[i][i]);
    }

    printf("\nLe polynôme d'interpolation est :");
    for(int m =(Q-1), i=0 ; i<P; i++)
    {
            printf("%fx^%d", Tab1[i], m);
            if(i!=(P-1))
                printf("+");
            m--;
    }
}

void interpolation_lagrange(float *x, float *y, int nbPt)
{
    //char poly[1000] = "";
    //char chaine[100];
    float c;
    printf("\nLe polynôme d'interpolation est: ");
    for(int i=0; i<nbPt; i++)
    {
        c=1;
        for(int j=0; j<nbPt; j++)
        {
            if(j != i)
            {
                if(x[j] < 0)
                {
                    printf("(x+%.9f)",(-1*x[j]));
                }
                else
                {
                    printf("(x-%.9f)", x[j]);
                }

                c *= (x[i] - x[j]);

            }


        }
        c = y[i]/c;
        printf("*%f", c);
        if(i != (nbPt-1))
        {
            printf("+");
        }




    }

}

void interpolation_newton(float *x, float *y, int nbPt)
{
    float A[100][100];
    int i, j;

    for(i=0; i<nbPt; i++)
    {
        A[i][0] = y[i];
    }

    for(i=1; i<nbPt; i++)
    {
        for(j=1; j<=i; j++)
        {
                A[i][j] = ((A[i][j-1] - A[i-1][j-1]) / (x[i] - x[i-j]));
        }
    }

    printf("\nLe polynôme d'interpolation est: ");
    for(i=0; i<nbPt; i++)
    {
        printf("%f", A[i][i]);
            for(j=0; j<i; j++)
            {
                printf("(x-%f)", x[j]);
            }
        if(i != (nbPt-1))
            printf("+");
    }



}

void matrice(float ptr[][100], size_t N, size_t M)
{
    size_t i,j;
    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            printf("%f\t", (ptr[i][j]));
        }
        printf("\n");
    }
}

void scan_interpo_point(float ptr[], float ptr2[], int N)
{
    size_t i;


    for(i = 0; i < N; i++)
    {
            printf("Entrez le point %d\n", (i+1));
            printf("x%d: ", (i+1));
            scanf("%f", &ptr[i]);
            printf("y%d: ", (i+1));
            scanf("%f", &ptr2[i]);

    }
}

void interpolation_moindres_carree(float *x, float *y, int nbPt)
{
    float P[100][100], B[100][100];
    int i, j, k;
    int p = nbPt - 2;

    for(i=0; i<(p+1); i++)
    {
        for(j=0; j<(p+1); j++)
        {
            for(k=0, P[i][j] = 0; k<(nbPt); k++)
            {
                P[i][j] += pow(x[k],(2*p-i-j));
            }
        }
    }

    for(i=0; i<(p+1); i++)
    {
        for(k=0, B[i][0] = 0; k<(nbPt); k++)
        {
            B[i][0] += (pow(x[k], p-i))*(y[k]);
        }
    }

    decomposition_choleski_interpo(P, B, (p+1), (p+1), (p+1), (p+1));


}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <conio.h>
#include <stdbool.h>
#include <ctype.h>
#include <locale.h>
#include "save_me.h"






void scanMatrice(float * ptr, size_t N, size_t M)
{
    size_t i, j;
    float k;

    for(i = 0; i < N; i++)
    {
        for(j = 0; j < M; j++)
        {
            printf("Taper l'element [%d;%d]: ", (i+1),(j+1));
            scanf("%f", &k);
            ptr[M * i + j] = k;
        }
    }
}

void printMatrice(float * ptr, size_t N, size_t M)
{
    size_t i,j;
    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            printf("%f\t", (ptr[M * i + j]));
        }
        printf("\n");
    }
}

void descente_triangulaire(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    int i,k;
    float t;
    float Tab[100];

    Tab[0]=(ptr2[Q * 0 + 0])/(ptr[M * 0 + 0]);
    for(i=1; i<P; i++)
    {
        for(k=0, t=0; k<i; k++)
        {
           t += (ptr[M * i + k]) * Tab[k];
        }
        Tab[i] = ((ptr2[Q * i + 0]) - t) / (ptr[M * i + i]);
    }

    printf("\nLes solutions sont: \n");
    for(i=0; i<P; i++)
    {
        printf("x%d= %f\n", (i+1), Tab[i]);
    }


}

void montee_triangulaire(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    int i,k;
    float t;
    float Tab[100];

    Tab[(P-1)]=(ptr2[Q * (P-1) + (Q-1)])/(ptr[M * (N-1) + (M-1)]);
    for(i=(P-2); i>-1; i--)
    {
        for(k=(i+1), t=0; k<(P); k++)
        {
           t += (ptr[M * i + k]) * Tab[k];
        }
        Tab[i] = ((ptr2[Q * i + 0]) - t) / (ptr[M * i + i]);
    }

    printf("\nLes solutions sont: \n");
    for(i=(P-1); i>-1; i--)
    {
        printf("x%d= %f\n", (i+1), Tab[i]);
    }

}

void elimination_gauss(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    int i, h, y, r;
    float pivot, t, w, k;
    for(i=0, pivot=0; i<N; i++)
    {
        pivot = ptr[M * i + i] ;
        //printf("Pivot= %d\n", pivot);
         //k=(ptr[M * i + (i-1)])/pivot;
        if(pivot!=0)
        {
         for(h=(i+1), w=0; h<N; h++)
         {
             k=(ptr[M * h + i])/pivot;
             for(y=0, t=0; y<M; y++)
             {
                 t = ptr[M * h + y] - k * ptr[M * i + y];
                 ptr[M * h + y] = t;

             }
             w = (ptr2[Q * h + 0]) - k * (ptr2[Q * i + 0]);
             ptr2[Q * h + 0] = w;

         }
        }
        else
        {
            //printf("Problème");
            float permut, permut1;
            /*int d;
            int invers_check = 0;
            for(int v=0; v<N; v++)
            {
                if(ptr[M * v + i] != 0)
                {
                    d = v;
                    //break;
                }
                else
                    invers_check++;

            }
            if(invers_check == N)
            {
                printf("\nDésolé la matrice est inversible.\n");
                i = N;
            }*/

            for(r=0; r<N; r++)
            {
                permut = ptr[M * i + r];
                ptr[M * i + r] = ptr[M * (i+1) + r];
                ptr[M * (i+1) + r] = permut;

            }
            permut1 = (ptr2[Q * i + 0]);
            ptr2[Q * i + 0] = ptr2[Q * (i+1) + 0];
            ptr2[Q * (i+1) + 0] = permut1;
            continue;
        }




    }

}

void decomposition_crout(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    float L[100][100];
    float U[100][100];
    int i, k ,j;
    float z, w, /*s,*/ e;
    for(i=0; i<N; i++)
    {
       L[i][0] = ptr[M * i + 0];


    }
    for(i=1; i<N; i++)
    {
        U[0][i] = (ptr[M * 0 + i]/L[0][0]);

    }
    for(i=1; i<(N); i++)
    {
        for(k=0, z=0; k<(i); k++)
        {
            z += L[i][k] * U[k][i];

        }
        L[i][i] = ptr[M * i + i] - z;



        for(j=(i+1); j<(N); j++)
        {
            for(k=0, w=0; k<(i); k++)
            {
                w += L[j][k] * U[k][i];

            }
            L[j][i] = ptr[M * j + i] - w;



            for(k=0, e=0; k<(i); k++)
            {
                e += L[i][k] * U[k][j];

            }
            U[i][j] = (ptr[M * i + j] - e) / (L[i][i]);


        }
    }
    /*for(k=0, s=0; k<(N); k++)
    {
        s = L[(N-1)][k] * U[k][(N-1)];
    }
    L[(N-1)][(N-1)] = ptr[M * (N-1) + (N-1)] - s;*/
    for(i=0; i<N; i++)
    {
        U[i][i] = 1;
    }

    /*printMatrice(&(L[0][0]), N, M);
    printf("\n");
    printMatrice(&(U[0][0]), N, M);*/

//    int i,j;
/*Affichage des matrices L et U */
    printf("\nAffichage de la matrice L\n");
    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            printf("%f\t", L[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("\nAffichage de la matrice U\n");
    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            printf("%f\t", U[i][j]);
        }
        printf("\n");
    }

    /*Descente triangulaire*/
    float t, Tab[100];
    Tab[0]=(ptr2[Q * 0 + 0])/(L[0][0]);
    for(i=1; i<P; i++)
    {
        for(k=0, t=0; k<i; k++)
        {
           t += (L[i][k]) * Tab[k];
        }
        Tab[i] = ((ptr2[Q * i + 0]) - t) / (L[i][i]);
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

    printf("\nLes solutions sont: \n");
    for(i=(P-1); i>-1; i--)
    {
        printf("x%d= %f\n", (i+1), Tab1[i]);
    }

}

void decomposition_choleski(float *ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    float L[100][100];
    float U[100][100];
    int i, k ,j;
    float w, z;

    L[0][0] = sqrt(ptr[M * 0 + 0]);

    for(i=1; i<N; i++)
    {
        L[i][0] = ptr[M * i + 0] / L[0][0];
    }

    for(k=1; k<N; k++)
    {
        for(j=0, w=0; j<k; j++)
        {
            w += pow(L[k][j], 2);
        }
        L[k][k] = sqrt(ptr[M * k + k] - w);

        for(i=(k+1); i<N; i++)
        {
            for(j=0, z=0; j<k; j++)
            {
                z += L[i][j] * L[k][j];
            }
            L[i][k] = (ptr[M * i + k] - z)/(L[k][k]);
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
    /*Affichage des matrices L et U */
    printf("\nAffichage de la matrice L\n");
    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            printf("%f\t", L[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("\nAffichage de la matrice U\n");
    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            printf("%f\t", U[i][j]);
        }
        printf("\n");
    }


    /*Descente triangulaire*/
    float t, Tab[100];
    Tab[0]=(ptr2[Q * 0 + 0])/(L[0][0]);
    for(i=1; i<P; i++)
    {
        for(k=0, t=0; k<i; k++)
        {
           t += (L[i][k]) * Tab[k];
        }
        Tab[i] = ((ptr2[Q * i + 0]) - t) / (L[i][i]);
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

    printf("\nLes solutions sont: \n");
    for(i=(P-1); i>-1; i--)
    {
        printf("x%d= %f\n", (i+1), Tab1[i]);
    }
}

void elimination_gauss_jordan(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    int r, i, j, k, y;
    float max, pivot, echange;
    for(j=0, r=-1; j<M; j++)
    {
        max = 0;
        for(i=0; i<N; i++)
        {
           if(max < fabs(ptr[M * (i) + j]))
           {
               max =  fabs(ptr[M * i + j]);
               k = i;
           }
        }
        pivot = ptr[M * k + j];
        if(pivot != 0)
        {
            r++;



            if(k != r)
            {
                /*Matrice A*/
               for(y=0; y<M; y++)
               {
                   echange = ptr[M * k + y];
                   ptr[M * k + y] = ptr[M * r + y];
                   ptr[M * r + y] = echange;
                }

               /*Matrice B*/
               echange = ptr2[Q * k + 0];
               ptr2[Q * k + 0] = ptr2[Q * r + 0];
               ptr2[Q * r + 0] = echange;
            }

            /*Normalisation*/
                for(y=0; y<M; y++)
                {
                    ptr[M * k + y] = ptr[M * k + y] / pivot;

                }
                ptr2[Q * k + 0] = ptr2[Q * k + 0] / pivot;


            for(i=0; i<N; i++)
            {
                if(i != r)
                {
                    for(y=0; y<M; y++)
                    {
                        ptr[M * i + y] = ptr[M * i + y] - (ptr[M * i + j] * ptr[M * r + y]);
                    }
                    ptr2[Q * i + 0] = ptr2[Q * i + 0] - (ptr[M * i + j] * ptr2[Q * r + 0]);
                }
            }

        }


    }

}

void inverse_matrice_carree(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    int r, i, j, k, y;
    float max, pivot, echange;
    for(j=0, r=-1; j<M; j++)
    {
        max = 0;
        for(i=0; i<N; i++)
        {
           if(max < fabs(ptr[M * (i) + j]))
           {
               max =  fabs(ptr[M * i + j]);
               k = i;
           }
        }
        pivot = ptr[M * k + j];
        if(pivot != 0)
        {
            r++;
            if(k != r)
            {

               for(y=0; y<M; y++)
               {
                   /*Matrice A*/
                   echange = ptr[M * k + y];
                   ptr[M * k + y] = ptr[M * r + y];
                   ptr[M * r + y] = echange;

                   /*Matrice B*/
                   echange = ptr2[Q * k + y];
                   ptr2[Q * k + y] = ptr2[Q * r + y];
                   ptr2[Q * r + y] = echange;
                }

            }
            for(y=0; y<M; y++)
            {
                ptr[M * k + y] = ptr[M * k + y] / pivot;
                ptr2[Q * k + y] = ptr2[Q * k + y] / pivot;
                //ptr[M * k + j] /= pivot;
            }


            for(i=0; i<N; i++)
            {
                if(i != r)
                {
                    for(y=0; y<M; y++)
                    {
                        ptr[M * i + y] = ptr[M * i + y] - (ptr[M * i + j] * ptr[M * r + y]);
                        ptr2[Q * i + y] = ptr2[Q * i + y] - (ptr[M * i + j] * ptr2[Q * r + y]);
                    }

                }
            }

        }


    }

}

void jacobi(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q, int *iteration)
{
    float x[100];
    for(int i=0; i<100; i++)
    {
        x[i] = 0;
    }
    int i, j, k, t, indi;
    float w, echange, max;

    for(i=0; i<N; i++)
    {
        if(ptr[M * i + i] == 0)
        {
            for(j=0, max=0; j<N; j++)
            {
                if(max < fabs(ptr[M * j + i]))
                {
                    max = fabs(ptr[M * j + i]);
                    indi = j;
                }
               /*if(j!=i && ptr[M * j + i] != 0)
               {
                   indi = j;
               }*/
            }
            for(j=0; j<N; j++)
            {
                echange = ptr[M * i + j];
                ptr[M * i + j] = ptr[M * indi + j];
                ptr[M * indi + j] = echange;
            }
            echange = ptr2[M * 0 + i];
            ptr2[M * 0 + i] = ptr2[M * 0 + indi];
            ptr2[M * 0 + indi] = echange;
        }
    }

    /*Vérification de la diagonale dominante*/
    float Tab[100][100];
    int s=0;
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            //Tab[i][j] = A[0][s++];
            Tab[i][j] = ptr[M * 0 + s];
            s++;

        }
        //printf("\n");
    }
    /*for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            printf("%f", Tab[i][j]);
        }
        printf("\n");
    }*/

    if(isDiagonalDominant(Tab, N))
    {
        for(t=0; t<*iteration; t++)
        {

            for(i=0; i<N; i++)
            {
                for(j=0, w=0; j<N; j++)
                {
                    if(j != i)
                    {
                        w += ptr[M * i + j] * x[j];
                    }
                }

                x[(i+N)] = (1 / ptr[M * i + i]) * (ptr2[Q * i + 0] - w);
            }

            for(i=0; i<N; i++)
            {
                x[i] = x[(i+N)];
            }

            printf("\n A la %de itération les solutions sont: \n", (t+1));
            for(k=0; k<N; k++)
            {
                printf("x%d= %f\n", (k+1), x[k]);
            }
        }
        printf("\n\n\n");
        printf("\nLes solutions sont: \n");
        for(k=0; k<N; k++)
        {
            printf("x%d= %f\n", (k+1), x[k]);
        }
    }
    else
    {
        printf("\nMême après de potentiels échanges de lignes, cette matrice n'est pas à diagonale dominante. Cette matrice ne converge pas pour cette méthode.\n");

    }
}

void gauss_seidel(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q, int *iteration)
{
    float x[100];
    for(int i=0; i<100; i++)
    {
        x[i] = 0;
    }
    int i, j, k, t, indi;
    float w, echange, max;



    for(i=0; i<N; i++)
    {
        if(ptr[M * i + i] == 0)
        {
            for(j=0, max=0; j<N; j++)
            {
                if(max < fabs(ptr[M * j + i]))
                {
                    max = fabs(ptr[M * j + i]);
                    indi = j;
                }

            }
            for(j=0; j<N; j++)
            {
                echange = ptr[M * i + j];
                ptr[M * i + j] = ptr[M * indi + j];
                ptr[M * indi + j] = echange;
            }
            echange = ptr2[M * 0 + i];
            ptr2[M * 0 + i] = ptr2[M * 0 + indi];
            ptr2[M * 0 + indi] = echange;
        }
    }

    /*for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            printf("%f\t", (ptr[M * i + j]));
        }
        printf("\n");
    }
    printf("\n");
    for(i=0; i<P; i++)
    {
        for(j=0; j<Q; j++)
        {
            printf("%f\t", (ptr2[M * i + j]));
        }
        printf("\n");
    }*/

    /*Vérification de la diagonale dominante*/
    float Tab[100][100];
    int s=0;
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            //Tab[i][j] = A[0][s++];
            Tab[i][j] = ptr[M * 0 + s];
            s++;

        }
        //printf("\n");
    }
    /*for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        {
            printf("%f", Tab[i][j]);
        }
        printf("\n");
    }*/

    if(isDiagonalDominant(Tab, N))
    {
        for(t=0; t<*iteration; t++)
        {
            for(i=0; i<N; i++)
            {
                for(j=0, w=0; j<N; j++)
                {
                    if(j != i)
                    {
                        w += ptr[M * i + j] * x[j];
                    }
                }

                x[i] = (1 / ptr[M * i + i]) * (ptr2[Q * i + 0] - w);
            }

            printf("\n A la %de itération les solutions sont: \n", (t+1));
            for(k=0; k<N; k++)
            {
                printf("x%d= %f\n", (k+1), x[k]);
            }
        }
        printf("\n\n\n");
        printf("\nLes solutions sont: \n");
        for(k=0; k<N; k++)
        {
            printf("x%d= %f\n", (k+1), x[k]);
        }
    }
    else
    {
        printf("Même après de potentiels échanges de lignes, cette matrice n'est pas à diagonale dominante. Cette matrice ne converge peut-être pas pour cette méthode.\n");
        char reponse;
        do
            {
                printf("Voulez-vous quand même continuer (O/N) ?\n");
                reponse = getche();
                reponse = toupper(reponse);
                fflush(stdin);
            }
            while((reponse!='O')&&(reponse!='N'));

            if(reponse=='O')
            {
                for(t=0; t<*iteration; t++)
                {
                    for(i=0; i<N; i++)
                    {
                        for(j=0, w=0; j<N; j++)
                        {
                            if(j != i)
                            {
                                w += ptr[M * i + j] * x[j];
                            }
                        }

                        x[i] = (1 / ptr[M * i + i]) * (ptr2[Q * i + 0] - w);
                    }

                    printf("\n A la %de itération les solutions sont: \n", (t+1));
                    for(k=0; k<N; k++)
                    {
                        printf("x%d= %f\n", (k+1), x[k]);
                    }
                }
                printf("\n\n\n");
                printf("\nLes solutions sont: \n");
                for(k=0; k<N; k++)
                {
                    printf("x%d= %f\n", (k+1), x[k]);
                }
            }
    }
}

void elimination_gauss_jordan_advanced(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    int r, i, j, k, y;
    float max, pivot, echange, passeur;

    for(j=0, r=0; j<N; j++)
    {
        /*Recherche du maximum dans la première pour le placer sur le r équivalent*/
        max = 0;
        for(i=r; i<M; i++)
        {
            if(max < fabs(ptr[M * i + j]))
            {
                max = fabs(ptr[M * i + j]);
                k = i;
            }
        }

        /*Si l'indice du maximum ne correspond pas, on effectue sur le champ l'echange de ligne*/
        if(k != r)
        {
            for(i=0; i<N; i++)
            {
                echange = ptr[M * k + i];
                ptr[M * k + i] = ptr[M * r + i];
                ptr[M * r + i] = echange;

            }
            for(i=0; i<Q; i++)
            {

                echange = ptr2[Q * k + i];
                ptr2[Q * k + i] = ptr2[Q * r + i];
                ptr2[Q * r + i] = echange;
            }

        }


        /*Division sur la ligne pour obtenir 1 sur la diagonale*/
        pivot = ptr[M * r + j];
        for(i=r; i<M; i++)
        {
            ptr[M * r + i] = ptr[M * r + i] / pivot;

        }
        for(i=0; i<Q; i++)
        {

            ptr2[Q * r + i] = ptr2[Q * r + i] / pivot;
        }

        /*Normalisation de la colonne correspondante*/
        for(i=0; i<N; i++)
        {
            passeur = ptr[M * i + r];
            if(i!=r)
            {
                for(y=0; y<N; y++)
                {

                        ptr[M * i + y] = ptr[M * i + y] - (passeur * ptr[M * j + y]);



                }
                for(y=0; y<Q; y++)
                {

                        ptr2[Q * i + y] = ptr2[Q * i + y] - (passeur * ptr2[Q * j + y]);

                }

            }
        }

        r++;

    }
}

float elimination_gauss_jordan_advanced_determinant(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q)
{
    int r, i, j, k, y, p=0;
    float max, pivot, echange, passeur, determinant, Tab[100];

    for(j=0, r=0; j<N; j++)
    {
        /*Recherche du maximum dans la première pour le placer sur le r équivalent*/
        max = 0;
        for(i=r; i<M; i++)
        {
            if(max < fabs(ptr[M * i + j]))
            {
                max = fabs(ptr[M * i + j]);
                k = i;
            }
        }

        /*Si l'indice du maximum ne correspond pas, on effectue sur le champ l'echange de ligne*/
        if(k != r)
        {
            p++;
            for(i=0; i<N; i++)
            {
                echange = ptr[M * k + i];
                ptr[M * k + i] = ptr[M * r + i];
                ptr[M * r + i] = echange;

            }
            for(i=0; i<Q; i++)
            {

                echange = ptr2[Q * k + i];
                ptr2[Q * k + i] = ptr2[Q * r + i];
                ptr2[Q * r + i] = echange;
            }

        }


        /*Division sur la ligne pour obtenir 1 sur la diagonale*/
        Tab[j] = ptr[M * r + j];
        pivot = ptr[M * r + j];
        for(i=r; i<M; i++)
        {
            ptr[M * r + i] = ptr[M * r + i] / pivot;

        }
        for(i=0; i<Q; i++)
        {

            ptr2[Q * r + i] = ptr2[Q * r + i] / pivot;
        }

        /*Normalisation de la colonne correspondante*/
        for(i=0; i<N; i++)
        {
            passeur = ptr[M * i + r];
            if(i!=r)
            {
                for(y=0; y<N; y++)
                {

                        ptr[M * i + y] = ptr[M * i + y] - (passeur * ptr[M * j + y]);



                }
                for(y=0; y<Q; y++)
                {

                        ptr2[Q * i + y] = ptr2[Q * i + y] - (passeur * ptr2[Q * j + y]);

                }

            }
        }

        r++;

    }
    /*for(i=0; i<N; i++)
    {
        printf("Tab[%d] = %lf\n", (i), Tab[i]);
        if(Tab[i] < 0 || Tab[i] > 0)
        {
            printf("Hello I'm not here\n");
        }
    }*/
    for(i=0, determinant=1; i<N; i++)
    {
        //printf("\nAvant %f *= %f", determinant, Tab[i]);
        determinant *= Tab[i];
        //printf("\nAprès %f *= %f", determinant, Tab[i]);

    }
    determinant = pow(-1,p) * determinant;
    //printf("\nLe déterminant de cette matrice est %f", determinant);
    return determinant;
}

void thomas_tridiagonal_method(float A[][100], float B[], int N, int M)
{
    float m;
    float X[100];
    int i;
    //forward elimination

    for (i=1;i<N;i++)
     {
      m = A[i][i-1]/A[i-1][i-1];
      A[i][i] = A[i][i] - (m * A[i-1][i]);
      B[i] = B[i] - (m * B[i-1]);
     }

  //backward substitution

     X[N-1]=B[N-1]/A[N-1][N-1];

     for(i=N-2;i>=0;i--)
    {
      X[i]=(B[i]-A[i][i+1]*X[i+1])/A[i][i];
    }
      printf("\nLes solutions sont:\n");
      for (i=0;i<N;i++)
      printf("x%d = %4.2f\n", (i+1), X[i]);
}

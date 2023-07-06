#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <conio.h>
#include <stdbool.h>
#include <ctype.h>
#include <locale.h>

void delay(int sec)
    {
     clock_t prec, suiv;
     int pause;
     pause = sec * CLOCKS_PER_SEC;
     prec = suiv = clock();
     while(suiv - prec < pause)
     {
         suiv=clock();
     }
    }
//Prototype de la procédure permettant de vérifier la loyauté des nombres entiers entrés
void number_check_int(char *texte, int *z)
    {
        _Bool VERIFY;
      do
                {

                printf("%s", texte);
                VERIFY=scanf("%i", z);
                if(!VERIFY)
                {
                    printf("Ce n'est pas un nombre réel.\n Veuillez réessayer\n");
                    fflush(stdin);
                    delay(2);
                }
                }while(!VERIFY);

    }

//Prototype de la procédure permettant de vérifier la loyauté des nombres décimaux entrés
void number_check_double(char *texte, double *z)
    {
        _Bool VERIFY;
      do
                {

                printf("%s", texte);
                VERIFY=scanf("%lf", z);
                if(!VERIFY)
                {
                    printf("Ce n'est pas un nombre réel.\n Veuillez réessayer\n");
                    fflush(stdin);
                    delay(2);
                }
                }while(!VERIFY);

    }

void number_check_float(char *texte, float *z)
{
    _Bool VERIFY;
      do
                {

                printf("%s", texte);
                VERIFY=scanf("%f", z);
                if(!VERIFY)
                {
                    printf("Ce n'est pas un nombre réel.\n Veuillez réessayer\n");
                    fflush(stdin);
                    delay(2);
                }
                }while(!VERIFY);
}
//Prototype de la fonction permttant de calculer un polynôme de degré 3
double result_fonction(double x, double *a, double *b, double *c, double *d)
{
    return ((*a*(pow(x,3))+(*b*(pow(x,2))+(*c*x)+*d)));
}



/*Prototype de la fonction qui scanne n'importe quelle polynôme*/
void polynome(double Tab[], int N)
    {
        int I;
        for(I=0; I<(N+1); I++)
        {
            printf("Entrer le %de coefficient", (I+1));
            number_check_double(": ", &Tab[I]);
            printf("\n");
        }
        printf("Le polynôme entré est :");
        int k=N;
        for(I=0; I<(N+1); I++)
        {
            printf("%fx^%d", Tab[I], k);
            k--;
            if(I!=N)
                printf(" + ");
        }

    }

double calc_polynome(double Tab[], double *x, int N)
{
        double result=0;
        /*int M=N;
        for(int I=0; I<(N+1); I++)
        {
            result += Tab[I]*pow(*x,M);
            M--;

        }*/
        //return result = pow(*x, 3) + pow(*x, 2) - 3**x -3;
        return result = *x - 2 * tan(*x) + 1;
        //return result = log(*x);
}

void table_initialisation(double Tab[])
{
    for(int i=0; i<100; i++)
    {
        Tab[i] = 0;
    }
}

void matric_recop(float *fin, int N, int M, float *debut)
{
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<M; j++)
        {
            fin[M * i + j] = debut[M * i + j];
        }
    }

}

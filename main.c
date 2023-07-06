#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <locale.h>
#include <conio.h>
#include <ctype.h>
#include <time.h>
#include <stdbool.h>
#include "control_fonction.h"
#include "equation_lineaire.h"
#include "interpolation_lineaire.h"
#include "equation_differentielle.h"
#include "save_me.h"





int main()
{
    system("color F0");
    setlocale(LC_CTYPE,"");
    int ans;
    char reponse='O';
    /*double vin = 2;
    float kev = polynome(&vin);
    printf("Kev = %f", kev);
    fflush(stdin);*/

    /*Boucle de base du programme*/
    while (reponse=='O')
    {
        fflush(stdin);
        printf("------------Bienvenue en Analyses Numériques-----------\n\n");
        printf("1-Résolution des équations non linéaires\n\n");
        printf("2-Résolution des systèmes d'équations linéaires\n\n");
        printf("3-Interpolation linéaire\n\n");
        printf("4-Equations différentielles\n\n");
        printf("Veuillez faire un choix: ");
        scanf("%d", &ans);

        switch(ans)
        {
            case 1:
                {
                    double Tab[100];
                    double e;
                    int N = 0;
                    int n_iteration;
                    printf("\n-------Bonjour et bienvenue dans la résolution des équations non linéaires-------\n");
                    /*number_check_int("Entrez le degré du polynôme :", &N);
                    printf("Veuillez entrer les coefficients du pôlynome dans l'ordre décroissant\n");
                    polynome(Tab, N);*/

                    number_check_double("\nVeuillez entrer la condition d'arrêt: ", &e);

                    number_check_int("\nVeuillez entrer le nombre maximal d'itérations: ", &n_iteration);

                    char answer='O';

                    /*Boucle de la partie équations non linéaires*/
                    while(answer=='O')
                    {
                    fflush(stdin);
                    int method;


                    printf("\nNous disposons de diverses méthodes.\n\n");
                    printf("1- Méthode de la bissection\n\n");
                    printf("2- Méthode de Lagrange\n\n");
                    printf("3- Méthode des points des fixes\n\n");
                    printf("4- Méthode de la sécante\n\n");
                    printf("5- Méthode de la corde 1\n\n");
                    printf("6- Méthode de la corde 2\n\n");
                    printf("7- Méthode de Newton\n\n");

                    printf("Veuillez choisir une méthode: ");
                    scanf("%i", &method);

                    switch(method)
                    {
                        case 1:
                            {

                                /*Méthode de la bissection*/
                                double x1, x2, xm, fx1, fx2, fxm;
                                int i;
                                printf("\n----Bienvenue dans l'algorithme de la bissection----\n");

                                printf("Veuillez entrer l'intervalle [x1 ; x2]\n");
                                number_check_double("Entrer x1: ", &x1);
                                number_check_double("Entrer x2: ", &x2);


                                for (i=0; i<n_iteration; i++)
                                {
                                    printf("-----%de itération-----\n", (i+1));
                                    xm=(x1+x2)/2;
                                    fx1 = calc_polynome(Tab, &x1, N);
                                    fx2 = calc_polynome(Tab, &x2, N);
                                    fxm = calc_polynome(Tab, &xm, N);

                                    if (((fabs((x2-x1)))/(2*fabs(xm)))<e)
                                    {
                                        printf("Les valeurs sont:\nx1= %f\nx2= %f\nxm= %f\nf(x1)= %f\nf(x2)= %f\nf(xm)= %f\n", x1, x2, xm, fx1, fx2, fxm);

                                        printf("La convergence est atteinte.\nLa racine xm est égal à %f. f(xm) est égal à %f\n", xm, fxm);
                                        printf("Arret.\n");
                                        i=n_iteration;
                                        //break;


                                    }
                                    else
                                    {
                                        printf("Les valeurs sont:\nx1= %f\nx2= %f\nxm= %f\nf(x1)= %f\nf(x2)= %f\nf(xm)= %f\n", x1, x2, xm, fx1, fx2, fxm);

                                        if ((fx1*fxm)<0)
                                        {
                                            x2=xm;
                                        }

                                        if ((fxm*fx2)<0)
                                        {
                                            x1=xm;
                                        }
                                        printf("\n");

                                        if ((fx1*fx2)>0)
                                        {
                                            printf("Il y'a un nombre paire de solutions dans l'intervalle y compris 0.\n");
                                            i=n_iteration;
                                            //break;
                                        }

                                        if (i==(n_iteration-1))
                                        {
                                            printf("Convergence non atteinte en %d itérations\n", n_iteration);
                                        }
                                    }

                                }

                                break;
                            }
                        case 2:
                            {
                                //table_initialisation(Tab);
                                /*Méthode de Lagrange*/
                                double x1, x2, xm, fx1, fx2, fxm;
                                int i;
                                printf("\n----Bienvenue dans l'algorithme de Lagrange----\n");


                                printf("Veuillez entrer l'intervalle [x1 ; x2]\n");
                                number_check_double("Entrer x1: ", &x1);
                                number_check_double("Entrer x2: ", &x2);


                                for (i=0; i<n_iteration; i++)
                                {
                                    printf("-----%de itération-----\n", (i+1));
                                    fx1 = calc_polynome(Tab, &x1, N);
                                    fx2 = calc_polynome(Tab, &x2, N);
                                    /*fx1=result_fonction(x1, &a, &b, &c, &d);
                                    fx2=result_fonction(x2, &a, &b, &c, &d);*/
                                    xm = x1 - ((x2-x1)/(fx2-fx1))*fx1;
                                    //xm=((x1*fx2-x2*fx1)/(fx2-fx1));
                                    fxm = calc_polynome(Tab, &xm, N);
                                    //fxm=result_fonction(xm, &a, &b, &c, &d);


                                    if (((fabs((x2-x1)))/(2*fabs(xm)))<e)
                                    {
                                        printf("Les valeurs sont:\nx1= %f\nx2= %f\nxm= %f\nf(x1)= %f\nf(x2)= %f\nf(xm)= %f\n", x1, x2, xm, fx1, fx2, fxm);

                                        printf("La convergence est atteinte.\nLa racine xm est égal à %f. f(xm) est égal à %f\n", xm, fxm);
                                        printf("Arret.\n");
                                        break;


                                    }
                                    else
                                    {
                                        printf("Les valeurs sont:\nx1= %f\nx2= %f\nxm= %f\nf(x1)= %f\nf(x2)= %f\nf(xm)= %f\n", x1, x2, xm, fx1, fx2, fxm);

                                        if ((fx1*fxm)<0)
                                        {
                                            x2=xm;
                                        }

                                        if ((fxm*fx2)<0)
                                        {
                                            x1=xm;
                                        }
                                        printf("\n");

                                        if ((fx1*fx2)>0)
                                        {
                                            printf("Il y'a un nombre paire de solutions dans l'intervalle y compris 0.\n");
                                            i=n_iteration;
                                            //break;
                                        }

                                        if (i==(n_iteration-1))
                                        {
                                            printf("Convergence non atteinte en %d itérations\n", n_iteration);
                                        }
                                    }

                                }
                                break;
                            }
                        case 3:
                            {
                                /*Méthode des points fixes*/
                                double x0, xk, gxk;
                                int i;
                                printf("\n----Bienvenue dans l'algorithme des points fixes----\n");

                                number_check_double("Veuillez entrer un x0, valeur estimée initiale du point fixe: ", &x0);

                                xk=x0;

                                for (i=0; i<n_iteration; i++)
                                {
                                    printf("-----%de itération-----\n", (i+1));

                                    //gxk = cos(xk);
                                    gxk = atan((xk+1)/2);
                                    //gxk = result_fonction(xk, &a, &b, &c, &d);
                                    //gxk = (1 / (1+xk));

                                    if ((fabs((gxk-xk))/fabs(gxk))<e)
                                    {

                                        printf("La convergence est atteinte.\nLa solution est: %lf\n", gxk);
                                        printf("Arret.\n");
                                        break;
                                    }
                                    else
                                    {
                                        printf("Les valeurs sont:\nx%d= %lf\ng(x%d)= %lf\n", i, xk, i, gxk);

                                        xk=gxk;

                                        if (i==(n_iteration-1))
                                        {
                                            printf("Convergence non atteinte en %d itérations\n", n_iteration);
                                        }
                                    }

                                }
                                break;
                            }
                        case 4:
                            {
                                /*Méthode de la sécante*/
                                double x0, x1, x2;
                                int i;
                                printf("\n----Bienvenue dans l'algorithme de la sécante----\n");

                                number_check_double("Veuillez entrer un x0, première valeur initiale: ", &x0);

                                number_check_double("Veuillez entrer un x1, deuxième valeur initiale: ", &x1);

                                for (i=0; i<n_iteration; i++)
                                {
                                    printf("-----%de itération-----\n", (i+1));

                                    //x2 = x1 - ((cos(x1)*(x1-x0))/(cos(x1)-cos(x0)));
                                    //x2= x1 - (result_fonction(x1, &a, &b, &c, &d)*(x1-x0))/((result_fonction(x1, &a, &b, &c, &d)-(result_fonction(x0, &a, &b, &c, &d))));
                                    x2 = x1 - (calc_polynome(Tab, &x1, N) * (x1-x0)) / ((calc_polynome(Tab, &x1, N) - (calc_polynome(Tab, &x0, N))));

                                    if ((fabs((x2-x1))/fabs(x2))<e)
                                    {

                                        printf("La convergence est atteinte.\nLa solution est: %lf\n", x2);
                                        printf("Arret.\n");
                                        break;
                                    }
                                    else
                                    {
                                        printf("Les valeurs sont:\nx%d= %lf\nx%d= %lf\nx%d= %lf\n", i, x0, (i+1), x1, (i+2), x2);

                                        x0=x1;
                                        x1=x2;

                                        if (i==(n_iteration-1))
                                        {
                                            printf("Convergence non atteinte en %d itérations\n", n_iteration);
                                        }
                                    }

                                }
                                break;
                            }
                        case 5:
                            {
                                /*Méthode de la corde 1*/
                                double x1, x2, xk, gxk;
                                int i;
                                printf("\n----Bienvenue dans l'algorithme de la corde 1----\n");
                                //printf("Soit f(x) = cos(x)\n");
                                printf("Veuillez entrer l'intervalle [x1 ; x2]\n");
                                number_check_double("Entrez x1: ", &x1);
                                number_check_double("Entrez x2: ", &x2);
                                number_check_double("Veuillez entrer un x0 qui appartient à l'intervalle entré ci-dessus pour commencer: ", &xk);

                                for (i=0; i<n_iteration; i++)
                                {
                                    printf("-----%de itération-----\n", (i+1));

                                    //gxk = xk - ((x2-x1)/(cos(x2)-cos(x1)))*cos(xk);

                                    //gxk = xk - ((x2-x1)/(result_fonction(x2, &a, &b, &c, &d)-result_fonction(x1, &a, &b, &c, &d)))*result_fonction(xk, &a, &b, &c, &d);

                                    gxk = xk - ((x2-x1) / (calc_polynome(Tab, &x2, N) - calc_polynome(Tab, &x1, N))) * calc_polynome(Tab, &xk, N);

                                    if ((fabs((gxk-xk))/fabs(gxk))<e)
                                    {
                                        printf("La convergence est atteinte.\nLa solution est: %lf\n", gxk);
                                        printf("Arret.\n");
                                        break;
                                    }
                                    else
                                    {
                                        printf("Les valeurs sont:\nx%d= %lf\ng(x%d)= %lf\n", i, xk, i, gxk);

                                        xk=gxk;

                                        if (i==(n_iteration-1))
                                        {
                                            printf("Convergence non atteinte en %d itérations\n", n_iteration);
                                        }
                                    }

                                }

                                break;
                            }
                        case 6:
                            {
                                //Méthode de la corde 2
                                double x0, xk, gxk;
                                int i;
                                printf("\n----Bienvenue dans l'algorithme de la corde 2----\n");
                                //printf("Soit g(x) = cos(x)\n");
                                number_check_double("Veuillez entrer un x0, valeur estimée initiale du point fixe: ", &x0);

                                xk=x0;

                                for (i=0; i<n_iteration; i++)
                                {
                                    printf("-----%de itération-----\n", (i+1));

                                    //gxk = xk - ((cos(xk))/(-sin(x0)));

                                    //gxk = xk - ((result_fonction(xk, &a, &b, &c, &d))/(2*x0+1));
                                    gxk = xk - ((calc_polynome(Tab, &xk, N)) / (1-2*(1+pow(tan(xk), 2))));

                                    if ((fabs((gxk-xk))/fabs(gxk))<e)
                                    {

                                        printf("La convergence est atteinte.\nLa solution est: %lf\n", gxk);
                                        printf("Arret.\n");
                                        break;
                                    }
                                    else
                                    {
                                        printf("Les valeurs sont:\nx%d= %lf\ng(x%d)= %lf\n", i, xk, i, gxk);

                                        xk=gxk;

                                        if (i==(n_iteration-1))
                                        {
                                            printf("Convergence non atteinte en %d itérations\n", n_iteration);
                                        }
                                    }

                                }
                                break;
                            }
                        case 7:
                            {
                                //Méthode de Newton
                                double x0, xk, gxk;
                                int i;
                                printf("\n----Bienvenue dans l'algorithme de Newton----\n");
                                //printf("Soit g(x) = cos(x)\n");

                                number_check_double("Veuillez entrer un x0, valeur estimée initiale du point fixe: ", &x0);

                                xk=x0;

                                for (i=0; i<n_iteration; i++)
                                {
                                    printf("-----%de itération-----\n", (i+1));

                                    //gxk = xk - ((cos(xk))/(-sin(xk)));

                                    //gxk = xk - ((result_fonction(xk, &a, &b, &c, &d))/(2*xk-2));
                                    gxk = xk - ((calc_polynome(Tab, &xk, N)) / (1-2*(1+pow(tan(xk), 2))));


                                    if ((fabs((gxk-xk))/fabs(gxk))<e)
                                    {

                                        printf("La convergence est atteinte.\nLa solution est: %lf\n", gxk);
                                        printf("Arret.\n");
                                        break;
                                    }
                                    else
                                    {
                                        printf("Les valeurs sont:\nx%d= %lf\ng(x%d)= %lf\n", i, xk, i, gxk);

                                        xk=gxk;

                                        if (i==(n_iteration-1))
                                        {
                                            printf("Convergence non atteinte en %d itérations\n", n_iteration);
                                        }
                                    }

                                }
                                break;
                            }
                        default:
                            {
                                printf("Méthode incorrecte. Veuillez choisir une méthode correcte.\n\n");
                                delay(3);
                                system("cls");
                                continue;

                            }

                    //Boucle de reprise de la partie équations non linéaires
                    }
                    do
                        {
                            printf("\n\nVoulez-vous recommencer la partie équations non linéaires (O/N) ? ");
                            answer = getche();
                            answer = toupper(answer);
                            fflush(stdin);
                        }
                        while((answer!='O')&&(answer!='N'));

                        if(answer=='N')
                        {
                            answer='N';
                        }
                        if(answer=='O')
                        {
                            //system("cls");
                        }
                    }
                    break;
                }

            case 2:
                {
                    float A[100][100], B[100][100];
                    float Arecop[100][100], Brecop[100][100];
                    int N, M, P, Q;
                    printf("\n--------------Bienvenue dans la partie des systèmes d'équations linéeaires-------------\n\n");
                    number_check_int("Entrer le nombre de lignes de la matrice A: ", &N);
                    number_check_int("Entrer le nombre de colonnes de la matrice A: ", &M);

                    number_check_int("Entrer le nombre de lignes de la matrice B: ", &P);
                    number_check_int("Entrer le nombre de colonnes de la matrice B: ", &Q);



                    scanMatrice(&(A[0][0]), N, M);
                    printf("\n");
                    printMatrice(&(A[0][0]), N, M);
                    printf("\n");
                    scanMatrice(&(B[0][0]), P, Q);
                    printMatrice(&(B[0][0]), P, Q);
                    printf("\n");
                    matric_recop(&Arecop[0][0], N, M, &A[0][0]);
                    matric_recop(&Brecop[0][0], P, Q, &B[0][0]);
                    printMatrice(&Arecop[0][0], N, M);
                    printf("\n");
                    printMatrice(&Brecop[0][0], P, Q);
                    char answer='O';

                    /*Boucle de la partie systèmes d'équations algébriques*/
                    while(answer=='O')
                    {
                    fflush(stdin);

                   int method;
                   printf("\n1-Elimination de GAUSS\n\n");
                   printf("2-Décomposition LU de Crout\n\n");
                   printf("3-Décomposition LU factorisation de Choleski\n\n");
                   printf("4-Elimination de GAUSS JORDAN\n\n");
                   printf("5-Méthode de Jacobi\n\n");
                   printf("6-Méthode de Gauss-Seidel\n\n");
                   printf("7-Méthode de Thomas pour matrices tridiagonles\n\n");
                   printf("Veuillez faire un choix: ");
                   scanf("%d", &method);

                   switch(method)
                   {
                        case 1:
                            {
                                printf("\n\n----Bienvenue dans l'élimination de Gauss----\n");
                                matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                matric_recop(&B[0][0], P, Q, &Brecop[0][0]);

                                float Tab[100][100];
                                int s=0;
                                for(int i=0; i<N; i++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        Tab[i][j] = A[0][s++];
                                    }
                                    //printf("\n");
                                }



                                if(isInvertible(Tab, N))
                                {
                                    matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                    matric_recop(&B[0][0], P, Q, &Brecop[0][0]);
                                    elimination_gauss(&(A[0][0]), &(B[0][0]), N, M, P, Q);
                                    montee_triangulaire(&(A[0][0]), &(B[0][0]), N, M, P, Q);
                                }else
                                {
                                    printf("\nLa matrice n'est pas inversible. Il y a 0 ou une infinité de solutions. La méthode de GAUSS ne peut pas être utilisée.\n");
                                }

                                break;
                            }
                        case 2:
                            {
                                printf("\n\n----Bienvenue dans la décomposition LU de Crout----\n");
                                matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                matric_recop(&B[0][0], P, Q, &Brecop[0][0]);

                                float Tab[100][100];
                                int s=0;
                                for(int i=0; i<N; i++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        Tab[i][j] = A[0][s++];
                                    }
                                    //printf("\n");
                                }



                                if(isInvertible(Tab, N))
                                {
                                    matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                    matric_recop(&B[0][0], P, Q, &Brecop[0][0]);
                                    decomposition_crout(&(A[0][0]), &(B[0][0]), N, M, P, Q);
                                }else
                                {
                                    printf("\nLa matrice n'est pas inversible. Il y a 0 ou une infinité de solutions. La méthode de Crout ne peut pas être utilisée.\n");
                                }


                                break;
                            }
                        case 3:
                            {
                                printf("\n\n----Bienvenue dans la décomposition LU de Choleski----\n");
                                matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                matric_recop(&B[0][0], P, Q, &Brecop[0][0]);

                                float Tab[100][100];
                                int s=0;
                                for(int i=0; i<N; i++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        Tab[i][j] = A[0][s++];
                                    }
                                    //printf("\n");
                                }



                                if(isInvertible(Tab, N))
                                {
                                    /*for(int i=0; i<N; i++)
                                    {
                                        for(int j=0; j<N; j++)
                                        {
                                            printf("%f", Tab[i][j]);
                                        }
                                        printf("\n");
                                    }*/
                                    if(isSymetric(Tab, N))
                                    {
                                        if(defined_positive(Tab, N))
                                        {
                                            matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                            matric_recop(&B[0][0], P, Q, &Brecop[0][0]);
                                            decomposition_choleski(&A[0][0], &B[0][0], N, M, P, Q);
                                        }else
                                        {
                                            printf("La matrice n'est pas définie positive. La méthode de Choleski ne peut pas être utilisée\n");
                                        }

                                    }
                                    else
                                    {
                                        printf("La matrice n'est pas symétrique. La méthode de Choleski ne peut pas être utilisée\n");
                                    }

                                }else
                                {
                                    printf("\nLa matrice n'est pas inversible. Il y a 0 ou une infinité de solutions. La méthode de Choleski ne peut pas être utilisée.\n");
                                }



                                //break;
                            }
                        case 4:
                            {
                                printf("\n\n----Bienvenue dans l'élimination de GAUSS JORDAN----\n");
                                matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                matric_recop(&B[0][0], P, Q, &Brecop[0][0]);

                                elimination_gauss_jordan_advanced_determinant(&(A[0][0]), &(B[0][0]), N, M, P, Q);
                                printf("\n");
                                printf("Affichage de la matrice identité résultante\n");
                                printMatrice(&(A[0][0]), N, M);
                                printf("\n");
                                printf("Affichage de la matrice B conteanat directement les solutions\n");
                                printMatrice(&(B[0][0]), P, Q);
                                break;
                            }
                        case 5:
                            {
                                printf("\n\n----Bienvenue dans la méthode de Jacobi----\n");
                                matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                matric_recop(&B[0][0], P, Q, &Brecop[0][0]);

                                int iteration;
                                printf("Veuillez entrer le nombre d'itération: ");
                                scanf("%d", &iteration);

                                float Tab[100][100];
                                int s=0;
                                for(int i=0; i<N; i++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        Tab[i][j] = A[0][s++];
                                    }
                                    //printf("\n");
                                }



                                if(isInvertible(Tab, N))
                                {
                                    matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                    matric_recop(&B[0][0], P, Q, &Brecop[0][0]);
                                    jacobi(&(A[0][0]), &(B[0][0]), N, M, P, Q, &iteration);
                                }else
                                {
                                    printf("\nLa matrice n'est pas inversible. Il y a 0 ou une infinité de solutions. La méthode de Jacobi ne peut pas être utilisée.\n");
                                }



                            }break;
                        case 6:
                            {
                                printf("\n\n----Bienvenue dans la méthode de Gauss-Seidel----\n");
                                matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                matric_recop(&B[0][0], P, Q, &Brecop[0][0]);

                                int iteration;
                                printf("Veuillez entrer le nombre d'itération: ");
                                scanf("%d", &iteration);

                                float Tab[100][100];
                                int s=0;
                                for(int i=0; i<N; i++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        Tab[i][j] = A[0][s++];
                                    }
                                    //printf("\n");
                                }



                                if(isInvertible(Tab, N))
                                {
                                    matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                    matric_recop(&B[0][0], P, Q, &Brecop[0][0]);
                                    gauss_seidel(&(A[0][0]), &(B[0][0]), N, M, P, Q, &iteration);
                                }else
                                {
                                    printf("\nLa matrice n'est pas inversible. Il y a 0 ou une infinité de solutions. La méthode de Gauss-Seidel ne peut pas être utilisée.\n");
                                }




                            }break;

                        case 7:
                            {
                                printf("\n\n----Bienvenue dans la Méthode de Thomas pour les matrices tridiagonles----\n");
                                matric_recop(&A[0][0], N, M, &Arecop[0][0]);
                                matric_recop(&B[0][0], P, Q, &Brecop[0][0]);

                                float Tab[100][100];
                                float Tab1[100];
                                int s=0;
                                for(int i=0; i<N; i++)
                                {
                                    for(int j=0; j<N; j++)
                                    {
                                        Tab[i][j] = A[0][s++];
                                    }
                                    //printf("\n");
                                }
                                for(int i=0, s=0; i<P; i++)
                                {


                                        Tab1[i] = B[0][s++];

                                    //printf("\n");
                                }

                                if(isTridiagonal(Tab, N))
                                {
                                    //printf("Hoyé\n");
                                    thomas_tridiagonal_method(Tab, Tab1, N, M);
                                }else
                                {
                                    printf("\nLa matrice n'est pas tridiagonale, la méthode de Thomas ne peut pas être utilisée\n");
                                }




                            }break;
                        default:
                            {
                                printf("Méthode incorrecte. Veuillez réessayer.\n\n");
                                delay(3);
                                system("cls");
                                continue;

                            }

                   }
                   do
                        {
                            printf("\n\nVoulez-vous recommencer la partie systèmes d'équations algébriques (O/N) ? ");
                            answer = getche();
                            answer = toupper(answer);
                            fflush(stdin);
                        }
                        while((answer!='O')&&(answer!='N'));

                        if(answer=='N')
                        {
                            answer='N';
                        }
                        if(answer=='O')
                        {
                            //system("cls");
                        }



                    }
                    break;
                }
            case 3:
                {
                    char answer='O';
                    int method, nbPt;
                    float X[100], Y[100];
                    printf("\n------------Bienvenue dans la partie interpolation linéaire--------------\n\n");
                    printf("Entrer le nombre de points à interpoler: ");
                    scanf("%d", &nbPt);
                    scan_interpo_point(X, Y, nbPt);

                    while(answer =='O')
                    {
                        fflush(stdin);
                        printf("\n1-Méthode de Lagrange\n\n");
                        printf("2-Méthode de Newton\n\n");
                        printf("3-Approximation au sens des moindres carrés\n\n");
                        printf("Veuillez faire un choix: ");
                        scanf("%d", &method);

                        switch(method)
                        {
                        case 1:
                            {
                                printf("\nBienvenue dans la résolution de Lagrange\n");
                                interpolation_lagrange(&X[0], &Y[0], nbPt);



                            }break;

                        case 2:
                            {
                                printf("\nBienvenue dans la résolution de Newton\n");
                                interpolation_newton(&X[0], &Y[0], nbPt);


                            }break;

                        case 3:
                            {
                                printf("\nBienvenue dans la résolution par la méthode des moindres carrés\n");
                                interpolation_moindres_carree(&X[0], &Y[0], nbPt);

                            }break;

                        default:
                            {
                                printf("Méthode incorrecte. Veuillez réessayer.\n\n");
                                delay(3);
                                system("cls");
                                continue;
                            }
                        }





                        do
                        {
                            printf("\n\nVoulez-vous recommencer la partie interpolation linéaire (O/N) ? ");
                            answer = getche();
                            answer = toupper(answer);
                            fflush(stdin);
                        }
                        while((answer!='O')&&(answer!='N'));

                        if(answer=='N')
                        {
                            answer='N';
                        }
                        if(answer=='O')
                        {
                            //system("cls");
                        }

                    }
                }break;
            case 4:
                {
                    char answer='O';
                    int method;
                    double h, t0, y0;
                    int iteration;
                    double hfinal =0, t0final=0, y0final=0;
                    int iterationfinal=0;
                    printf("\n-----------------Bienvenue dans la partie interpolation linéaire-------------------\n\n");
                    number_check_double("Veuillez entrer un pas de temps h: ", &h);
                    number_check_double("Veuillez entrer la condition initiale t0: ", &t0);
                    number_check_double("y0: ", &y0);
                    number_check_int("Veuillez entrer le nombre d'itérations: ", &iteration);

                    recopie(hfinal, t0final, y0final, iterationfinal, h, t0, y0, iteration);


                    while(answer == 'O')
                    {
                        fflush(stdin);

                        printf("\n1-Méthode d'Euler\n\n");
                        printf("2-Méthode de Runge-Kunta, méthode du point milieu\n\n");
                        printf("3-Méthode de Runge-Kunta, méthode d'Euler modifiée\n\n");
                        printf("Veuillez faire un choix: ");
                        scanf("%d", &method);

                        switch(method)
                        {
                        case 1:
                            {
                                recopie(h, t0, y0, iteration, hfinal, t0final, y0final, iterationfinal);

                                printf("\nBienvenue dans la résolution d'Euler\n");


                                euler(t0, y0, h, iteration);





                            }break;

                        case 2:
                            {
                                recopie(h, t0, y0, iteration, hfinal, t0final, y0final, iterationfinal);

                                printf("\nBienvenue dans la résolution de Runge-Kutta, méthode du point milieu\n");

                                runge_kutta_milieu(t0, y0, h, iteration);

                            }break;

                        case 3:
                            {
                                printf("\nBienvenue dans la résolution de Runge-Kutta, méthode d'Euler modifié\n");
                                recopie(h, t0, y0, iteration, hfinal, t0final, y0final, iterationfinal);
                                runge_kutta_euler(t0, y0, h, iteration);
                            }break;

                        default:
                            {
                                printf("Méthode incorrecte. Veuillez réessayer.\n\n");
                                delay(3);
                                system("cls");
                                continue;
                            }

                        }

                       do
                        {
                            printf("\n\nVoulez-vous recommencer la partie équations différentielles (O/N) ? ");
                            answer = getche();
                            answer = toupper(answer);
                            fflush(stdin);
                        }
                        while((answer!='O')&&(answer!='N'));

                        if(answer=='N')
                        {
                            answer='N';
                        }
                        if(answer=='O')
                        {
                            //system("cls");
                        }
                    }

                }break;

            default:
                {
                    printf("Type de résolution incorrecte. Veuillez réessayer\n\n");
                    delay(3);
                    system("cls");
                    continue;
                }


        }



        //Boucle de reprise du programme entier
        do
            {
                printf("\n\nVoulez-vous recommencer le programme (O/N) ? ");
                reponse = getche();
                reponse = toupper(reponse);
                fflush(stdin);
            }
            while((reponse!='O')&&(reponse!='N'));

            if(reponse=='O')
            {
                //system("cls");
            }

    }
    return 0;
}

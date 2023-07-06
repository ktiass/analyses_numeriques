void delay(int sec);

//Prototype de la proc�dure permettant de v�rifier la loyaut� des nombres entiers entr�s
void number_check_int(char *texte, int *z);

//Prototype de la proc�dure permettant de v�rifier la loyaut� des nombres d�cimaux double entr�s
void number_check_double(char *texte, double *z);

//Prototype de la proc�dure permettant de v�rifier la loyaut� des nombres d�cimaux float entr�s
void number_check_float(char *texte, float *z);

//Prototype de la fonction permttant de calculer un polyn�me de degr� 3
double result_fonction(double x, double *a, double *b, double *c, double *d);

/*Prototype de la fonction qui scanne n'importe quelle polyn�me*/
void polynome(double Tab[], int N);

double calc_polynome(double Tab[], double *x, int N);

void table_initialisation(double Tab[]);

void matric_recop(float *fin, int N, int M, float *debut);

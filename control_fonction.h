void delay(int sec);

//Prototype de la procédure permettant de vérifier la loyauté des nombres entiers entrés
void number_check_int(char *texte, int *z);

//Prototype de la procédure permettant de vérifier la loyauté des nombres décimaux double entrés
void number_check_double(char *texte, double *z);

//Prototype de la procédure permettant de vérifier la loyauté des nombres décimaux float entrés
void number_check_float(char *texte, float *z);

//Prototype de la fonction permttant de calculer un polynôme de degré 3
double result_fonction(double x, double *a, double *b, double *c, double *d);

/*Prototype de la fonction qui scanne n'importe quelle polynôme*/
void polynome(double Tab[], int N);

double calc_polynome(double Tab[], double *x, int N);

void table_initialisation(double Tab[]);

void matric_recop(float *fin, int N, int M, float *debut);

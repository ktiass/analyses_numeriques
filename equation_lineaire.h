void scanMatrice(float * ptr, size_t N, size_t M);


void printMatrice(float * ptr, size_t N, size_t M);


void descente_triangulaire(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void montee_triangulaire(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void elimination_gauss(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void decomposition_crout(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void decomposition_choleski(float *ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void elimination_gauss_jordan(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void inverse_matrice_carree(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void jacobi(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q, int *iteration);


void gauss_seidel(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q, int *iteration);


void elimination_gauss_jordan_advanced(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


float elimination_gauss_jordan_advanced_determinant(float * ptr, float *ptr2, size_t N, size_t M, size_t P, size_t Q);


void thomas_tridiagonal_method(float A[][100], float B[], int N, int M);

void interpolation_lagrange(float *x, float *y, int nbPt);

void interpolation_moindres_carree(float *x, float *y, int nbPt);

void matrice(float ptr[][100], size_t N, size_t M);

void decomposition_choleski_interpo(float ptr[][100], float ptr2[][100], size_t N, size_t M, size_t P, size_t Q);

void interpolation_newton(float *x, float *y, int nbPt);

void scan_interpo_point(float ptr[], float ptr2[], int N);

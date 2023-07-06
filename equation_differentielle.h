void scanFonction(float y[], float t[], int n, int m);

float calcul_fonction(float t0, float y0);

void euler(float t0, float y0, double h, int iteration);

void runge_kutta_milieu(float t0, float y0, double h, int iteration);

void runge_kutta_euler(float t0, float y0, double h, int iteration);

void recopie(double hfinal, double t0final, double y0final, int iterationfinal, double h, double t0, double y0, int iteration);

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void getCofactor(float mat[][100], float temp[][100], int p, int q, int n)
		{
		int i = 0, j = 0;

		// Looping for each
		// element of the matrix
		for (int row = 0; row < n; row++)
		{
			for (int col = 0; col < n; col++)
			{
				// Copying into temporary matrix
				// only those element which are
				// not in given row and column
				if (row != p && col != q)
				{
					temp[i][j++] = mat[row][col];
					//temp[n * i + (j+1)] = mat[n * row + col];

					// Row is filled, so increase
					// row index and reset col index
					if (j == n - 1)
					{
						j = 0;
						i++;
					}
				}
			}
		}
		}

float determinantOfMatrix(float mat[][100],int n)
		{
		float D = 0; // Initialize result

		// Base case : if matrix
		// contains single element
		if (n == 1)
			return mat[0][0];
			//return mat[n * 0 + 0];

		// To store cofactors
		//int [][]temp = new int[N][N];
		float temp[100][100];


		// To store sign multiplier
		int sign = 1;

		// Iterate for each
		// element of first row
		for (int f = 0; f < n; f++)
		{
			// Getting Cofactor of mat[0][f]

			getCofactor(mat, temp, 0, f, n);
			D += sign * mat[0][f] * determinantOfMatrix(temp, n - 1);
			//D += sign * mat[n * 0 + f] * determinantOfMatrix(temp, n - 1);

			// terms are to be added
			// with alternate sign
			sign = -sign;
		}

		return D;
		}
bool isInvertible(float mat[][100], int n)
		{
		    bool iv = false;
		    printf("Le déterminant est %f\n", determinantOfMatrix(mat, n));
			if (determinantOfMatrix(mat, n) != 0)
				{
				    printf("La matrice est inversible\n");
				    iv = true;
				}
			else
				{
				    printf("La matrice n'est pas inversible\n");
				    iv = false;
				}
            return iv;

		}

void cut_matrice(float temp[][100], float mat[][100], int c)
{
    for(int i=0; i<c; i++)
    {
        for(int j=0; j<c; j++)
        {
            temp[i][j] = mat[i][j];
        }
    }
}

bool defined_positive(float mat[][100], int n)
{
    float temp[100][100];
    int i;
    for(i=1; i<=n; i++)
    {
        cut_matrice(temp, mat, i);
        printf("Déterminant = %f", determinantOfMatrix(temp, i));
        if(determinantOfMatrix(temp, i) < 0)
        {
            return false;
        }
    }
    return true;

}

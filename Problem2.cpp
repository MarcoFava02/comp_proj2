#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <string>
#include <armadillo>
#include "Project2.hpp"

int main() {

    //Lets set conditions of the problem
    int N = 6;
    int x_max = 1;
    int x_min = 0;
    double n = N+1; //number of steps
    double h = (x_max-x_min)/n; //stepsize
    double a = -1./(h * h); // sub/super diagonal
    double d = 2./(h * h);  // diagonal

    //Use the headfile to make the tridiagonal matrix
    arma::mat A = make_tridiag_matrix(N,a,d);

    //Lets find the eigenvalues and eigenvectors using the function eig_sym of armadillo
    //Create a vector for the values and a matrix for the vectors without defining the dimensions.
    arma::vec eigenvalues ;
    arma::mat eigenvectors ;

    // Use the function
    arma::eig_sym(eigenvalues, eigenvectors, A);

    //As it is said, we need to scale each eigen vector to the unit norm with the function normalise
    arma::vec eigenvalues_norm = arma::normalise (eigenvalues,1);
    //However we cannot do that with a matrix, we have to pick each column(vector), normalize it and add it to a new matrix

    //Create the normalized matrix which will be filled in the loop. But no specify the dimensions
    //because then when filling the empty values will be a 6*12 matrix.
    arma::mat eigenvectors_norm;
    for (int i = 0; i <= N-1; i++)
    {
        // Pick the full column i of the matrix A
        arma::vec column_vector =  eigenvectors.col(i);
        //Normalize the previous vector
        arma::vec column_vector_norm = arma::normalise (column_vector,1);
        //Add it to the normalized matrix
        eigenvectors_norm .insert_cols(i, column_vector_norm);
    }


    // Analytical solution. Now let fix the dimensions to avoid confusions when compiling
    arma::vec eigenvalues_analyt = arma::vec(N).fill(0.);
    arma::mat eigenvectors_analyt_norm; 
    // El bucle debe ser de 0 a N-1
    for (int j = 0; j <= N-1; j++)
    {
        eigenvalues_analyt(j) = d + 2 * a * std::cos(((j + 1) * arma::datum::pi) / (N + 1));

        // Crea un vector que contiene cada eigenvector para que luego se añada como columna en la matriz final
        arma::vec eigenvectors_analyt_values = arma::vec(N).fill(0.);
        for (int i = 0; i <= N-1; i++)
        {
            eigenvectors_analyt_values(i) = std::sin(((i + 1) * (j + 1) * arma::datum::pi) / (N + 1));
        }
        //Lets normalize each vectorthe vector
        arma::vec eigenvectors_analyt_values_norm = arma::normalise(eigenvectors_analyt_values, 1);
        //Now we are going to replace each define column of the matrix with the vector
        eigenvectors_analyt_norm.insert_cols (j,eigenvectors_analyt_values_norm);
    }


    //Lets normalize the values
    arma::vec eigenvalues_analyt_norm = arma::normalise(eigenvalues_analyt, 1);

    //Shows in the terminal the eigenvalues and the eigenvectors done with the armadillos function and also with a analytical solution
    //They are already scaled/normalized
    std::cout << "The inicial matrix A is: " << std::endl;
    std::cout << A << std::endl;

    std::cout << "The normalized eigenvalues with Armadillo are:" << std::endl;
    std::cout << eigenvalues_norm << std::endl;

    std::cout << "The analytically normalized eigenvalues  are:" << std::endl;
    std::cout << eigenvalues_analyt_norm << std::endl;

    std::cout << "The normalized eigenvectors with Armadillo are:" << std::endl;
    std::cout << eigenvectors_norm << std::endl;

    std::cout << "The analytically normalized eigenvectors are:" << std::endl;
    std::cout << eigenvectors_analyt_norm << std::endl;

}


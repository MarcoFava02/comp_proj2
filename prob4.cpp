#include <iostream>
#include <armadillo>
#include <assert.h>
#include <math.h>
#include <cmath>

#define n 6 // !!this is the size of the matrix, so it is the number of steps minus one (in the text of the problem is called N)!!

#include "funcprob4.h"

int main(){
    //arma::mat A = arma::mat("1.0 0.0 0.0 0.5 0.0; 0.0 1.0 -0.7 0.0 9.0; 0.0 -0.7 1.0 3.0 0.0; 0.5 0.0 3.0 1.0 0.0; 0.0 9.0 0.0 0.0 1.0");
    //arma::mat A = arma::mat("1. 0. 0. 0.5; 0. 1. -0.7 0.; 0. -0.7 1. 0.; 0.5 0. 0. 1.");
    arma::vec eigenvalues(n);
    arma::mat eigenvectors(n,n);
    int iterations = 0;
    int maxiter = 100000;
    bool converged = true;
    double eps = 1e-10; //tolerance value

    double h = 1./(n+1);  // h is the stepsize we use to discretize the problem (since the size of the matrix is n, the number of stepsize is n+1)
    double a = -1/h/h; // a is the value of the elements of the subdiagonals
    double d = 2/h/h;  // d is the value of the elements of the diagonal
    arma::mat A = arma::mat(n,n); //initialize NxN matrix A 

    //check that the size of the matricies are correct
    int m = A.n_rows;               // m is the number of rows of the matrix A
    assert(m == A.n_cols && n>1 && n==m);

    //inizialize the matrix A
    make_tridiag(A, a, d);
    //std::cout << A;
    
    //diagonalize the matrix A to find the eigenvalues and eigenvectors
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    //print the results
    if(converged){
        std::cout << "The algorith converged after " << iterations << " iterations and returned eigenvalues:\n" << eigenvalues << "with eigenvectors:\n" << eigenvectors << std::endl;
    }else if(!converged){
        std::cout << "The algorith did not converged and returned eigenvalues:\n" << eigenvalues << "with eigenvectors:\n" << eigenvectors << std::endl;
    }

    return 0;
}

#ifndef funcprob4
#define funcprob4

#include <iostream>
#include <armadillo>
#include <assert.h>
#include <math.h>
#include <cmath>


// Determine the the max off-diagonal element of a symmetric matrix A
// - Saves the matrix element indicies to k and l 
// - Returns absolute value of A(k,l) as the function return value
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged); //(should A be constant???)

// Function that makes A a tridiagonal matrix with upper and lower diag equal to a and central diag equal to d
void make_tridiag(arma::mat& A, const double a, const double d);

// Function for problem 5 that makes A a tridiagonal matrix with upper and lower diag equal to a and central diag equal to d but other elements the same
void make_tridiag_prob5(arma::mat& A, const double a, const double d);




void make_tridiag(arma::mat& A, const double a, const double d){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j){
                A(i,j) = d;
            } 
            else if(j==i+1 || j==i-1){
                A(i,j) = a;
            } 
            else{
                A(i,j) = 0;
            }
        }
}
}




void make_tridiag_prob5(arma::mat& A, const double a, const double d){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(i==j){
                A(i,j) = d;
            } 
            else if(j==i+1 || j==i-1){
                A(i,j) = a;
            } 
        }
}
}

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    arma::mat M;
    arma::mat P;
    double tau, s, c, t;

    //the starting matricies are A and R but we are going to change the values in the matricies M and P, 
    //then copy them in A and R (as if M=A(m+1) and A=A(m), P=R(m+1) and R=R(m))
    M = A;
    P = R;

    //compute tau,t,s,c
    tau = (A(l,l) - A(k,k))/(2*A(k,l));
    if(tau>=0){
        t = 1/(tau + sqrt(1+pow(tau,2)));
    }else if(tau<0){
        t = -1/(-tau + sqrt(1+pow(tau,2)));
    }
    c = 1/sqrt(1+pow(t,2));
    s = c*t;

    //trasform the new matricies in M and P
    M(k,k) = A(k,k)*pow(c,2) - 2*A(k,l)*c*s + A(l,l)*pow(s,2);
    M(l,l) = A(l,l)*pow(c,2) + 2*A(k,l)*c*s + A(k,k)*pow(s,2);
    M(k,l) = M(l,k) = 0;
    
    for(int i=0; i<n; i++){
        if(i!=k && i!=l){
            M(i,k) = A(i,k)*c - A(i,l)*s;
            M(k,i) = M(i,k);
            M(i,l) = A(i,l)*c + A(i,k)*s;
            M(l,i) = M(i,l);
        }
        P(i,k) = R(i,k)*c - R(i,l)*s;
        P(i,l) = R(i,l)*c + R(i,k)*s;
    }
    
    //copy M and P in A and R
    A = M;
    R = P;

}


void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged){
    arma::mat R = arma::mat(n,n, arma::fill::eye);
    int k, l;
    double a;

    //find the first pair (k,l) for the maximum value of the matrix
    a = max_offdiag_symmetric(A,k,l);

    //start the while loop that stops when all the non diagonal elements are smaller than epsilon 
    //or if the number of iterations reaches maxiter
    while(abs(a)>eps){
        iterations += 1;
        //stop while loop if iterations reaches maxiter
        if(iterations > maxiter){
            converged = false;
            break; 
        }else{
            jacobi_rotate(A,R,k,l);
            a = max_offdiag_symmetric(A,k,l);
        }
    }

    //find and write eigenvalues
    for(int i=0; i<n; i++){
        eigenvalues(i) = A(i,i);
    }

    //the eigenvectors matrix is just the last R matrix
    eigenvectors = R;
}


double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){
    //int n = A.n_rows;               // n is the number of rows of the matrix A
    //assert(n == A.n_cols && n > 1);   // here we check if A is square and larger than 1x1
    // If the conditions are not respected (?), then the program is killed with an assertion error
    
    //cout << A.is.square() << endl; this was an alternative way using armadillo library, although it's incomplete
    //if you want we can delete this two lines of comments

    double maxval;
    maxval = A(0,1);   // first thing first, we start with the first off-diagonal value of A as a 'test' max value
    int i, j;   // k is the index for the rows of the maximum value, l is for columns
    k = i =0;
    l = j = 1;
    
    do{     // while loop to overwrite maxval and the indeces k,l if there's a bigger element in A
        j=j+1;      // we check the condition and we get to the next element in the same row
        if(abs(A(i,j)) > abs(maxval)){ // here we check if the absolute value of element A(i,j) is bigger
            maxval = A(i,j);
            k = i;  // reminder: in armadillo the indices of a n-vector go from 0 to n-1
            l = j;  
        }
        
        if(j==(n-1)) { // if condition to get all the elements of the upper triangle: if we get to the end of the row...
            i=i+1;  // ...then we get to the next row...
            j=i;  // ...with the right index for the column 
        }
    }while(i<n-1); // the algorithm goes until we get to the row i=n-1
    
    return maxval;
}
 

 #endif
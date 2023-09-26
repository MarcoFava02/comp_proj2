#include <iostream>
#include <armadillo>
#include <assert.h>
#include <math.h>
#include <cmath>

#define n 9 // !!this is the size of the matrix, so it is the number of steps minus one (in the text of the problem is called N)!!
#define filename1 "eigenvec1.txt" //the eigenvectors will be saved here
#define filename2 "eigenvec2.txt"
#define filename3 "eigenvec3.txt"

#include "funcprob4.h"



//function to find the three eigenvectors correspondig to the three lowest eigenvalues
void best_eigenvec(const arma::vec& eigenvalues, const arma::mat& eigenvectors, arma::mat& bestvec);


int main(){
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

    //write the three eigenvectors corresponding to the three lowest eigenvalues on a file
    arma::mat bestvec(n,3);
    best_eigenvec(eigenvalues, eigenvectors, bestvec);
    std::ofstream ofile;

    //vector one
    ofile.open(filename1);
    ofile <<0<<","<<0<<std::endl; //also add point (0,0)
    for(int i=1;i<n+1;i++){
        ofile << i*h << "," << bestvec(i-1,0) << std::endl;
    }
    ofile <<1<<","<<0<<std::endl; //also add point (1,0)
    ofile.close();

    //vector two
    ofile.open(filename2);
    ofile <<0<<","<<0<<std::endl;
    for(int i=1;i<n+1;i++){
        ofile << i*h << "," << bestvec(i-1,1) << std::endl;
    }
    ofile <<1<<","<<0<<std::endl;
    ofile.close();

    //vector three
    ofile.open(filename3);
    ofile <<0<<","<<0<<std::endl;
    for(int i=1;i<n+1;i++){
        ofile << i*h << "," << bestvec(i-1,2) << std::endl;
    }
    ofile <<1<<","<<0<<std::endl;
    ofile.close();

    //print the results
    if(converged){
        std::cout << "The algorith converged after " << iterations << " iterations and returned in a file, called "<<filename1<<", the three eigenvectors corresponding to the three lowest eigenvalues." << std::endl;
    }else if(!converged){
        std::cout << "The algorith did not converged and returned eigenvalues:\n" << eigenvalues << "with eigenvectors:\n" << eigenvectors << std::endl;
    }

    return 0;
}




void best_eigenvec(const arma::vec& eigenvalues, const arma::mat& eigenvectors, arma::mat& bestvec){
    arma::vec bestval(3);

    for(int i=0;i<3;i++){
        bestval(i) = eigenvalues(i);
        for(int j=0;j<n;j++){
            bestvec(j,i) = eigenvectors(j,i);
        }
        for(int j=i;j<n;j++){

            switch (i){
            case 0:
                if(eigenvalues(j)<bestval(i)){
                    bestval(i) = eigenvalues(j);
                    for(int k=0;k<n;k++){
                        bestvec(k,i) = eigenvectors(k,j);
                    }
                }
            case 1:
                if(eigenvalues(j)<bestval(i) && bestval(i-1)!=eigenvalues(j)){
                    bestval(i) = eigenvalues(j);
                    for(int k=0;k<n;k++){
                        bestvec(k,i) = eigenvectors(k,j);
                    }
                }
            case 2:
                if(eigenvalues(j)<bestval(i) && bestval(i-1)!=eigenvalues(j) && bestval(i-2)!=eigenvalues(j)){
                    bestval(i) = eigenvalues(j);
                    for(int k=0;k<n;k++){
                        bestvec(k,i) = eigenvectors(k,j);
                    }
                }
            }
        }
    }
    std::cout << "the lowest eigenvalues are: \n" <<bestval<< "and the best eigenvectors are: \n" <<bestvec<< std::endl;
}
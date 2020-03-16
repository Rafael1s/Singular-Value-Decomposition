#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <algorithm>

class svd {
public:
	// ctor
	svd();
	void Factorization(std::vector<std::vector<double>> matrix, std::vector<std::vector<double>>& s,
		std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v);
	void VerifyFactorisition(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& s,
		std::vector<std::vector<double>>& v);
	int GetRank(std::vector<std::vector<double>>& s);
	void GetPseudoinverse(std::vector<std::vector<double>>&  s,
		std::vector<std::vector<double>>&  u, std::vector<std::vector<double>>&  v,
		std::vector<std::vector<double>>&  pseudoinv);
	void VerifyPseudoinverse(std::vector<std::vector<double>>& A, 
		std::vector<std::vector<double>>& A_pesudoinv);

private:
	// M_target is (A^t * A - eigenvalue*I)
	void print_m(std::vector<std::vector<double>> matr, std::string operation);
	void prod_matrix(std::vector<std::vector<double>> A,
		std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& Prod);
	void transpose_matr(std::vector<std::vector<double>> A,
		std::vector<std::vector<double>>& A_t);
	void GaussJordanElimination(std::vector<std::vector<double>> matrix, std::vector<double>& eigen_vec);
	void ConjugateFor_M_target(std::vector<double> eigen_vec, std::vector<std::vector<double>>& h_matrix);
	void InverseConjugateFor_M_target(std::vector<double> eigen_vec, std::vector<std::vector<double>>& ih_matrix);
	void InverseDiagMatrix(std::vector<std::vector<double>> matrix,
		std::vector<std::vector<double>>& inv_matrix);
	void ReduceMatrix(std::vector<std::vector<double>> matrix,
		std::vector<std::vector<double>>& r_matrix, std::size_t new_size);
	void GetEigenvalsEigenvecs(std::vector<std::vector<double>> matrix,
		std::vector<double>& eigenvalues, std::vector<std::vector<double>>& eigen_vecs, std::size_t eig_count);
	void print_eigenv(std::vector<double> eigenv, std::string operation);

private:
	bool track;
    std::vector<std::vector<double>> matrix_i;
};


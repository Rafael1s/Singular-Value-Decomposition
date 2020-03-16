

#include "svd.h"

svd::svd() {
	track = true;
}

void svd::Factorization(std::vector<std::vector<double>> A, std::vector<std::vector<double>>& s,
		std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v)
	{

		std::vector<std::vector<double>> A_t, A_t_x_A, u_1, v_1, s_inverse, av_matrix;
		std::vector<double> eigenvalues;

		print_m(A, "getting matrix");
		transpose_matr(A, A_t);

		prod_matrix(A_t, A, A_t_x_A);
		GetEigenvalsEigenvecs(A_t_x_A, eigenvalues, v_1, 0);

		transpose_matr(v_1, v);

		s.resize(A.size());

		const std::size_t e_size = eigenvalues.size();
		for (std::size_t index = 0; index < e_size; index++)
		{
			s[index].resize(e_size);
			s[index][index] = eigenvalues[index];
		}

		InverseDiagMatrix(s, s_inverse);

		prod_matrix(A, v, av_matrix);
		prod_matrix(av_matrix, s_inverse, u);

		if (track) {
			std::cout << "\nS = \n"; print_m(s, "diag matrix");
			std::cout << "\nU = \n"; print_m(u, "left factor");
			std::cout << "\nV = \n"; print_m(v, "right factor");
		}

	}



// Inverse Diagonal Matrix that is not necessary non-degenerate
void svd::InverseDiagMatrix(std::vector<std::vector<double>> A,
		std::vector<std::vector<double>>& inv_A)
	{
	    const std::size_t m_size = A.size();

		inv_A.resize(m_size);

		for (std::size_t index = 0; index < m_size; index++)
			inv_A[index].resize(A[index].size());

		for (std::size_t index = 0; index < m_size; index++) {
			if (A[index][index] != 0.) {
				inv_A[index][index] = 1.0/A[index][index];
			}
			else {
				inv_A[index][index] = 0.;
			}
		}
	}


void svd::prod_matrix(std::vector<std::vector<double>> A,
		std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& prod){
		prod.resize(A.size());
		for (std::size_t row = 0; row < A.size(); row++){
			prod[row].resize(A[row].size());
			std::fill(prod[row].begin(), prod[row].end(), 0);
		}

		std::size_t m_size = A.size();

		for (std::size_t row = 0; row < m_size; row++)
			for (std::size_t col = 0; col < m_size; col++){
				for (std::size_t k = 0; k < m_size; k++)
					prod[row][col] += A[row][k] * B[k][col];
			}
}

void svd::transpose_matr(std::vector<std::vector<double>> A,
		std::vector<std::vector<double>>& A_t){

	std::size_t m_size = A.size();

	A_t.resize(m_size);
	for (auto& m : A)
		m.resize(m_size);

	for (std::size_t row = 0; row < m_size; row++)
			A_t[row].resize(A[row].size());


	for (std::size_t row = 0; row < m_size; row++)
		for (std::size_t col = 0; col < m_size; col++)
			A_t[row][col] = A[col][row];

	// print_m(A_t, "transpose");
}

int RandomNumber() { return (std::rand() % 10000); }

// Recursive computation of eigenvalues and eigen-vectors
// Compute eigenvalue[i] and eigen_vecs[i], where i = eig_count
void svd::GetEigenvalsEigenvecs(std::vector<std::vector<double>> matrix,
	std::vector<double>& eigenvalues, 
	std::vector<std::vector<double>>& eigen_vecs, 
	std::size_t eig_count){
	
    	const std::size_t m_size = matrix.size();
		std::vector<double> vec; vec.resize(m_size);
		std::generate(vec.begin(), vec.end(), RandomNumber);

		if (eigenvalues.size() == 0 && eigen_vecs.size() == 0)
		{
			eigenvalues.resize(m_size);
			eigen_vecs.resize(eigenvalues.size());
			matrix_i = matrix;
		}

		std::vector<double> m; m.resize(m_size);
		std::vector<double> m_temp; m_temp.resize(m_size);

		double lambda_old = 0;

		// Power Iteration algoritm for finding eigenvalues of the symmetric matrix A^t*A  
		// The solution is obtained as  (A^t*A)^\infinity * arbitrary vector

		std::size_t index = 0; bool is_eval = false;
		while (is_eval == false)
		{
			// m - will be eigenvector
			for (std::size_t row = 0; row < m_size; row++)
				m[row] = 0.f;

			for (std::size_t row = 0; row < m_size; row++)
			{
				for (std::size_t col = 0; col < m_size; col++)
					m[row] += matrix[row][col] * vec[col];
			}

			for (std::size_t col = 0; col < m_size; col++)
				vec[col] = m[col];

			if (index > 0)
			{
				// finish compute eigenvalue if lambda almost const
				double lambda = (index > 0) ? (m[0] / m_temp[0]) : m[0];
				is_eval = (fabs(lambda - lambda_old) <= 10e-10) ? true : false;

				eigenvalues[eig_count] = lambda;
				lambda_old = lambda;
			}

			for (std::size_t row = 0; row < m_size; row++)
				m_temp[row] = m[row];

			index++;
		}

		if (track)
			print_eigenv(eigenvalues, "is_eval == true");

		std::vector<std::vector<double>> matrix_new;

		if (m_size > 1)
		{
			std::vector<std::vector<double>> M_target;
			M_target.resize(m_size);

			for (std::size_t row = 0; row < m_size; row++)
				M_target[row].resize(m_size);

			// M_target is (A^t * A - eigval*I)
			for (std::size_t row = 0; row < m_size; row++)
				for (std::size_t col = 0; col < m_size; col++)
					M_target[row][col] = (row == col) ? \
					(matrix[row][col] - eigenvalues[eig_count]) : matrix[row][col];

			// Get eigen_vecs[i]
			std::vector<double> eigen_vec;
			GaussJordanElimination(M_target, eigen_vec);

			//  Matrix H - Conjugate for matrix A^t*A
			std::vector<std::vector<double>> H;
			ConjugateFor_M_target(eigen_vec, H);

			std::vector<std::vector<double>> H_matrix_prod;
			prod_matrix(H, matrix, H_matrix_prod);

			// inverse matrix for H: inv_H
			std::vector<std::vector<double>> inv_H;
			InverseConjugateFor_M_target(eigen_vec, inv_H);

			std::vector<std::vector<double>> inv_H_matrix_prod;

			// Here, we get inv_H_matrix_prod = H * M * H^-1
			prod_matrix(H_matrix_prod, inv_H, inv_H_matrix_prod);

			// matrix_new = H * M * H^-1 is equalent to  A_t_x_A without first row and first col
			ReduceMatrix(inv_H_matrix_prod, matrix_new, m_size - 1);

		}

		if (m_size <= 1)
		{
			for (std::size_t index = 0; index < eigenvalues.size(); index++)
			{
				double lambda = eigenvalues[index];

				// M_target is (A^t * A - eigval*I)
				std::vector<std::vector<double>> M_target;
				M_target.resize(matrix_i.size());

				for (std::size_t row = 0; row < matrix_i.size(); row++)
					M_target[row].resize(matrix_i.size());

				std::size_t mi_size = matrix_i.size();

				for (std::size_t row = 0; row < mi_size; row++)
					for (std::size_t col = 0; col < mi_size; col++)
						M_target[row][col] = (row == col) ? \
						(matrix_i[row][col] - lambda) : matrix_i[row][col];

				eigen_vecs.resize(matrix_i.size());
				GaussJordanElimination(M_target, eigen_vecs[index]);

				// Normalize eigen vectors
				double eigsum_sq = 0;

				for (std::size_t v = 0; v < eigen_vecs[index].size(); v++)
					eigsum_sq += std::pow(eigen_vecs[index][v], 2.0);

				for (std::size_t v = 0; v < eigen_vecs[index].size(); v++)
					eigen_vecs[index][v] /= std::sqrt(eigsum_sq);

				// Essentially eigenvalues[index] should be positive for symmetric matrix,
				// However, negative values or very small values 
				// may occur due to the accuracy of the calculations.
				if (eigenvalues[index] < 10e-5)
					eigenvalues[index] = 0;

				eigenvalues[index] = std::sqrt(eigenvalues[index]);
			}

			return;
		}

		GetEigenvalsEigenvecs(matrix_new, eigenvalues, eigen_vecs, eig_count + 1);

		return;
}

// Discard first row and first column
void svd::ReduceMatrix(std::vector<std::vector<double>> matrix,
		std::vector<std::vector<double>>& reduced_matr, std::size_t new_size){
		if (new_size > 1)
		{
			reduced_matr.resize(new_size);
			std::size_t index_d = matrix.size() - new_size;
			std::size_t row = index_d, row_n = 0;
			while (row < matrix.size())
			{
				reduced_matr[row_n].resize(new_size);
				std::size_t col = index_d, col_n = 0;
				while (col < matrix.size())
					reduced_matr[row_n][col_n++] = matrix[row][col++];

				row++; row_n++;
			}
		}

		else if (new_size == 1)
		{
			reduced_matr.resize(new_size);
			reduced_matr[0].resize(new_size);
			reduced_matr[0][0] = matrix[1][1];
		}
}

void svd::ConjugateFor_M_target(std::vector<double> eigen_vec,
		std::vector<std::vector<double>>& h_matrix){
		h_matrix.resize(eigen_vec.size());

		for (std::size_t row = 0; row < eigen_vec.size(); row++)
			h_matrix[row].resize(eigen_vec.size());

		h_matrix[0][0] = 1.0 / eigen_vec[0];

		for (std::size_t row = 1; row < eigen_vec.size(); row++)
			h_matrix[row][0] = -eigen_vec[row] / eigen_vec[0];

		for (std::size_t row = 1; row < eigen_vec.size(); row++)
			h_matrix[row][row] = 1;
}

void svd::InverseConjugateFor_M_target(std::vector<double> eigen_vec, std::vector<std::vector<double>>& ih_matrix){
		ih_matrix.resize(eigen_vec.size());
		for (std::size_t row = 0; row < eigen_vec.size(); row++)
			ih_matrix[row].resize(eigen_vec.size());

		ih_matrix[0][0] = eigen_vec[0];
		for (std::size_t row = 1; row < eigen_vec.size(); row++)
			ih_matrix[row][0] = -eigen_vec[row];

		for (std::size_t row = 1; row < eigen_vec.size(); row++)
			ih_matrix[row][row] = 1;
}


// Gauss Jordan elimination algorithm --> reduce the matrix to  row echelon form.
// 1 0 0 v1 
// 0 1 0 v2
// 0 0 1 v3
// 0 0 0  0
// The solution of Av 0 is (-v1, -v2, -v3, 1)^t
// The eigenvector v is optaibed up to factor, so the last element may be set to 1

// Get eigen_vector for given M = A^t*A*v - \lambda*v
void svd::GaussJordanElimination(std::vector<std::vector<double>> M, std::vector<double>& eigen_vec) {

		for (std::size_t s = 0; s < M.size() - 1; s++)
		{
			const double diag_elem = M[s][s];
			if (diag_elem != 0 && diag_elem != 1)
			{
				for (std::size_t col = s; col < M[s].size(); col++)
					M[s][col] /= diag_elem;
			}


			// Move the column with element M[s][s] to the end of matrix
			if (diag_elem == 0)
			{
				for (std::size_t col = s; col < M[s].size(); col++)
					std::swap(M[s][col], M[s + 1][col]);
			}

			// Gauss–Jordan elimination process
			for (std::size_t row = 0; row < M.size(); row++)
			{
				const double element = M[row][s];
				if (row != s)
				{
					// Eliminate element M[row][s]. Subtract the line 's' from line 'row != s' 
					for (std::size_t col = s; col < M[row].size(); col++)
						M[row][col] = M[row][col] - M[s][col] * element;
				}
			}
		}

		std::size_t row = 0;
		while (row < M.size())
			eigen_vec.push_back(-M[row++][M.size() - 1]);

		eigen_vec[eigen_vec.size() - 1] = 1;
}

void svd::VerifyFactorisition(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& s,
	std::vector<std::vector<double>>& v) {

	std::vector<std::vector<double>> prod_US, prod_USVt, v_transp;

	prod_matrix(u, s, prod_US);
	transpose_matr(v, v_transp);
	prod_matrix(prod_US, v_transp, prod_USVt);

	if (track)
	     print_m(prod_USVt, "VerifyFactorisition: Product usv^t should be equal to initial matrix");

}

void svd::VerifyPseudoinverse(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& A_pesudoinv) {

	std::vector<std::vector<double>> prod;

	// Product should be Identity matrix or with several zeros on the diagaonal  
	prod_matrix(A, A_pesudoinv, prod);

	if (track)
		print_m(prod, "VerifyPseudoinverse: Product A*A_pinv should be initial matrix");
}


int svd::GetRank(std::vector<std::vector<double>>& s) {
	int rank = 0;

	for (std::size_t row = 0; row < s.size(); row++)
		if (s[row][row] != 0) {
			rank++;
		}

	return rank;
}

void svd::GetPseudoinverse(std::vector<std::vector<double>>&  s, 
	std::vector<std::vector<double>>&  u, std::vector<std::vector<double>>&  v,
	std::vector<std::vector<double>>&  pseudoinv){

	std::vector<std::vector<double>> s_inv, prod_VS_inv, u_transp;
	InverseDiagMatrix(s, s_inv);

	prod_matrix(v, s_inv, prod_VS_inv);

	transpose_matr(u, u_transp);
	prod_matrix(prod_VS_inv, u_transp, pseudoinv);

	int rank = GetRank(s);

	if(rank < (int)s.size())
	    print_m(pseudoinv, "Pseudoinverse Matrix");
	else
		print_m(pseudoinv, "Inverse Matrix");

}

void svd::print_m(std::vector<std::vector<double>>	M, std::string operation){

	    std::cout <<  "print_m: " << operation.c_str() << ": " << std::endl;

		for (std::size_t row = 0; row < M.size(); row++)
		{
			for (std::size_t col = 0; col < M[row].size(); col++)
				std::cout << std::fixed << M[row][col] << " ";

			std::cout << "\n";
		}

		std::cout << "\n";
}

void svd::print_eigenv(std::vector<double> eigenv, std::string operation) {

	std::cout << "print_eigenv: " << operation.c_str() << ": " << std::endl;

	for (std::size_t row = 0; row < eigenv.size(); row++)
	{
			std::cout << std::fixed << eigenv[row] << " ";
	}

	std::cout << "\n";
}

#include "svd.h"
#include "read_matrix.h"

int main(int argc, char* argv[])
{
	std::string rep = "y";
	std::size_t matrix_size = 0;
	std::vector<std::vector<double>> u, v;
	std::vector<std::vector<double>> matrix, s, pseudoinverse;

	read_matrix readm;
	svd fsvd;
	std::cout << "Singular Value Decomposition (SVD):\n\n";

	while (rep != "n" && rep != "N") {

		std::string filename = readm.get_file_name();
		readm.get_matrix_by_file_name(filename, matrix);

		fsvd.Factorization(matrix, s, u, v);

		// Verification, print u*s*v^t. It should be = matrix
		fsvd.VerifyFactorisition(u, s, v);

		int rank = fsvd.GetRank(s);
		std::cout << "Rank: " << rank << std::endl;

		fsvd.GetPseudoinverse(s, u, v, pseudoinverse);

		// The product should be Identical/ or with several zeros on the diagaonal  
		fsvd.VerifyPseudoinverse(matrix, pseudoinverse);

		std::cout << "Continue (Y/N) ?" << std::endl;
		std::cin >> rep;
	}

	return 0;
}

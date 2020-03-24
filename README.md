## Singular Value Decomposition

This package contains 2 C++ classes: **svd**  and **read_matrix**  and test file with the main function.

### Class svd

This class provides factorization of any square matrix **M** by **SVD algorithm** for   
the square matrices:   _M = USV^t_,  where _S_ is the diagonal matrix, and _U_, _V_ are two orthonormal matrices.  

 This class contains 5 public functions   
   
_void **Factorization**(std::vector<std::vector<double>> matrix, std::vector<std::vector<double>>& s,
                    std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v);_
    
_void **VerifyFactorisition**(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& s,
		std::vector<std::vector<double>>& v);_   
    
_int **GetRank**(std::vector<std::vector<double>>& s);_

_void **GetPseudoinverse**(std::vector<std::vector<double>>&  s,
		std::vector<std::vector<double>>&  u, std::vector<std::vector<double>>&  v,
		std::vector<std::vector<double>>&  pseudoinv);_
    
_void **VerifyPseudoinverse**(std::vector<std::vector<double>>& A, 
		std::vector<std::vector<double>>& A_pesudoinv);_

and a number of private functions, see **svd.h**.

### Tracking

The constructor svd contains the boolean variable _track_.    
The default value of _track_ is _true_ that provides the extended log on the console.     
Set _track = false_ to reduce the log.     


### Class read_matrix

This class contains 2 public functions:
   
_std::string **get_file_name();**_

_void **get_matrix_by_file_name**(
      std::string filename, std::vector<std::vector<double>>& matrix);_
	
and several private functions, see read_matrix.h

### Briefly about SVD algorithm

#### 1. Eigenvalues of the symmetric matrix (Power Iteration algorithm)

Let _A_ be the input square matrix. First we find eigenvalues of the symmetric matrix   
_B = A^t * A_. The first loop in the function _GetEigenvalsEigenvecs_ finds   
the maximal eigenvalue of _B_ (in magnitude).  For explanations of this algorithm, see

1. [qr-algorithm](http://madrury.github.io/jekyll/update/statistics/2017/10/04/qr-algorithm.html),
2. [Power itaration](https://en.wikipedia.org/wiki/Power_iteration)

The function _GetEigenvalsEigenvecs_ is recursive, so the remaining eigenvalues will be calculated
by the following recursive calls.  The function _ReduceMatrix_ transforms the matrix _B_ to the matrix
_B_ of the m_size - 1, where m_size is the size if current matrix _B_.   
   
#### 2. Gauss-Jordan Elimination algorithm

The matrix _M = B - \lambda*I_ is named the traget matrix. For each eigenvalue \lambda,
the target matrix _B_ -  is transformed into an upper triangular matrix, _so called row ecehlon form_,   seee

3. [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination)

see function _GaussJordanElimination_.

#### 3. Eigenalues of target matrix 

Further, for each eigenvalue \lambda,  the eigenvector vector V(\lambda) of the target matrix _M_ will be found. 
This is performed by Gauss-Jordan Elimination algorithm, 



	








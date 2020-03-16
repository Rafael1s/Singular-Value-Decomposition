## Singular Value Decomposition

This package contains 2 C++ classes: **svd**  and **read_matrix**  and test file with the main function.

### Class svd

This class provides factorization of any square matrix **M** by **SVD algorithm** for the square matrices:   

          M = USV^t,    
where _S_ is the diagonal matrix, and _U_, _V_ are two orthonormal matrices.  

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

**Tracking**

The constructor svd contains the boolean variable _track_.    
The default value of _track_ is _true_ that provides the extended log on the console.     
Set _track = false_ to reduce the log.     







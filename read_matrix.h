#pragma once

#include <string>
#include <vector>

class read_matrix {
public: 
	
	std::string get_file_name(); // insert filename without extension
	void get_matrix_by_file_name(std::string filename, std::vector<std::vector<double>>& matrix);

private:
	// This is the input matrix
	std::vector<std::vector<double>> A; 

	void load(const std::string& input_stream);
	void get_matrix(std::vector<std::vector<double>>& matrix);
};

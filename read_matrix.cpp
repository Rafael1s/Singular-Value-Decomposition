#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <string>

#include "read_matrix.h"

void read_matrix::get_matrix_by_file_name(
	std::string filename,
	std::vector<std::vector<double>>& matrix) {

	load(filename);
	get_matrix(matrix);
}

std::string read_matrix::get_file_name() {

	std::string file_raw, filename;
	
	// The file name shoul be with extension '.txt'
	// Input file name without extension
	std::cout << "Enter name of file containg matrix: ";
	std::cin >> file_raw;
	filename = file_raw.append(".txt");

	std::cout << "Filname: " << filename << std::endl;
	return filename;
}

void read_matrix::load(const std::string& input_stream)
{
	std::ifstream file(input_stream);
	if (!file)
	{
		std::cerr << "Error opening matrix file.\n";
		return;
	}

	int row, col;
	file >> row; 
	col = row;

	std::cout << "row: " << row << std::endl;

	if (row < 1) 
	{
		std::cerr << "Matrix sizes are out of bounds.\n";
		return;
	}

	// Size the vector using the values read from the file.
	// std::vector<std::vector<double>>  A;
	A.resize(row);
	for (auto& m : A)
		m.resize(col);

	// Read the input file.
	for (auto& row : A)
		for (auto& inner : row)
			file >> inner;

}

void read_matrix::get_matrix(std::vector<std::vector<double>>& out_A) {
	
	const std::size_t m_size = A.size();

	out_A.resize(m_size);
	for (auto& m : out_A)
		m.resize(m_size);

	
	std::cout << "get_matrix: " << std::endl;

	for (std::size_t row = 0; row < m_size; row++)
		for (std::size_t col = 0; col < m_size; col++) {
				out_A[row][col] = A[row][col];
		}
    
}
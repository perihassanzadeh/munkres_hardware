/*
 *   Copyright (c) 2007 John Weaver
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 */

/*
 * Some example code.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
#include <limits>

constexpr int NORMAL = 0;
constexpr int STAR   = 1;
constexpr int PRIME  = 2;
constexpr size_t MAX_SIZE = 256;

// Function declarations
void replace_infinites(double matrix[][MAX_SIZE], size_t size);
void minimize_along_direction(double matrix[][MAX_SIZE], size_t size, bool over_columns);
bool find_uncovered_in_matrix(const double matrix[][MAX_SIZE], const int mask_matrix[][MAX_SIZE],
                            const bool row_mask[], const bool col_mask[],
                            const double item, size_t &row, size_t &col, size_t size);
int step1(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size);
int step2(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size);
int step3(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size,
          size_t &saverow, size_t &savecol);
int step4(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size,
          size_t saverow, size_t savecol);
int step5(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size);

// Function implementations
void replace_infinites(double matrix[][MAX_SIZE], size_t size) {
    double max = matrix[0][0];
    constexpr auto infinity = std::numeric_limits<double>::infinity();

    // Find the greatest value in the matrix that isn't infinity
    for (size_t row = 0; row < size; row++) {
        for (size_t col = 0; col < size; col++) {
            if (matrix[row][col] != infinity) {
                if (max == infinity) {
                    max = matrix[row][col];
                } else {
                    max = std::max<double>(max, matrix[row][col]);
                }
            }
        }
    }

    // A value higher than the maximum value present in the matrix
    if (max == infinity) {
        max = 0;
    } else {
        max++;
    }

    for (size_t row = 0; row < size; row++) {
        for (size_t col = 0; col < size; col++) {
            if (matrix[row][col] == infinity) {
                matrix[row][col] = max;
            }
        }
    }
}

void minimize_along_direction(double matrix[][MAX_SIZE], size_t size, bool over_columns) {
    for (size_t i = 0; i < size; i++) {
        double min = over_columns ? matrix[0][i] : matrix[i][0];

        for (size_t j = 1; j < size && min > 0; j++) {
            min = std::min<double>(
                min,
                over_columns ? matrix[j][i] : matrix[i][j]);
        }

        if (min > 0) {
            for (size_t j = 0; j < size; j++) {
                if (over_columns) {
                    matrix[j][i] -= min;
                } else {
                    matrix[i][j] -= min;
                }
            }
        }
    }
}

bool find_uncovered_in_matrix(const double matrix[][MAX_SIZE], const int mask_matrix[][MAX_SIZE],
                            const bool row_mask[], const bool col_mask[],
                            const double item, size_t &row, size_t &col, size_t size) {
    for (row = 0; row < size; row++) {
        if (!row_mask[row]) {
            for (col = 0; col < size; col++) {
                if (!col_mask[col]) {
                    if (matrix[row][col] == item) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int step1(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size) {
    for (size_t row = 0; row < size; row++) {
        bool found_star = false;
        for (size_t col = 0; col < size && !found_star; col++) {
            if (matrix[row][col] == 0) {
                bool has_star = false;
                for (size_t nrow = 0; nrow < row && !has_star; nrow++) {
                    if (mask_matrix[nrow][col] == STAR) {
                        has_star = true;
                    }
                }
                if (!has_star) {
                    mask_matrix[row][col] = STAR;
                    found_star = true;
                }
            }
        }
    }
    return 2;
}

int step2(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size) {
    size_t covercount = 0;
    for (size_t row = 0; row < size; row++) {
        for (size_t col = 0; col < size; col++) {
            if (mask_matrix[row][col] == STAR) {
                col_mask[col] = true;
                covercount++;
            }
        }
    }

    if (covercount >= size) {
        return 0;
    }
    return 3;
}

int step3(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size,
          size_t &saverow, size_t &savecol) {
    if (find_uncovered_in_matrix(matrix, mask_matrix, row_mask, col_mask, 0, saverow, savecol, size)) {
        mask_matrix[saverow][savecol] = PRIME;
    } else {
        return 5;
    }

    for (size_t ncol = 0; ncol < size; ncol++) {
        if (mask_matrix[saverow][ncol] == STAR) {
            row_mask[saverow] = true;
            col_mask[ncol] = false;
            return 3;
        }
    }

    return 4;
}

int step4(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size,
          size_t saverow, size_t savecol) {
    // Find alternating sequence of starred and primed zeros
    size_t row = saverow;
    size_t col = savecol;
    bool found_star = false;
    bool found_prime = false;

    // Find starred zero in column
    for (size_t i = 0; i < size; i++) {
        if (mask_matrix[i][col] == STAR) {
            row = i;
            found_star = true;
            break;
        }
    }

    if (found_star) {
        // Find primed zero in row
        for (size_t j = 0; j < size; j++) {
            if (mask_matrix[row][j] == PRIME) {
                col = j;
                found_prime = true;
                break;
            }
        }
    }

    // Update masks
    if (found_star && found_prime) {
        mask_matrix[saverow][savecol] = STAR;
        mask_matrix[row][col] = NORMAL;
    }

    // Clear primes and reset masks
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            if (mask_matrix[i][j] == PRIME) {
                mask_matrix[i][j] = NORMAL;
            }
        }
    }

    for (size_t i = 0; i < size; i++) {
        row_mask[i] = false;
        col_mask[i] = false;
    }

    return 2;
}

int step5(double matrix[][MAX_SIZE], int mask_matrix[][MAX_SIZE],
          bool row_mask[], bool col_mask[], size_t size) {
    double h = std::numeric_limits<double>::max();
    
    // Find minimum uncovered value
    for (size_t row = 0; row < size; row++) {
        if (!row_mask[row]) {
            for (size_t col = 0; col < size; col++) {
                if (!col_mask[col]) {
                    if (h > matrix[row][col] && matrix[row][col] != 0) {
                        h = matrix[row][col];
                    }
                }
            }
        }
    }

    // Add h to covered rows
    for (size_t row = 0; row < size; row++) {
        if (row_mask[row]) {
            for (size_t col = 0; col < size; col++) {
                matrix[row][col] += h;
            }
        }
    }

    // Subtract h from uncovered columns
    for (size_t col = 0; col < size; col++) {
        if (!col_mask[col]) {
            for (size_t row = 0; row < size; row++) {
                matrix[row][col] -= h;
            }
        }
    }

    return 3;
}

void solve_munkres(double* matrix, size_t rows, size_t columns) {
    const size_t size = std::max(rows, columns);
    
    // Local arrays for processing
    double local_matrix[MAX_SIZE][MAX_SIZE];
    int mask_matrix[MAX_SIZE][MAX_SIZE];
    bool row_mask[MAX_SIZE];
    bool col_mask[MAX_SIZE];
    size_t saverow = 0, savecol = 0;
    
    // Copy input matrix to local array
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            local_matrix[i][j] = matrix[i * columns + j];
        }
    }

    // Initialize masks
    for (size_t i = 0; i < size; i++) {
        row_mask[i] = false;
        col_mask[i] = false;
    }

    // Prepare matrix values
    replace_infinites(local_matrix, size);
    minimize_along_direction(local_matrix, size, rows >= columns);
    minimize_along_direction(local_matrix, size, rows < columns);

    // Follow the steps
    int step = 1;
    while (step) {
        switch (step) {
            case 1:
                step = step1(local_matrix, mask_matrix, row_mask, col_mask, size);
                break;
            case 2:
                step = step2(local_matrix, mask_matrix, row_mask, col_mask, size);
                break;
            case 3:
                step = step3(local_matrix, mask_matrix, row_mask, col_mask, size, saverow, savecol);
                break;
            case 4:
                step = step4(local_matrix, mask_matrix, row_mask, col_mask, size, saverow, savecol);
                break;
            case 5:
                step = step5(local_matrix, mask_matrix, row_mask, col_mask, size);
                break;
        }
    }

    // Store results back to input matrix
    for (size_t row = 0; row < rows; row++) {
        for (size_t col = 0; col < columns; col++) {
            if (mask_matrix[row][col] == STAR) {
                matrix[row * columns + col] = 0;
            } else {
                matrix[row * columns + col] = -1;
            }
        }
    }
}

int main(int argc, char *argv[]) 
{
	int nrows = 101;
	int ncols = 101;
	
	if ( argc == 3 ) {
		nrows = atoi(argv[1]);
		ncols = atoi(argv[2]);
	}

	std::ifstream myfile("10x10_assignment.txt");
	if (!myfile.is_open()) {
		std::cerr << "Failed to open file" << std::endl;
		return 1;
	}

	std::string temp;
	int rowz, colz;
	int line = 1;
	std::vector<double> vec;

	// Read dimensions
	if (!std::getline(myfile, temp)) {
		std::cerr << "Failed to read number of rows" << std::endl;
		return 1;
	}
	rowz = std::stoi(temp);

	if (!std::getline(myfile, temp)) {
		std::cerr << "Failed to read number of columns" << std::endl;
		return 1;
	}
	colz = std::stoi(temp);

	// Read matrix data
	while (std::getline(myfile, temp)) {
		std::stringstream ss(temp);
		std::string ele;
		while (ss >> ele) {
			vec.push_back(std::stod(ele));
		}
	}

	// Verify we have enough elements
	if (vec.size() != rowz * colz) {
		std::cerr << "Matrix dimensions don't match data size" << std::endl;
		return 1;
	}

	// Create a flat array for the matrix
	double* matrix = new double[rowz * colz];
	int item = 0;
	for (int i = 0; i < rowz; i++) {
		for (int j = 0; j < colz; j++) {
			matrix[i * colz + j] = vec[item];
			std::cout << vec[item] << " ";
			item++;
		}
		std::cout << std::endl;
	}

	// Display begin matrix state.
	std::cout << "\nInitial matrix state:" << std::endl;
	for (int row = 0; row < rowz; row++) {
		for (int col = 0; col < colz; col++) {
			std::cout.width(2);
			std::cout << matrix[row * colz + col] << ",";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	// Apply Munkres algorithm to matrix.
	solve_munkres(matrix, rowz, colz);

	// Display solved matrix.
	std::cout << "Solved matrix state:" << std::endl;
	for (int row = 0; row < rowz; row++) {
		for (int col = 0; col < colz; col++) {
			std::cout.width(2);
			std::cout << matrix[row * colz + col] << ",";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;

	// Verify solution
	for (int row = 0; row < rowz; row++) {
		int rowcount = 0;
		for (int col = 0; col < colz; col++) {
			if (matrix[row * colz + col] == 0)
				rowcount++;
		}
		if (rowcount != 1)
			std::cerr << "Row " << row << " has " << rowcount << " columns that have been matched." << std::endl;
	}

	for (int col = 0; col < colz; col++) {
		int colcount = 0;
		for (int row = 0; row < rowz; row++) {
			if (matrix[row * colz + col] == 0)
				colcount++;
		}
		if (colcount != 1)
			std::cerr << "Column " << col << " has " << colcount << " rows that have been matched." << std::endl;
	}

	delete[] matrix;
	return 0;
}

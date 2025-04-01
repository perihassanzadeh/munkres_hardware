//Peri Hassanzadeh Munkres Algorithm Modified for Hardware Acceleration 

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

double global_matrix[MAX_SIZE][MAX_SIZE];
int global_mask_matrix[MAX_SIZE][MAX_SIZE];
bool global_row_mask[MAX_SIZE];
bool global_col_mask[MAX_SIZE];
size_t global_size;
size_t global_saverow;
size_t global_savecol;
double global_item; 

// Function declarations
void replace_infinites(void);
void minimize_along_direction(bool over_columns);
bool find_uncovered_in_matrix(size_t &row, size_t &col);  
int step1(void);
int step2(void);
int step3(void);
int step4(void);
int step5(void);


void replace_infinites(void) {
    double max = global_matrix[0][0];
    constexpr auto infinity = std::numeric_limits<double>::infinity();

    // Find the greatest value in the matrix that isn't infinity
    for (size_t row = 0; row < global_size; row++) {
        for (size_t col = 0; col < global_size; col++) {
            if (global_matrix[row][col] != infinity) {
                if (max == infinity) {
                    max = global_matrix[row][col];
                } else {
                    max = std::max<double>(max, global_matrix[row][col]);
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

    for (size_t row = 0; row < global_size; row++) {
        for (size_t col = 0; col < global_size; col++) {
            if (global_matrix[row][col] == infinity) {
                global_matrix[row][col] = max;
            }
        }
    }
}

void minimize_along_direction(bool over_columns) {
    for (size_t i = 0; i < global_size; i++) {
        double min = over_columns ? global_matrix[0][i] : global_matrix[i][0];

        for (size_t j = 1; j < global_size && min > 0; j++) {
            min = std::min<double>(
                min,
                over_columns ? global_matrix[j][i] : global_matrix[i][j]);
        }

        if (min > 0) {
            for (size_t j = 0; j < global_size; j++) {
                if (over_columns) {
                    global_matrix[j][i] -= min;
                } else {
                    global_matrix[i][j] -= min;
                }
            }
        }
    }
}

bool find_uncovered_in_matrix(size_t &row, size_t &col) {
    for (row = 0; row < global_size; row++) {
        if (!global_row_mask[row]) {
            for (col = 0; col < global_size; col++) {
                if (!global_col_mask[col]) {
                    if (global_matrix[row][col] == global_item) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

int step1(void) {
    for (size_t row = 0; row < global_size; row++) {
        bool found_star = false;
        for (size_t col = 0; col < global_size && !found_star; col++) {
            if (global_matrix[row][col] == 0) {
                bool has_star = false;
                for (size_t nrow = 0; nrow < row && !has_star; nrow++) {
                    if (global_mask_matrix[nrow][col] == STAR) {
                        has_star = true;
                    }
                }
                if (!has_star) {
                    global_mask_matrix[row][col] = STAR;
                    found_star = true;
                }
            }
        }
    }
    return 2;
}

int step2(void) {
    size_t covercount = 0;
    for (size_t row = 0; row < global_size; row++) {
        for (size_t col = 0; col < global_size; col++) {
            if (global_mask_matrix[row][col] == STAR) {
                global_col_mask[col] = true;
                covercount++;
            }
        }
    }

    if (covercount >= global_size) {
        return 0;
    }
    return 3;
}

int step3(void) {
    global_item = 0;  // Set the global item before calling find_uncovered_in_matrix
    if (find_uncovered_in_matrix(global_saverow, global_savecol)) {
        global_mask_matrix[global_saverow][global_savecol] = PRIME;
    } else {
        return 5;
    }

    for (size_t ncol = 0; ncol < global_size; ncol++) {
        if (global_mask_matrix[global_saverow][ncol] == STAR) {
            global_row_mask[global_saverow] = true;
            global_col_mask[ncol] = false;
            return 3;
        }
    }

    return 4;
}

int step4(void) {
    // Find alternating sequence of starred and primed zeros
    size_t row = global_saverow;
    size_t col = global_savecol;
    bool found_star = false;
    bool found_prime = false;

    // Find starred zero in column
    for (size_t i = 0; i < global_size; i++) {
        if (global_mask_matrix[i][col] == STAR) {
            row = i;
            found_star = true;
            break;
        }
    }

    if (found_star) {
        // Find primed zero in row
        for (size_t j = 0; j < global_size; j++) {
            if (global_mask_matrix[row][j] == PRIME) {
                col = j;
                found_prime = true;
                break;
            }
        }
    }

    // Update masks
    if (found_star && found_prime) {
        global_mask_matrix[global_saverow][global_savecol] = STAR;
        global_mask_matrix[row][col] = NORMAL;
    }

    // Clear primes and reset masks
    for (size_t i = 0; i < global_size; i++) {
        for (size_t j = 0; j < global_size; j++) {
            if (global_mask_matrix[i][j] == PRIME) {
                global_mask_matrix[i][j] = NORMAL;
            }
        }
    }

    for (size_t i = 0; i < global_size; i++) {
        global_row_mask[i] = false;
        global_col_mask[i] = false;
    }

    return 2;
}

int step5(void) {
    double h = std::numeric_limits<double>::max();
    
    // Find minimum uncovered value
    for (size_t row = 0; row < global_size; row++) {
        if (!global_row_mask[row]) {
            for (size_t col = 0; col < global_size; col++) {
                if (!global_col_mask[col]) {
                    if (h > global_matrix[row][col] && global_matrix[row][col] != 0) {
                        h = global_matrix[row][col];
                    }
                }
            }
        }
    }

    // Add h to covered rows
    for (size_t row = 0; row < global_size; row++) {
        if (global_row_mask[row]) {
            for (size_t col = 0; col < global_size; col++) {
                global_matrix[row][col] += h;
            }
        }
    }

    // Subtract h from uncovered columns
    for (size_t col = 0; col < global_size; col++) {
        if (!global_col_mask[col]) {
            for (size_t row = 0; row < global_size; row++) {
                global_matrix[row][col] -= h;
            }
        }
    }

    return 3;
}

void solve_munkres(double* matrix, size_t rows, size_t columns) {
    global_size = std::max(rows, columns);
    
    // Copy input matrix to global matrix
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            global_matrix[i][j] = matrix[i * columns + j];
        }
    }

    // Initialize masks
    for (size_t i = 0; i < global_size; i++) {
        global_row_mask[i] = false;
        global_col_mask[i] = false;
    }

    // Prepare matrix values
    replace_infinites();
    minimize_along_direction(rows >= columns);
    minimize_along_direction(rows < columns);

    // Follow the steps
    int step = 1;
    while (step) {
        switch (step) {
            case 1:
                step = step1();
                break;
            case 2:
                step = step2();
                break;
            case 3:
                step = step3();
                break;
            case 4:
                step = step4();
                break;
            case 5:
                step = step5();
                break;
        }
    }

    // Store results back to input matrix
    for (size_t row = 0; row < rows; row++) {
        for (size_t col = 0; col < columns; col++) {
            if (global_mask_matrix[row][col] == STAR) {
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

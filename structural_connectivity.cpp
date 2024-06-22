
//**************************************************************************************
//**************************************************************************************

	#include <iostream>
	#include <fstream>
	#include <vector>
	#include <cmath>
	#include "random_number.cpp"

	using namespace std;

//**************************************************************************************
//**************************************************************************************

	// Define the size of each population
	const int GPe_Num = 460;
	const int STN_Num = 136;

	const float PI = 4.0 * atan(1.0);	// Definition of pi number.
	long iseed = -80L;			// Seed for random numbers.

	float P_GG = 0.05;			// GPe-GPe connection probability.
	float P_GS = 0.02;			// GPe-STN connection probability.
	float P_SG = 0.05;			// STN-GPe connection probability.

//**************************************************************************************
//**************************************************************************************

int main() {

	int counter = -1;
	float random_num;
	float norm_const = 0;

	// Write to a text file
	ofstream outfile1("structural_connectivity.txt");

	// Create an adjacency matrix
	float adj_matrix[STN_Num + GPe_Num][STN_Num + GPe_Num];

	// Populate the adjacency matrix with intra-population connectivity: STN --> STN
        for (int i = 0; i < STN_Num; i++) for (int j = 0; j < STN_Num; j++) adj_matrix[i][j] = 0;

	// Populate the adjacency matrix with intra-population connectivity: GPe --> GPe
	for (int i = STN_Num; i < STN_Num + GPe_Num; i++){
		for (int j = STN_Num; j < STN_Num + GPe_Num; j++){
			if (i == j) adj_matrix[i][j] = 0;
			else{
				random_num = ran2(&iseed);
				if (random_num <= P_GG) adj_matrix[i][j] = 1;
				else adj_matrix[i][j] = 0;
			}
		}
	}

	// Populate the adjacency matrix with inter-population connectivity: GPe --> STN
	for (int i = STN_Num; i < STN_Num + GPe_Num; i++){
		for (int j = 0; j < STN_Num; j++){
			random_num = ran2(&iseed);
			if (random_num <= P_GS) adj_matrix[i][j] = 2;
			else adj_matrix[i][j] = 0;
		}
	}


	// Populate the adjacency matrix with inter-population connectivity: STN --> GPe
	for (int i = 0; i < STN_Num; i++){
		for (int j = STN_Num; j < STN_Num + GPe_Num; j++){
			random_num = ran2(&iseed);
			if (random_num <= P_SG) adj_matrix[i][j] = 3;
			else adj_matrix[i][j] = 0;
		}
	}
	
//**************************************************************************************
//**************************************************************************************

	// Print the adjacency matrix
	for (int i = 0; i < STN_Num + GPe_Num; i++) for (int j = 0; j < STN_Num + GPe_Num; j++) outfile1 << i << '\t' << j << '\t' << adj_matrix[i][j] << endl;
	outfile1.close();

return 0;
}

//**************************************************************************************
//**************************************************************************************


//**********************************************************************************************
//**********************************************************************************************

#include <bits/stdc++.h>
#include <iostream>
#include <vector> 
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <array>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cstdio>
#include <assert.h>
#include <cmath>

#define REAL 0
#define IMAG 1

using namespace std; 


//**********************************************************************************************
//**********************************************************************************************

// Delta Function 
float Delta_Function(int x, float dt){
    
	float flag;

	if (x == 0) flag = 1.0 / dt;
	else flag = 0.0;               
	return(flag);
}

//**********************************************************************************************
//**********************************************************************************************

// Step Function 
float Step_Function(float x){
    
	float flag;

	if (x >= 0.0) flag = 1.0;
	else flag = 0.0;               
	return(flag);
}

//**********************************************************************************************
//**********************************************************************************************

// Mean (Vector)
float Mean_Vector(vector<float> &vec, int n){ 
	
	float sum = 0.0; 

	for (int i = 0; i < n; i++) sum = sum + vec[i]; 
	return sum / n; 
} 

// Mean (Array) 
float Mean_Array(float arr[], int n){ 

	float sum = 0.0; 

	for (int i = 0; i < n; i++) sum = sum + arr[i]; 
	return sum / n; 
} 


//**********************************************************************************************
//**********************************************************************************************

// Standard Deviation (Vector)
float Standard_Deviation_Vector(vector<float> &vec, int n){ 
	
	float sum = 0.0;

	for (int i = 0; i < n; i++) sum = sum + (vec[i] - Mean_Vector(vec, n)) * (vec[i] - Mean_Vector(vec, n)); 
	return sqrt(sum / (n - 1)); 
} 

// Standard Deviation (Array)
float Standard_Deviation_Array(float arr[], int n){ 

	float sum = 0.0;

	for (int i = 0; i < n; i++) sum = sum + (arr[i] - Mean_Array(arr, n)) * (arr[i] - Mean_Array(arr, n)); 
	return sqrt(sum / (n - 1)); 
} 


//**********************************************************************************************
//**********************************************************************************************

// Coefficient of Variation (Vector)
float Coefficient_Variation(vector<float> &vec, int n){ 

	return Standard_Deviation_Vector(vec, n) / Mean_Vector(vec, n); 
} 

// Coefficient of Variation (Array)
float Coefficient_Variation_Array(float arr[], int n){ 

	return Standard_Deviation_Array(arr, n) / Mean_Array(arr, n); 

} 

//**********************************************************************************************
//**********************************************************************************************

// Pearson Correlation Coefficient (Vector)
float Pearson_Correlation_Vector(vector<float> &vec1, vector<float> &vec2, int n1, int n2){ 

	float sum = 0.0;
	int n;

	if (n1 != n2){
		n = 0;
		cout << "Error: Vectors dimentions do not match!" << endl;
	}

	else n = n1;

	for (int i = 0; i < n; i++) sum = sum + ((vec1[i] - Mean_Vector(vec1, n1)) * (vec2[i] - Mean_Vector(vec2, n2))); 

	return sum / (Standard_Deviation_Vector(vec1, n1) * Standard_Deviation_Vector(vec2, n2)); 
} 

// Pearson Correlation Coefficient (Array)
float Pearson_Correlation_Array(float arr1[], float arr2[], int n1, int n2){ 

	float sum = 0.0;
	int n;

	if (n1 != n2){
		n = 0;
		cout << "Error: Arrays dimentions do not match!" << endl;
	}

	else n = n1;

	for (int i = 0; i < n; i++) sum = sum + ((arr1[i] - Mean_Array(arr1, n1)) * (arr2[i] - Mean_Array(arr2, n2))); 

	return sum / (Standard_Deviation_Array(arr1, n1) * Standard_Deviation_Array(arr2, n2)); 
} 

//******************************************************************************************
//******************************************************************************************

//**********************************************************************************************
//**********************************************************************************************


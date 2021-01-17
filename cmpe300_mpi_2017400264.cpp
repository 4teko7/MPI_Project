/*
Student Name: Bilal Tekin
Student Number: 2017400264
Compile Status: Compiling
Program Status: Working
Notes: Program is working as expected. Just compile, give parameters and run.
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <bits/stdc++.h>
using namespace std;
ifstream input; // Opens the input file
void printSingleMatrix(int matrix[], int size, int rank);
void printResult(int matrix[], int size);
vector<pair<float, float>> getMaxsAndMinsInMatrix(vector<vector<float>> matrix);
pair<int, int> findIndexOfHitAndMiss(vector<vector<float>> matrix, int index);

int main(int argc, char *argv[]) {
    int rank; // rank of the current processor
    int size; // total number of processors

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // gets the rank of the current processor
    MPI_Comm_size(MPI_COMM_WORLD, &size); // gets the total number of processors

    // ****************************************** //
    string P, N, A, M, T; //Process Number, Data Line Number, Feature Number, Iteration Number, Top Features Number
    string secondLine; // Second line of input file
    string line; //for every line of input file
    input.open(argv[1]);
    getline(input, P);  // Get Process number
    getline(input, secondLine); // Get second line
    stringstream secondLineStream(secondLine);
    secondLineStream >> N >> A >> M >> T;
    int dataLineNumber = stoi(N); //Data Line Number
    int slaveNumber = stoi(P) - 1; //Slave Number
    int amountOfLineForSlaves = dataLineNumber / slaveNumber; //Amount Of Line For Slaves
    float arr[dataLineNumber * (stoi(A) + 1) + amountOfLineForSlaves * (stoi(A) + 1)]; //All Data is put inside this array
    float pref[amountOfLineForSlaves * (stoi(A) + 1)]; //This array for data that is sent to slaves
    int weightIds[stoi(T) * slaveNumber]; //Weight Array
    if (rank == 0) {
        for (int j = 0; j < dataLineNumber * (stoi(A) + 1) + amountOfLineForSlaves * (stoi(A) + 1); j++) //Filling data array with 0's 
            arr[j] = 0;

        int i = amountOfLineForSlaves * (stoi(A) + 1); //This is for passing Master processor
        while (getline(input, line)) { //Reading all lines and putting them into data array.
            stringstream ss(line);
            do {
                if (ss.eof())
                    break;
                string variable;
                ss >> variable;
                try {
                    arr[i++] = stof(variable);
                }
                catch (exception e) {
                    break;
                }
            } while (ss);
        }
        input.close();
    }

    // Distributing all data to slave processors
    MPI_Scatter(arr, amountOfLineForSlaves * (stoi(A) + 1), MPI_FLOAT, pref, amountOfLineForSlaves * (stoi(A) + 1), MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (rank != 0) { //If the process is slave this will work.
        vector<pair<float, float>> W;
        for (int i = 0; i < stoi(A); i++) { // filling Weight array with 0 and corresponding index.
            W.push_back(make_pair(0, i));
        }

        vector<vector<float>> matrix; //Vector for data of the slave processor
        vector<float> lineMatrix; //Every line of the vectore matrix
        int count = 0;
        for (int i = 0; i < amountOfLineForSlaves; i++) {
            for (int j = 0; j < (stoi(A) + 1); j++) {
                lineMatrix.push_back(pref[count++]);
            }
            matrix.push_back(lineMatrix);
            lineMatrix.clear();
        }
        vector<pair<float, float>> maxsAndMins = getMaxsAndMinsInMatrix(matrix); //Every min and max of the array will be found via this method.

        //These loops will calculate the hits and misses and with the relief algorithm, the Weight vector will be filled.
        for (int i = 0; i < stoi(M); i++) {
            pair<int, int> hitAndMissIndexes = findIndexOfHitAndMiss(matrix, i); //Finding hit and miss of the corresponding target
            for (int j = 0; j < stoi(A); j++) {
                float hitPart = 0;
                float missPart = 0;
                if (hitAndMissIndexes.first != -1) { // If hit isn't found, hitPart will remain zero.
                
                    //With relief algorithm hitPart is calculated
                    hitPart = (fabs(matrix[i][j] - matrix[hitAndMissIndexes.first][j]) / (maxsAndMins[j].second - maxsAndMins[j].first)) / stoi(M);
                }
                if (hitAndMissIndexes.second != -1) { // If miss isn't found, missPart will remain zero

                    //With relief algorithm missPart is calculated
                    missPart = (fabs(matrix[i][j] - matrix[hitAndMissIndexes.second][j]) / (maxsAndMins[j].second - maxsAndMins[j].first)) / stoi(M);
                }
                W[j].first = W[j].first - hitPart + missPart; // Weight is updated.
            }
        }

        sort(W.rbegin(), W.rend()); //Sorting the Weight in descending order

        for (int i = 0; i < stoi(T); i++) {
            weightIds[i] = W[i].second; //Ids of highest features are put into weightIds array
        }

        sort(weightIds, weightIds + stoi(T)); //Sorting wightIds array according to index.
        printSingleMatrix(weightIds, stoi(T), rank); // printing all processors weightIds

        // Sending the weightIds that the slave processors found to master
        MPI_Send(
            /* data         = */ &weightIds,
            /* count        = */ stoi(T),
            /* datatype     = */ MPI_INT,
            /* destination  = */ 0,
            /* tag          = */ 0,
            /* communicator = */ MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD); // synchronizing processes

    if (rank == 0) { // Master process will execute this if.
        int result[slaveNumber * stoi(T)]; //Result of master will be written into this array
        int count = 0;
        for (int i = 1; i <= slaveNumber; i++) {

            // Receiving all weightIds from slave processors
            MPI_Recv(
                /* data         = */ &weightIds,
                /* count        = */ stoi(T),
                /* datatype     = */ MPI_INT,
                /* source       = */ i,
                /* tag          = */ 0,
                /* communicator = */ MPI_COMM_WORLD,
                /* status       = */ MPI_STATUS_IGNORE);

            for (int j = 0; j < stoi(T); j++) {
                result[count++] = weightIds[j]; //Filling result array with weightIds coming from slave processors
            }
        }
        sort(result, result + stoi(T) * slaveNumber); // Sorting all the ids that is in the result array
        printResult(result, stoi(T) * slaveNumber); // Printing the Master Processor result.
    }

    // ****************************************** //

    MPI_Finalize();

    return 0;
}

//Printing The Master Processor Result.
void printResult(int matrix[], int size) {
    cout << "Master P0" << " : ";
    int pre = -1;
    for (int k = 0; k < size; k++) {
        if (pre == matrix[k])
            continue;
        pre = matrix[k];
        cout << matrix[k] << " ";
    }
    cout << endl;
}

//Printing the result of all slave processors
void printSingleMatrix(int matrix[], int size, int rank) {
    cout << "Slave P" << rank << " : ";
    for (int k = 0; k < size; k++) {
        cout << matrix[k] << " ";
    }
    cout << endl;
}

// Finding index of Hit and Miss from the slave data matrix
pair<int, int> findIndexOfHitAndMiss(vector<vector<float>> matrix, int index) {
    float minOfHit = INT16_MAX;
    float minOfMiss = INT16_MAX;
    float sumOfHit = INT16_MAX;
    float sumOfMiss = INT16_MAX;
    int indexOfHit = -1;
    int indexOfMiss = -1;

    for (int i = 0; i < matrix.size(); i++) {
        float sumOfHitTmp = 0;
        float sumOfMissTmp = 0;
        for (int j = 0; j < matrix[i].size(); j++) {
            if (i != index && matrix[index][matrix[index].size() - 1] == matrix[i][matrix[i].size() - 1]) { // Checking if the line of matrix is hit
                sumOfHitTmp += fabs(matrix[index][j] - matrix[i][j]);
            }
            else if (i != index && matrix[index][matrix[index].size() - 1] != matrix[i][matrix[i].size() - 1]) { //Checking if the line of matrix is Miss
                sumOfMissTmp += fabs(matrix[index][j] - matrix[i][j]);
            }
        }
        if (sumOfHit > sumOfHitTmp && sumOfHitTmp != 0) { //Updating Hit and Hit Index
            sumOfHit = sumOfHitTmp;
            indexOfHit = i;
        }
        if (sumOfMiss > sumOfMissTmp && sumOfMissTmp != 0) { //Updating Miss and Miss Index
            sumOfMiss = sumOfMissTmp;
            indexOfMiss = i;
        }
    }

    return make_pair(indexOfHit, indexOfMiss); // Return pair of hit index and miss index
}

// Returns all maxes and mins of columns in the matrix of slave processors
vector<pair<float, float>> getMaxsAndMinsInMatrix(vector<vector<float>> matrix) {
    vector<pair<float, float>> maxsAndMins;
    float minInCol = INT16_MAX;
    float maxInCol = INT16_MIN;
    for (int i = 0; i < matrix[0].size() - 1; i++) {
        for (int j = 0; j < matrix.size(); j++) {
            if (minInCol > matrix[j][i])
                minInCol = matrix[j][i];
            if (maxInCol < matrix[j][i])
                maxInCol = matrix[j][i];
        }
        maxsAndMins.push_back(make_pair(minInCol, maxInCol));
        minInCol = INT16_MAX;
        maxInCol = INT16_MIN;
    }

    return maxsAndMins;
}
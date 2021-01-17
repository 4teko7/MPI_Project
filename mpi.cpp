#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include<bits/stdc++.h> 
using namespace std;
vector<string> matrix;
ifstream input;
void printMatrix();
void printArray(float arr[], int size, int lineSize);
void printArrayTwo(vector<vector<float>> matrix);
void printSingleMatrix(int matrix[], int size,int rank);
void printVectorOfPairs(vector<pair<float,float>> vectorOfPairs);
void printVector(vector<float> vect,int rank);
vector<pair<float,float>> getMaxsAndMinsInMatrix(vector<vector<float>> matrix);
pair<int, int> findIndexOfMiss(vector<vector<float>> matrix, int index);
int main(int argc, char *argv[])
{
    int rank; // rank of the current processor
    int size; // total number of processors

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // gets the rank of the current processor
    MPI_Comm_size(MPI_COMM_WORLD, &size); // gets the total number of processors


    // ****************************************** //
        string P,N,A,M,T;
        string secondLine;
        string line;
        // FILE *cin = fopen(argv[1], "r");
        input.open(argv[1]);
        getline(input, P); // P=number of processors
        getline(input, secondLine); // N=number of processors
        stringstream secondLineStream(secondLine);
        secondLineStream >> N >> A >> M >> T;
        // cout << P << " " << N << " " << A << " " << M << " " << T << endl;
        int dataLineNumber = stoi(N);
        int slaveNumber = stoi(P) - 1;
        int amountOfLineForSlaves = dataLineNumber / slaveNumber;
        float arr[dataLineNumber * (stoi(A) + 1) + amountOfLineForSlaves * (stoi(A) + 1)];
        float pref[amountOfLineForSlaves * (stoi(A) + 1)];
    
    if(rank == 0){
        

        for(int j = 0; j< dataLineNumber * (stoi(A) + 1) + amountOfLineForSlaves * (stoi(A) + 1);j++)
            arr[j]=0;

        int i = amountOfLineForSlaves * (stoi(A) + 1);
        while (getline(input, line))
        {

            stringstream ss(line);
            do
            {
                if (ss.eof())
                    break;
                string variable;
                ss >> variable;
                try
                {
                    // cout << stof(variable) << " ";
                    arr[i++] = stof(variable);
                }
                catch (exception e)
                {
                    cout << "hata var" << endl;
                    break;
                }
            } while (ss);
                // cout << endl;

        }
            input.close();
        
    }
            MPI_Scatter(arr,amountOfLineForSlaves * (stoi(A) + 1),MPI_FLOAT,pref,amountOfLineForSlaves * (stoi(A) + 1),MPI_FLOAT,0,MPI_COMM_WORLD);


    
    if (rank != 0) {
        vector<pair<float, float>> W;
        // float W[stoi(A)]
        for (int i = 0; i < stoi(A); i++) {
            W.push_back(make_pair(0,i));
        }

        vector<vector<float>> matrix;
        vector<float> lineMatrix;
        int count = 0;
        for(int i = 0; i<amountOfLineForSlaves; i++){
            for (int j = 0; j < (stoi(A) + 1); j++) {
                lineMatrix.push_back(pref[count++]);
            }
            matrix.push_back(lineMatrix);
            lineMatrix.clear();
            
        }
        vector<pair<float,float>> maxsAndMins = getMaxsAndMinsInMatrix(matrix);
        
        for (int i = 0; i < stoi(M); i++) {
            pair<int, int> hitAndMissIndexes = findIndexOfMiss(matrix,i);
            for (int j = 0; j < stoi(A); j++) {
                float hitPart = 0;
                float missPart = 0;
                if(hitAndMissIndexes.first != -1){
                    hitPart = (abs(matrix[i][j] - matrix[hitAndMissIndexes.first][j]) / (maxsAndMins[j].second - maxsAndMins[j].first)) / stoi(M);
                }
                if(hitAndMissIndexes.second != -1) {
                    missPart = (abs(matrix[i][j] - matrix[hitAndMissIndexes.second][j]) / (maxsAndMins[j].second - maxsAndMins[j].first)) / stoi(M);
                }
                W[j].first = W[j].first - hitPart + missPart;
            }
        }
        
        sort(W.rbegin(), W.rend());
        // vector<float> weightIds;
        int weightIds[stoi(T)];
        for (int i = 0; i < stoi(T); i++) {
            // weightIds.push_back(W[i].second);
            weightIds[i] = W[i].second;
        }

        sort(weightIds, weightIds + stoi(T));

        // sort(weightIds.begin(), weightIds.end());
        
        // cout << "************** RANK: " << rank << "*****************" << endl;
        printSingleMatrix(weightIds,stoi(T),rank);

        // printVector(weightIds,rank);
        // printVectorOfPairs(W);
        // cout<<endl;
    }


    

    // int pref[N]; // stores preferences of each // local disk on processors

    // If it's master processor, reads from input file

    // sends data from root array arr to pref array on each processor
    // MPI_Scatter(arr,N,MPI_INT,pref,N,MPI_INT,0,MPI_COMM_WORLD);



    // int masterSignal = 1;
    // while(masterSignal){

    //     if(rank!= 0){
    //         int i = 0;
    //         for(; i<N;i++){
    //             printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,pref[i] );
    //         }
    //     }

    //     if(rank==0){
    //         masterSignal=0;
    //     }

    //     MPI_Bcast(&masterSignal, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast


    // }


    // ****************************************** //


    MPI_Barrier(MPI_COMM_WORLD); // synchronizing processes
    MPI_Finalize();

    return 0;
}



void printMatrix(){
    for(int k=0; k<matrix.size(); k++){
            cout << matrix[k]<<" ";
    }
    cout<<endl;
}

void printSingleMatrix(int matrix[], int size,int rank){
    cout << "RANK: " << rank << " => ";
    for(int k=0; k<size; k++){
        cout << matrix[k] <<" ";
    }
    cout<<endl;
}
void printVector(vector<float> vect,int rank){
    cout << "RANK: " << rank << " => ";
    for(int k=0; k<vect.size(); k++){
        cout << vect[k] <<" ";
    }
    cout<<endl;
}

void printVectorOfPairs(vector<pair<float,float>> vectorOfPairs){
    for(int k=0; k<vectorOfPairs.size(); k++){
        cout << "FIRST: " << vectorOfPairs[k].first << " - SECOND: " << vectorOfPairs[k].second <<"\n";
    }
}

void printArray(float arr[], int size, int lineSize){
    cout << "*********************************************" << endl;
    for(int k=0; k < size; k++){
            cout << arr[k] << " ";
        if((k + 1) % lineSize == 0)
            cout<<endl;
        }

    cout << "*********************************************" << endl;

}


void printArrayTwo(vector<vector<float>> matrix){
    cout << "*********************************************" << endl;
    for(int i = 0; i<matrix.size();i++){
        for(int j = 0; j< matrix[i].size();j++){
            cout << matrix[i][j] << " ";
        }
        cout<<endl;
    }

    cout << "*********************************************" << endl;

}

pair<int, int> findIndexOfMiss(vector<vector<float>> matrix, int index){
    float minOfHit = INT16_MAX;
    float minOfMiss = INT16_MAX;
    float sumOfHit = INT16_MAX;
    float sumOfMiss = INT16_MAX;
    int indexOfHit = -1;
    int indexOfMiss = -1;

    for (int i = 0; i < matrix.size(); i++) {
        float sumOfHitTmp = 0;
        float sumOfMissTmp = 0;
        for (int j = 0; j < matrix.size(); j++) {
            if(i != index && matrix[index][matrix[index].size() - 1] == matrix[i][matrix[index].size() - 1]) { //HIT
                sumOfHitTmp += abs(matrix[index][j] -  matrix[i][j]);
            } else if(i != index && matrix[index][matrix[index].size() - 1] != matrix[i][matrix[index].size() - 1]){ //MISS
                sumOfMissTmp += abs(matrix[index][j] -  matrix[i][j]);
            }
        }
        if(sumOfHit > sumOfHitTmp && sumOfHitTmp != 0) {
            sumOfHit = sumOfHitTmp;
            indexOfHit = i;
        }
        if(sumOfMiss > sumOfMissTmp && sumOfMissTmp != 0){
            sumOfMiss = sumOfMissTmp;
            indexOfMiss = i;
        }
    }

    return make_pair(indexOfHit,indexOfMiss);
    
}

vector<pair<float,float>> getMaxsAndMinsInMatrix(vector<vector<float>> matrix) {
    vector<pair<float,float>> maxsAndMins;
    float minInCol = INT16_MAX;
    float maxInCol = INT16_MIN;
    for (int i = 0; i < matrix[0].size() - 1; i++) {
        for (int j = 0; j < matrix.size(); j++) {
            if(minInCol > matrix[j][i]) minInCol = matrix[j][i];
            if(maxInCol < matrix[j][i]) maxInCol = matrix[j][i];
        }
        maxsAndMins.push_back(make_pair(minInCol,maxInCol));
        minInCol = INT16_MAX;
        maxInCol = INT16_MIN;
    }

    return maxsAndMins;
}
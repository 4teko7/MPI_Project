#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
vector<string> matrix;
ifstream input;
void printMatrix();
void printArray(float arr[], int size, int lineSize);
void printArrayTwo(vector<vector<float>> matrix);
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
        // for(int i = 0; i<amountOfLineForSlaves * (stoi(A) + 1); i++){
      
        //     cout << "rank : " << rank << " - " << pref[i] << endl;
        // }

        vector<vector<float>> matrix;
        vector<float> lineMatrix;
        int count = 0;
        for(int i = 0; i<amountOfLineForSlaves; i++){
            for (int j = 0; j < (stoi(A) + 1); j++) {
                lineMatrix.push_back(pref[count++]);
            }
            matrix.push_back(lineMatrix);
            lineMatrix.clear();
            
            // printf("Process Numb %d and %d th element of my list is %f\n",rank,i+1,pref[i] );
        }

        for (int i = 0; i < stoi(M); i++) {
            pair<int, int> hitAndMissIndexes = findIndexOfMiss(matrix,i);
            // cout << hitAndMissIndexes.first << " " << hitAndMissIndexes.second << endl;
        }
        

        printArrayTwo(matrix);

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
                cout << matrix[k]<<"\n";
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
    int indexOfHit = 0;
    int indexOfMiss = 0;


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
        if(sumOfHit > sumOfHitTmp) {
            sumOfHit = sumOfHitTmp;
            indexOfHit = i;
        }
        if(sumOfMiss > sumOfMissTmp){
            sumOfMiss = sumOfMissTmp;
            indexOfMiss = i;
        }
    }

    return make_pair(indexOfHit,indexOfMiss);
    
}
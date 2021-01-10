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
        float arr[dataLineNumber * (stoi(A) + 1) + (stoi(A) + 1)];
        float pref[amountOfLineForSlaves * (stoi(A) + 1)];

    if(rank == 0){
        

        for(int j = 0; j<(stoi(A) + 1);j++)
            arr[j]=0;



        int i = amountOfLineForSlaves * (stoi(A) + 1);
        while (getline(input, line))
        {
            // vector<float> matrix;

            stringstream ss(line);
            do
            {
                if (ss.eof())
                    break;
                string variable;
                ss >> variable;
                try
                {
                    // matrix.push_back(stoi(variable));
                    cout << stoi(variable) << " ";
                    arr[i++] = stof(variable);
                }
                catch (exception e)
                {
                    break;
                }
            } while (ss);
                cout << endl;

        }
        
        // printArray(arr,dataLineNumber * (stoi(A) + 1) + (stoi(A) + 1),(stoi(A) + 1));
        // printArray(arr);
        input.close();



        // for (int i = 0; i < slaveNumber; i++){
        //     for (int j = 0; j < amountOfLineForSlaves; j++){
        //         vector<string> linesOfSlave;
        //         linesOfSlave.push_back(matrix[j]);
        //     }

    }

            MPI_Scatter(arr,amountOfLineForSlaves * (stoi(A) + 1),MPI_FLOAT,pref,amountOfLineForSlaves * (stoi(A) + 1),MPI_FLOAT,0,MPI_COMM_WORLD);

    
    if (rank != 0) {
        float W[A];
        for (int i = 0; i < A; i++)
            W[i] = 0;

        for (int i = 1; i < M; i++){
                        
        }
        
        
        int i = 0;
            for(; i<amountOfLineForSlaves * (stoi(A) + 1);i++){
                printf("Process Numb %d and %d th element of my list is %f\n",rank,i+1,pref[i] );
            }

                cout << P << " " << N << " " << A << " " << M << " " << T << endl;

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
#include <iostream>
#include <vector>
#include <algorithm> // optional (for std::shuffle if you want)
#include <ctime>
#include <random>
#include <cstdlib>
#include <fstream>


using namespace std;


void randomize(vector<int>& array);
void mergeSort(vector<int>& array, int l, int r);
void merge(vector<int>& array, int l, int m, int r);


int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <N>\n";
    return 1;
    }
    int N = stoi(argv[1]);

    vector<int> myArray(N);
    random_device rd;            
    mt19937 gen(rd());           
    uniform_int_distribution<> dist(1, 100);

  
    for (int &val: myArray){
        val = dist(gen);
    }

    cout << "Initializing the array to sort." << "\n"; 


     if (myArray.size() <= 100){
        for (int i : myArray){
            cout << i << " ";
        }
    }
    cout << "\n";
    auto start = chrono::high_resolution_clock::now();
    mergeSort(myArray, 0, myArray.size()-1);
    auto end = chrono::high_resolution_clock::now();

    cout << "The array has now been sorted. It took " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms.\n";

    ofstream outfile("results.txt", ios::app);
    outfile << N << " " 
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << "\n";
    outfile.close();

    
    for (int i : myArray) {
        cout << i << " ";
    }
    cout << "\n";

    return 0;
}

void mergeSort(vector<int>& array, int l, int r){
    
    if (l >= r){
            return;
    }
        int m = l + (r - l) / 2;
    
        mergeSort(array, l, m);
        mergeSort(array, m+1, r);

        merge(array, l, m , r);
    
}

void merge(vector<int>& array, int l, int m, int r){

    int i, j, k;

    int n1 = m - l + 1;
    int n2 = r - m;

    vector<int> leftArr(n1);
    vector<int> rightArr(n2);

    for (i = 0; i < n1; i++){
        leftArr[i] = array[l + i];
    }
    for (j= 0; j < n2; j++){
        rightArr[j] = array[m + 1 + j];
    }

    i = 0;
    j = 0;
    k = l;

    while (i < n1 && j < n2){
        if (leftArr[i] <= rightArr[j]){
            array[k] = leftArr[i];
            i++;
        } else {
            array[k] = rightArr[j];
            j++;
        }
        k++;
    }
    while ( i < n1 ){
        array[k] = leftArr[i];
        i++;
        k++;
    }
    while ( j < n2){
        array[k] = rightArr[j];
        j++;
        k++;
    }

}





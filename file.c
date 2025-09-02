#include <stdio.h>
#include <stdlib.h>
#include <time.h>



void randomize(int array[], int N);
void mergeSort(int array[], int l, int r);
void merge(int array[], int l, int m, int r);


int main() {
    srand(time((NULL)));
    int N;
    printf("How many elements are in your array?\n");
    scanf("%d", &N);
    int myArray[N];
    int length = sizeof(myArray)/sizeof(myArray[0]);

    for (int i = 0; i < N; i++){
        myArray[i] = (i + 1) * 2;
    }

    printf("Randomize the order to resort with merge sort.\n");
    randomize(myArray, N);

    for (int i = 0; i < N; i++){
        printf("%d\n", myArray[i]);
    }

    mergeSort(myArray, 0, length - 1);
    printf("Now, it's been sorted.\n");

    for (int i = 0; i < N; i++){
        printf("%d\n", myArray[i]);
    }

    return 0;
}

void randomize(int array[], int length){
    
    for (int i = length - 1; i > 0; i--){
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
    
}

void mergeSort(int array[], int l, int r){
    
    if (l < r){
        int m = l + (r - l) / 2;
    
        mergeSort(array, l, m);
        mergeSort(array, m+1, r);

        merge(array, l, m , r);
    }
  
}

void merge(int array[], int l, int m, int r){

    int i, j, k;

    int n1 = m - l + 1;
    int n2 = r - m;

    int leftArr[n1];
    int rightArr[n2];

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





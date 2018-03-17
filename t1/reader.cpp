#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#define MIN(x, y) ((x < y) ? x : y)

bool debug = true;

// ------ UTILS MALLOC / FREE -- ----
int pendingAdressesCount = 0;
void* mallocLogging(size_t size){
    void* mallocResult = malloc(size);
    if(debug) printf("Allocated \"%d\". Address: \"%p\" \n", size, mallocResult);
    pendingAdressesCount++;
    return mallocResult;
}

void freeLogging(void* var){
    free(var);
    if(debug) printf("Deallocated address: \"%p\" \n", var);
    pendingAdressesCount--;
}

// ------ MATRIX UTILS ------
typedef struct {
    long long  *array;
    long long **matrix;
    int iMax;
    int jMax;
}Image;

Image* allocateImage(int iMax, int jMax){
    long long *array;
    array = (long long*) mallocLogging(sizeof(long long)*iMax*jMax);
    
    long long **matrixPointers;
    matrixPointers = (long long**) mallocLogging(sizeof(long long*)*iMax);
    
    for(int i = 0; i < iMax; i++){
        matrixPointers[i] = &array[i*jMax] ;
    }
    Image *image = (Image*) mallocLogging(sizeof(Image));
    image->array = array;
    image->matrix = matrixPointers;
    image->iMax = iMax;
    image->jMax = jMax; 
    return image;
}

void freeImage(Image* image){
    freeLogging(image->matrix);
    freeLogging(image->array);
    freeLogging(image);
}

// ------ IMAGE UTILS ------

Image* readImage(char* filename){
    if(debug) printf("Reading %s\n", filename);
    
    FILE* file = fopen(filename, "r" );
    if (!file) {
        printf("Error: Unable to open file %s.\n\n", filename);
        exit(1);
    }
    char version[3];
    int col_len, row_len, max_gray;
    fgets(version, sizeof(version), file);
    fscanf(file, "%d", &col_len);
    fscanf(file, "%d", &row_len);
    fscanf(file, "%d", &max_gray);

    if (debug) printf("%d %d %d ",col_len, row_len, max_gray);

    Image* image = allocateImage(row_len, col_len);

    long long number;
    for(int i = 0; i < image->iMax;i++){
        for(int j = 0; j < image->jMax;j++){
            fscanf(file, "%lld", &number);
            if (number < 0){
                printf("Error: number from PGM is lower than zero. Maybe PGM file max gray scale is greater than long long?");
                exit(1);
            }
            image->matrix[i][j] = number;
        }    
    }

    fclose(file);
    return image;
}

void checkOverflow(long long greater, long long lower){
    if (lower > greater){
        printf("Error: Overflow has happened during process\n");
        exit(1);
    }
}

long long pow(long long base, int exponent){
    if (debug) printf("%lld^%d=", base, exponent);
    long long result = 1;
    for(int i = 0; i < exponent; i++){
        result = result * base;
    }

    if (debug) printf("%lld; ", result);
    return result;
}

Image* generateIntegralImage(Image* source, int powExponent){
    Image* integralImage = allocateImage(source->iMax, source->jMax);

    long long zero = 0;
    long long accumulated = 0;
    for(int i = 0; i < source->iMax;i++){
        accumulated += pow(source->matrix[i][0],powExponent);
        integralImage->matrix[i][0] = accumulated;    
        checkOverflow(integralImage->matrix[i][0], zero);
    }

    accumulated = 0;
    for(int j = 0; j < source->jMax;j++){
        accumulated += pow(source->matrix[0][j], powExponent);
        integralImage->matrix[0][j] = accumulated;
        checkOverflow(integralImage->matrix[0][j], zero);
    }

    for(int i = 1; i < source->iMax;i++){
        for(int j = 1; j < source->jMax;j++){
            long long original = pow(source->matrix[i][j], powExponent);
            long long upper = integralImage->matrix[i][j-1];
            long long left = integralImage->matrix[i-1][j];
            long long intersected = integralImage->matrix[i-1][j-1];
            long long result = original - intersected + upper + left;

            integralImage->matrix[i][j] = result;
            checkOverflow(integralImage->matrix[i][j], 0); 
        }    
    }

    return integralImage;
}

void printImage(Image* image){
    printf("Size: %d x %d\n", image->iMax, image->jMax);
    for (int i = 0; i < image->iMax; i++){
        printf("%lld: [", i);
        for (int j = 0; j < image->jMax; j++){
            printf("%lld\t", image->matrix[i][j]);
        }
        printf("]\n");
    }
}

// ------ STATISTIC UTILS ------



// ------ MAIN UTILS ------

void printEnd(){
    if(debug){
        printf("\n------------------\n");
        printf("Finished program. Pending pointers: %d \n", pendingAdressesCount);
        printf("------------------\n");
    }
}

void printStart(char * argv[]){
    if(debug){
        printf("\n------------------\n");
        printf("Starting program.\n");
        printf("The argument supplied is %s, %s\n", argv[1], argv[2]);
        printf("------------------\n");
    }
}

long readTSize(char *arg){
    long tSize = strtol(arg, NULL, 10);
    if(tSize == 0){
        printf("Error: Invalid argument 2. It should be a long");
        exit(1);
    }
    if(debug) printf("T = %ld\n", tSize);
    return tSize;
}

int main(int argc, char * argv[]){
    if( argc != 3 ) {
        printf("Call this program using 2 arguments. Filename and T-Size. Ex: './program filename.pgm 50.\n");
        exit(1);
    }

    printStart(argv);
    
    long tSize = readTSize(argv[2]);
    
    Image* source = readImage(argv[1]);
    printImage(source);
    
    Image* sumIntegralImage = generateIntegralImage(source, 1);
    printImage(sumIntegralImage);

    Image* pow2IntegralImage = generateIntegralImage(source, 2);
    printImage(pow2IntegralImage);

    

    freeImage(source);
    freeImage(sumIntegralImage);
    freeImage(pow2IntegralImage);
    
    printEnd();
}


void varianceAccessingTwice(Image *source, long tSize){
    long long lowestVariance = 9999999999999999999999; //long long highest value
    long long iLowestVariance, jLowestVariance = -1;
    long windowSize = tSize*tSize;
    for(int i = 0; i < source->iMax - (tSize -1); i++){
        for(int j = 0; j < source->jMax - (tSize -1); j++){
            // Calculate Average
            long long windowSum = 0;
            for(int iWindow = i; iWindow < tSize; iWindow++){
                for(int jWindow = j; jWindow < tSize; jWindow++){
                    windowSum += source->array[iWindow][jWindow];
                }
            }
            long long windowAvg = sumWindow / windowSize;

            // Calculate Variance
            long long sumToVar = 0;
            for(int iWindow = i; iWindow < tSize; iWindow++){
                for(int jWindow = j; jWindow < tSize; jWindow++){
                    sumToVar += pow(source->array[iWindow][jWindow] - windowAvg, 2);
                }
            }
            long long varianceResult = sumToVar / windowSize;

            if(varianceResult < lowestVariance){
                lowestVariance = varianceResult;
                iLowestVariance = i;
                jLowestVariance = j;
            }
        }
    }
}

void varianceAccessingOnce(Image *source, long tSize){
    long long lowestVariance = 9999999999999999999999; //long long highest value
    long long iLowestVariance, jLowestVariance = -1;
    long windowSize = tSize*tSize;
    for(int i = 0; i < source->iMax - (tSize -1); i++){
        for(int j = 0; j < source->jMax - (tSize -1); j++){
            // Accumulate  
            long long windowSum = 0;
            long long sumToVar = 0;
            for(int iWindow = i; iWindow < tSize; iWindow++){
                for(int jWindow = j; jWindow < tSize; jWindow++){
                    sumToVar += pow(source->array[iWindow][jWindow], 2);
                    windowSum += source->array[iWindow][jWindow];
                }
            }
            // Calculate Average
            long long windowAvg = sumWindow / windowSize;
            // Calculate Variance
            long long varianceResult = sumToVar - windowSize * pow(windowAvg, 2) 
            
            if(varianceResult < lowestVariance){
                lowestVariance = varianceResult;
                iLowestVariance = i;
                jLowestVariance = j;
            }
        }
    }
}


void varianceIntegralImage(Image *source, long tSize){
    Image* sumIntegralImage = generateIntegralImage(source, 1);
    printImage(sumIntegralImage);

    Image* pow2IntegralImage = generateIntegralImage(source, 2);
    printImage(pow2IntegralImage);

    long long lowestVariance = 9999999999999999999999; //long long highest value
    long long iLowestVariance, jLowestVariance = -1;
    long windowSize = tSize*tSize;
    for(int i = 0; i < source->iMax - (tSize -1); i++){
        for(int j = 0; j < source->jMax - (tSize -1); j++){
            long long sumToVarA = pow2IntegralImage->array[i-1][j-1]
            long long sumToVarB = pow2IntegralImage->array[i-1][(tSize -1)]
            long long sumToVarC = pow2IntegralImage->array[(tSize -1)][(tSize -1)]
            long long sumToVarD = pow2IntegralImage->array[(tSize -1)][j-1]
            long long sumToVar = sumToVarC - sumToVarB - sumToVarD + sumToVarA 
                
            long long sumToAvgA = sumIntegralImage->array[i-1][j-1]
            long long sumToAvgB = sumIntegralImage->array[i-1][(tSize -1)]
            long long sumToAvgC = sumIntegralImage->array[(tSize -1)][(tSize -1)]
            long long sumToAvgD = sumIntegralImage->array[(tSize -1)][j-1]
            long long sumToAvg = sumToAvgC - sumToAvgB - sumToAvgD + sumToAvgA
            long long windowAvg = sumToAvg / windowSize
                
            long long varianceResult = sumToVar - windowSize * pow(windowAvg, 2) 
            
            if(varianceResult < lowestVariance){
                lowestVariance = varianceResult;
                iLowestVariance = i;
                jLowestVariance = j;
            }
        }
    }
    
    freeImage(sumIntegralImage);
    freeImage(pow2IntegralImage);
}

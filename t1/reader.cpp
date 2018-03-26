#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#define MIN(x, y) ((x < y) ? x : y)

bool debugVerbose = false;
bool debugSimple = false;

// ------------------------------------------ UTILS MALLOC / FREE -- ----
int pendingAdressesCount = 0;
void* mallocLogging(size_t size){
    void* mallocResult = malloc(size);
    if(debugVerbose) printf("Allocated \"%d\". Address: \"%p\" \n", size, mallocResult);
    pendingAdressesCount++;
    return mallocResult;
}

void freeLogging(void* var){
    free(var);
    if(debugVerbose) printf("Deallocated address: \"%p\" \n", var);
    pendingAdressesCount--;
}

// ------------------------------------------ MATRIX UTILS ------------------------------------------
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

// ------------------------------------------ IMAGE UTILS ------------------------------------------
Image* readImage(char* filename){
    if(debugVerbose) printf("Reading %s\n", filename);

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

    if (debugVerbose) printf("%d %d %d ",col_len, row_len, max_gray);

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

double pow(double base, int exponent){
    double result = 1;
    for(int i = 0; i < exponent; i++){
        result = result * base;
    }

    return result;
}

Image* generateIntegralImage(Image* source, int powExponent){
    Image* integralImage = allocateImage(source->iMax, source->jMax);

    for(int i = 0; i < source->iMax;i++){
        for(int j = 0; j < source->jMax;j++){
            long long original = pow(source->matrix[i][j], powExponent);
            long long upper = j==0 ? 0 : integralImage->matrix[i][j-1];
            long long left = i==0 ? 0 : integralImage->matrix[i-1][j];
            long long intersected = i==0 || j==0 ? 0 : integralImage->matrix[i-1][j-1];
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
        printf("%d: [", i);
        for (int j = 0; j < image->jMax; j++){
            printf("%lld\t", image->matrix[i][j]);
        }
        printf("]\n");
    }
}

// ------------------------------------------ STATISTIC UTILS ------------------------------------------
typedef struct {
    int iLowestVar;
    int jLowestVar;
    double lowestVariance;
    long tSize;
    double windowAverage;
}VarianceResult;

VarianceResult* getVarianceAccessingTwice(Image *source, long tSize){
    double lowestVariance = 9999999999999; //long long highest value
    int iLowestVariance, jLowestVariance = -1;
    double windowSize = tSize*tSize;
    double windowSum, windowAvg, sumToVar, variance, windowAverage;

    for(int i = 0; i < source->iMax - (tSize -1); i++){
        for(int j = 0; j < source->jMax - (tSize -1); j++){
            // Calculate Average
            windowSum = 0;
            for(int iWindow = i; iWindow < i + tSize; iWindow++){
                for(int jWindow = j; jWindow < j + tSize; jWindow++){
                    windowSum += source->matrix[iWindow][jWindow];
                    checkOverflow(windowSum, 0); 

                }
            }
            windowAvg = windowSum / windowSize;
            // Calculate Variance
            sumToVar = 0;
            for(int iWindow = i; iWindow < i + tSize; iWindow++){
                for(int jWindow = j; jWindow < j + tSize; jWindow++){
                    sumToVar += pow(source->matrix[iWindow][jWindow] - windowAvg, 2);
                    checkOverflow(sumToVar, 0); 
                }
            }
            variance = sumToVar / windowSize;
            
            if(variance < lowestVariance){
                lowestVariance = variance;
                windowAverage = windowAvg;
                iLowestVariance = i;
                jLowestVariance = j;
            }
        }
    }

    VarianceResult *varianceResult = (VarianceResult*) mallocLogging(sizeof(VarianceResult));
    
    varianceResult->windowAverage = windowAverage;
    varianceResult->iLowestVar = iLowestVariance;
    varianceResult->jLowestVar = jLowestVariance;
    varianceResult->lowestVariance = lowestVariance; 
    varianceResult->tSize = tSize;
    
    return varianceResult;
}

VarianceResult* getVarianceAccessingOnce(Image *source, long tSize){
    double lowestVariance = 9999999999999; //long long highest value
    int iLowestVariance, jLowestVariance = -1;
    double windowSize = tSize*tSize;

    double windowSum, sumToVar, windowAvg, variance, windowAverage; 
    for(int i = 0; i < source->iMax - (tSize -1); i++){
        for(int j = 0; j < source->jMax - (tSize -1); j++){
            // Accumulate  
            windowSum = 0;
            sumToVar = 0;
            for(int iWindow = i; iWindow < i + tSize; iWindow++){
                for(int jWindow = j; jWindow < j + tSize; jWindow++){
                    sumToVar += pow(source->matrix[iWindow][jWindow], 2);
                    windowSum += source->matrix[iWindow][jWindow];
                }
            }
            // Calculate Average
            windowAvg = windowSum / windowSize;
            // Calculate Variance
            variance = (sumToVar - windowSize * pow(windowAvg, 2)) / windowSize ;
            
            if(variance < lowestVariance){
                lowestVariance = variance;
                windowAverage = windowAvg;
                iLowestVariance = i;
                jLowestVariance = j;
            }
        }
    }

    VarianceResult *varianceResult = (VarianceResult*) mallocLogging(sizeof(VarianceResult));
    
    varianceResult->windowAverage = windowAverage;
    varianceResult->iLowestVar = iLowestVariance;
    varianceResult->jLowestVar = jLowestVariance;
    varianceResult->lowestVariance = lowestVariance; 
    varianceResult->tSize = tSize;
    
    return varianceResult;
}

VarianceResult* getVarianceUsingIntegralImage(Image *source, long tSize){
    Image* sumIntegralImage = generateIntegralImage(source, 1);
    if(debugVerbose) printImage(sumIntegralImage);

    Image* pow2IntegralImage = generateIntegralImage(source, 2);
    if(debugVerbose) printImage(pow2IntegralImage);

    double lowestVariance = 9999999999999; //long long highest value
    int iLowestVariance, jLowestVariance = -1;
    double windowSize = tSize*tSize;

    double pow2ToVarA, pow2ToVarB, pow2ToVarC, pow2ToVarD, 
        pow2ToVar, sumToAvgA, sumToAvgB, sumToAvgC, sumToAvgD, 
        windowSum, windowAvg, variance, windowAverage;
    for(int i = 0; i < source->iMax - (tSize -1); i++){
        for(int j = 0; j < source->jMax - (tSize -1); j++){
            pow2ToVarA = (i == 0 || j == 0 ? 0 : pow2IntegralImage->matrix[i-1][j-1] );
            pow2ToVarB = (i == 0 ? 0 : pow2IntegralImage->matrix[i-1][(j-1 + tSize)] );
            pow2ToVarC = pow2IntegralImage->matrix[(i-1 + tSize)][(j-1 + tSize)];
            pow2ToVarD = (j == 0 ? 0 : pow2IntegralImage->matrix[(i-1 + tSize)][j-1] );
            pow2ToVar = pow2ToVarC - pow2ToVarB - pow2ToVarD + pow2ToVarA; 

            sumToAvgA = (i == 0 || j == 0 ? 0 : sumIntegralImage->matrix[i-1][j-1] );
            sumToAvgB = (i == 0 ? 0 : sumIntegralImage->matrix[i-1][(j-1 + tSize)] );
            sumToAvgC = sumIntegralImage->matrix[(i-1 + tSize)][(j-1 + tSize)];
            sumToAvgD = (j == 0 ? 0 : sumIntegralImage->matrix[(i-1 + tSize)][j-1] );
            windowSum = sumToAvgC - sumToAvgB - sumToAvgD + sumToAvgA;
            windowAvg = windowSum / windowSize;
                
            variance = (pow2ToVar - windowSize * pow(windowAvg, 2)) / windowSize ;

            if(variance < lowestVariance){
                lowestVariance = variance;
                windowAverage = windowAvg;
                iLowestVariance = i;
                jLowestVariance = j;
            }
        }
    }
    
    freeImage(sumIntegralImage);
    freeImage(pow2IntegralImage);

    VarianceResult *varianceResult = (VarianceResult*) mallocLogging(sizeof(VarianceResult));
    
    varianceResult->windowAverage = windowAverage;
    varianceResult->iLowestVar = iLowestVariance;
    varianceResult->jLowestVar = jLowestVariance;
    varianceResult->lowestVariance = lowestVariance; 
    varianceResult->tSize = tSize;
    
    return varianceResult;
}

// ------------------------------------------ MAIN UTILS ------------------------------------------
void printEnd(){
    if(debugSimple){
        printf("\n------------------------------------------\n");
        printf("Finished program. Pending pointers: %d \n", pendingAdressesCount);
        printf("------------------------------------------\n");
    }
}

void printStart(char * argv[]){
    if(debugSimple){
        printf("\n------------------------------------------\n");
        printf("Starting program.\n");
        printf("The argument supplied is %s, %s\n", argv[1], argv[2]);
        printf("------------------------------------------\n");
    }
}

long readTSize(char *arg){
    long tSize = strtol(arg, NULL, 10);
    if(tSize == 0 || tSize < -1){
        printf("Error: Invalid argument 2. It should be a positive long or -1 to automatic\n");
        exit(1);
    }
    if(debugVerbose) printf("T = %ld\n", tSize);
    return tSize;
}

typedef struct {
    VarianceResult* varianceResult;
    double cpuTimeUsed;
}ClockedVarianceResult;

void freeClockedVarianceResult(ClockedVarianceResult* result){
    freeLogging(result->varianceResult);
    freeLogging(result);
}

ClockedVarianceResult* runCalculatingTime(VarianceResult* (*f)(Image*, long), Image* source, long tSize ){
    clock_t start, end;
    double cpuTimeUsed;
    start = clock();

    VarianceResult *result = (*f)(source, tSize);

    end = clock();
    cpuTimeUsed = ((double) (end - start)) / CLOCKS_PER_SEC;

    ClockedVarianceResult* clockedResult = (ClockedVarianceResult*) mallocLogging(sizeof(ClockedVarianceResult)); 
    clockedResult->varianceResult = result;
    clockedResult->cpuTimeUsed = cpuTimeUsed;
    if(debugVerbose){
        Image *target = allocateImage(tSize,tSize);
        for(int i = 0; i < tSize; i++){ 
            for(int j = 0; j < tSize; j++){
                target->matrix[i][j] = source->matrix[result->iLowestVar + i][result->jLowestVar + j];
            }   
        }
        printImage(target);
        freeImage(target);
    }
    return clockedResult;
}

void printResult(ClockedVarianceResult* result, char const* algorithmName){
    printf("%s:\t %lf", algorithmName, result->cpuTimeUsed);
    if(debugSimple) {
        printf("\t %lf \t %d \t %d \t %f\n", 
            result->varianceResult->lowestVariance, 
            result->varianceResult->iLowestVar, 
            result->varianceResult->jLowestVar, 
            result->varianceResult->windowAverage);
    } else {
        printf(" segundos\n");
    }
}

int runAll(Image* source, long tSize){
    ClockedVarianceResult *resultTwice, *resultOnce, *resultIntegral; 
    
    resultTwice = runCalculatingTime(getVarianceAccessingTwice, source, tSize);
    resultOnce = runCalculatingTime(getVarianceAccessingOnce, source, tSize);
    resultIntegral = runCalculatingTime(getVarianceUsingIntegralImage, source, tSize);
    
    ClockedVarianceResult* result = resultTwice;
    printResult(result, "Percorrendo duas vezes");
    freeClockedVarianceResult(result);

    result = resultOnce;
    printResult(result, "Percorrendo uma vez   ");
    freeClockedVarianceResult(result);

    result = resultIntegral;
    printResult(result, "Imagens Integrais:    ");
    freeClockedVarianceResult(result);
}

/* *********************************************************
 *  Main with parameters
 * =========================================================
 * ./a.out images/desired.pgm 9
 * *********************************************************
 * ./a.out images/desired.pgm -1
 * *********************************************************/
int main(int argc, char * argv[]){
    if( argc != 3 ) {
        printf("Call this program using 2 arguments. Filename and T-Size. Ex: './program filename.pgm 50.\n");
        exit(1);
    }

    printStart(argv);
    long tSize = readTSize(argv[2]);

    Image* source = readImage(argv[1]);
    if (debugVerbose) printImage(source);
    
    ClockedVarianceResult *resultTwice, *resultOnce, *resultIntegral; 
    int runCount = 1;
    long* tSizes;
    int tSizeCount;
    if(tSize == -1){
        runCount = 5;
        tSizeCount = 8;
        tSizes = (long*) mallocLogging(sizeof(long)*tSizeCount);
        tSizes[0] = 25;
        tSizes[1] = 50;
        tSizes[2] = 75;
        tSizes[3] = 100;
        tSizes[4] = 125;
        tSizes[5] = 150;
        tSizes[6] = 175;
        tSizes[7] = 200;
    }else{
        tSizeCount = 1;
        tSizes = (long*) mallocLogging(sizeof(long)*tSizeCount);
        tSizes[0] = tSize ;
    }

    for(int run = 0; run < runCount; run++){
        for(int i = 0; i < tSizeCount; i++){
            printf("T = %ld\n", tSizes[i]);
            runAll(source, tSizes[i]);
        }
    }

    freeLogging(tSizes);
    freeImage(source);
    
    printEnd();

    return 0;
}

/* ===================================================================================
 * Relatório de Resultados
 * ===================================================================================
 * 
 * A proposta do trabalho foi avaliar as seguintes implementações de processamento de 
 * imagens:
 * 1. Implementado pela função \textit{getVarianceAccessingTwice}
		
		\item{} Implementado pela função \textit{getVarianceAccessingOnce}
		
		\item{} Implementado pela função \textit{getVarianceUsingIntegralImage}	
	\end{itemize}
	
 * 
 * 
 * 
 * 
 * 
 * ===================================================================================
 */
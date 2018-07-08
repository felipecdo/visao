#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>

#define MIN(x, y) ((x < y) ? x : y)

bool debugVerbose = false;
bool debugSimple = false;

/* ------------------------------------------ UTILS MALLOC / FREE -----------------------------------
 * Funções auxiliares para malloc e free. Elas servem para um ter um controle a mais 
 * das funções com gerenciamento de memória para garantir que não existem ponteiros 
 * pendentes
 */
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

// ------------------------------------------ MATRIX/ARRAY UTILS ------------------------------------------
typedef struct {
    long long  *array;
    long long **matrix;
    int iMax;
    int jMax;
}Image;

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

/* ------------------------------------------ IMAGE UTILS ------------------------------------------
 * Funções auxiliares para leitura da imagem
 */
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

/* 
 * Função muito simples para checagem da maioria dos overflows 
 * para exibir uma mensagem mais amigável. 
 */
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

/*
 * Função para geração de imagem integral genérica para qualquer expoente
 */
Image* generateIntegralImage(Image* source, int powExponent){
    Image* integralImage = allocateImage(source->iMax, source->jMax);

    for(int i = 0; i < source->iMax;i++){
        for(int j = 0; j < source->jMax;j++){
            long long original = pow(source->matrix[i][j], powExponent);
            // Validação de i e j para quando está nas bordas superiores e esquerda
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


// ------------------------------------------ VARIANCE UTILS ------------------------------------------
typedef struct {
    int iLowestVar;
    int jLowestVar;
    double lowestVariance;
    long tSize;
    double windowAverage;
}VarianceResult;

/*
 * Como primeira abordagem, utilizou-se a fórmula (2) da variância no enunciado, 
 * implementado pela função a seguir deste programa. Isto pois, segundo
 * a fórmula, precisamos primeiro calcular a média da janela para finalmente calcularmos a
 * variância, necessitando percorrer duas vezes o mesmo conjunto de elementos dentro da 
 * janela.
 */

VarianceResult* getVarianceAccessingTwice(Image *source, long tSize){
    double lowestVariance = 9999999999999; // um valor bem grande
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

/* 
 * Já na segunda abordagem, utilizando a fórmula (4) do enunciado, conseguimos em uma única 
 * visita à todos os itens da janela calcular a variância.
 */
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

/* 
 * E finalmente, como terceira abordagem utilizamos Imagens Integrais. Note que para que
 * o cálculo da variância ser possível foi necessário gerar duas imagens integrais: 
 * (1) sumIntegralImage que dado um ponto, continha a soma de todos os pixels à esquerda 
 * e acima dele e (2) pow2IntegralImage que contém lógica similar, porém calculando o 
 * valor da imagem original ao quadrado. Assim, foi possível com somente 8 acessos (4 em 
 * cada imagem integral) calcular a variância da janela na imagem original.
 */
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
} ClockedVarianceResult;

void freeClockedVarianceResult(ClockedVarianceResult* result){
    freeLogging(result->varianceResult);
    freeLogging(result);
}

/* 
 * Função que gera o resultado em tempo de CPU encapsulando o resultado 
 * da variância original
 */
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

/* *****************************************************************
 *  É possível executar o programa usando a seguinte configuração
 * No primeiro, como foi pedido no enunciado do EP
 * -----------------------------------------------------------------
 * ./a.out images/desired.pgm 9
 * -----------------------------------------------------------------
 * 
 * Neste segundo, podemos executar para geração das estatísticas
 * discutidas abaixo
 * -----------------------------------------------------------------
 * ./a.out images/desired.pgm -1
 * -----------------------------------------------------------------
 * *****************************************************************/
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

1. Resultados

Com o objetivo de estimar o nível de ruído em uma imagem, a proposta do trabalho foi 
avaliar algumas implementações e técnicas de processamento de imagens utilizando 
como base a variância dada uma região quadrada de largura t, a qual chamaremos de janela.

Para avaliar e comparar os resultados, cada algoritmo foi executado variando t de 25 a 200
com o passo de 25. Além disso, foi repetido o processo cinco vezes para evitar picos de 
sistema (Para reproduzir este comportamento basta chamar o programa passando t = -1, 
conforme explicado no main). Os resultados de performance abaixo contém a média de tempo 
gasto em segundos das 5 execuções.

Vamos inicialmente analisar as imagens disponibilizadas: fig01.pgm, fig02.pgm, fig03.pgm, 
fig04.pgm e fig05.pgm, imagens 512x512 da famosa "Lena", sendo que para cada imagem, foi 
inserido um ruído a mais. 
 
O resultado da menor variância independe de algoritmo e pode ser observados entre janelas na 
tabela a seguir:

-----------------------------------------------------------------------------------------
| Imagem    t=25    t=50    t=75    t=100   t=125   t=150   t=175   t=200
-----------------------------------------------------------------------------------------
| fig01     4.02    9.60    41.17   76.31   187.77  341.11  528.03  684.83
| fig02     16.90   23.27   50.69   82.15   181.73  319.87  487.81  630.58
| fig03     52.43   62.95   86.56   114.98  203.29  327.30  476.92  606.29
| fig04     111.21  127.31  148.73  175.38  253.24  363.92  496.35  612.68
| fig05     194.19  218.18  237.21  262.29  330.67  429.02  544.91  649.20
-----------------------------------------------------------------------------------------

Conseguimos perceber alguns detalhes neste resultado:
1. Para janelas com tamanho t até 150, pode-se observar que a existe uma clara escala de 
ruído onde fig01 é a que tem menos ruído (menor variância mínima) e a fig05 possui mais 
ruído (maior variância mínima). 
2. No entanto, conforme aumentamos t, a chance de termos objetos diferentes na janela 
se torna maior. Isso faz com que a menor variância mínima, na mesma figura mas para 
tamanhos de t diferentes, possam gerar variâncias maiores, observado quando comparamos 
t=200 da fig01 e da fig02, por exemplo.

Agora em termos de performance, vamos analisar fig01.pgm, a imagem com muito pouco ruído 
(os resultados são em segundos):

-----------------------------------------------------------------------------------------
| Algoritmo             t=25    t=50    t=75    t=100   t=125   t=150   t=175   t=200
-----------------------------------------------------------------------------------------
| Percorrendo 2 vezes   3.53    12.65   25.53   40.43   55.70   70.24   82.78   92.72
| Percorrendo 1 vez     2.66    9.60    19.44   30.71   42.36   53.40   63.03   70.64
| Imagens Integrais     0.02    0.02    0.02    0.02    0.02    0.02    0.02    0.02
-----------------------------------------------------------------------------------------

A tabela acima evidencia que até existe um ganho de performance de acessar somente uma vez
os itens do array ao invés de acessar duas vezes. Mas conforme aumentamos o tamanho da 
janela, mais tempo de máquina é gasto para conseguir encontrar a menor variância. Já quando
utilizamos imagens integrais, independente do tamanho da janela e ela consegue encontrar
a menor variância da imagem muito mais rápido.

Vamos observar agora a fig05.pgm, figura com mais ruído.

-----------------------------------------------------------------------------------------
| Algoritmo             t=25    t=50    t=75    t=100   t=125   t=150   t=175   t=200
-----------------------------------------------------------------------------------------
| Percorrendo 2 vezes   3.56    12.71   25.59   40.43   55.71   70.19   82.78   92.72
| Percorrendo 1 vez     2.69    9.65    19.43   30.72   42.41   53.40   63.02   70.27
| Imagens Integrais     0.02    0.02    0.02    0.02    0.02    0.02    0.02    0.02
-----------------------------------------------------------------------------------------

Fica evidente que ter ou não ruído não é um fator determinante para a execução do algoritmo
pois o tempo de execução é bem similar. No entanto, vamos observar uma figura de tamanho 
diferente: mountain.ascii.pgm (incluída no pacote de entrega na pasta images/) de 640x480.

-----------------------------------------------------------------------------------------
| Algoritmo             t=25    t=50    t=75    t=100   t=125   t=150   t=175   t=200
-----------------------------------------------------------------------------------------
| Percorrendo 2 vezes   5.21	18.73	38.07	60.67	84.47	108.18	131.28	145.99
| Percorrendo 1 vez     3.94	14.24	28.98	46.28	64.68	83.49	98.89	111.07
| Imagens Integrais     0.04	0.03	0.03	0.03	0.03	0.03	0.03	0.02
-----------------------------------------------------------------------------------------

Veja que o tempo de execução se apresenta maior na imagem em qualquer um dos algoritmos.
Interessante de se observar de que apesar do aumento, quando utilizamos imagens 
integrais, a execução ainda é bem rápida, cerca de 30 a 40 milissegundos.

2. Discussão e conclusões

Com a análise da variância foi possível perceber que ela pode ser utilizada como um 
critério de medição ruído para imagens de mesma natureza. No entanto, como pudemos 
observar, dificilmente poderemos utilizá-la para tamanho de janela distintos e imagens 
com naturezas muito diferentes devido o valor variar bastante de uma janela para outra.

Já sobre o tempo de execução, notamos que apesar do relativo ganho em percorrer a janela
somente uma vez ao invés de duas, a velocidade de execução quando utilizamos Imagens 
Integrais é muito superior. Apesar de construirmos duas imagens extras em memória com
valores maiores (long long), o resultado é extremamente rápido para o cálculo da variância.

Com o último experimento foi possível perceber também que o tamanho da imagem é um 
fator importante para qualquer um dos algoritmos apresentados vide que quando aumentamos
o número de pixels da imagem, qualquer um dos algoritmos demorou mais.

Portanto, pelos resultados apresentados fica clara a contribuição das imagens integrais 
na execução da medição de ruído através da mínima variância da imagem original. Sua 
superioridade em tempos de execução consegue processar até imagens grandes com boa 
velocidade de processamento e muito melhor do que quando percorremos os valores da janela.

3. Informações técnicas 

Todos os testes foram rodados utilizando a seguinte configuração de máquina:
Processador Intel i5-7300HQ CPU 2.5GHz x 4
12Gib de memória com 64-bits
Placa gráfica Intel® Kabylake GT2 

 * ===================================================================================
 */

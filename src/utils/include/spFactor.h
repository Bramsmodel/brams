#ifndef __spFactor__
#define __spFactor__
int ComplexRowColElimination(MatrixPtr, ElementPtr);
int FactorComplexMatrix(MatrixPtr);
int CreateInternalVectors(MatrixPtr);
int CountMarkowitz(MatrixPtr, RealVector, int);
int MarkowitzProducts(MatrixPtr, int);
ElementPtr SearchForPivot(MatrixPtr, int, int);
ElementPtr SearchForSingleton(MatrixPtr, int);
ElementPtr QuicklySearchDiagonal(MatrixPtr, int);
ElementPtr SearchDiagonal(MatrixPtr, int);
ElementPtr SearchEntireMatrix(MatrixPtr, int);
RealNumber FindLargestInCol(ElementPtr);
RealNumber FindBiggestInColExclude(MatrixPtr, ElementPtr, int);
int ExchangeRowsAndCols(MatrixPtr, ElementPtr, int);
int ExchangeColElements(MatrixPtr, int, ElementPtr, int, ElementPtr, int);
int ExchangeRowElements(MatrixPtr, int, ElementPtr, int, ElementPtr, int);
int RealRowColElimination(MatrixPtr, ElementPtr);
int ComplexRowColElimination(MatrixPtr, ElementPtr);
int UpdateMarkowitzNumbers(MatrixPtr, ElementPtr);
ElementPtr CreateFillin(MatrixPtr, int, int);
int MatrixIsSingular(MatrixPtr, int);
int ZeroPivot(MatrixPtr, int);
int WriteStatus(MatrixPtr, int);
int spcRowExchange(MatrixPtr, int, int);
int spcColExchange(MatrixPtr, int, int);
#endif

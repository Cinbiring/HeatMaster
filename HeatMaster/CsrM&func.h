#pragma once
#include <cusparse.h>

class CsrM
{
public:
    int rows, cols, nnz, nnz_temp; // wymiary m x n, nnz (number non-zero)
    int* row_offset, * d_row_offset, * colmnId, * d_colmnId;
    float* value, * d_value, * coeff_v;
    cusparseSpMatDescr_t coeff_Matrix;

    CsrM(int M, int N, float Fo, float BiF, float BiH, float Tamb, float TH);
    CsrM(int M, int N, float a, float b, float Fo, float Bi, float Tamb, float I_sk, float resistivity, float diff, float conduct, float dt);

    ~CsrM();
};

void Dense_coeff_M_podloga(int M, int N, float Fo, float BiF, float BiH, float Tamb, float TH, float** coeff_Matrix, float* coeff_vect);
void Dense_coeff_M_przewod(int M, int N, float a, float b, float Fo, float Bi, float Tamb, float I_sk, float resistivity, float diff, float conduct, float dt, float** coeff_Matrix, float* coeff_vect);
void mesh_przewod(float a, float b, int mesh_density, float* dx, int* M, int* N);
void mesh_przewod(float a, float b, int mesh_density, float* dx);
int samples(float dt,float t);

void Odczyt(int* n, float** A, int* M, int* N, float* dt);
void Odczyt_przewod(int* n, float** A, int* M, int* N, float* dt);
#pragma once
#include "CsrM&func.h" 
#include <cusparse.h>
#include <iostream>
#include <wx/wx.h>

CsrM::CsrM(int M, int N, float Fo, float BiF, float BiH, float Tamb, float TH) : nnz_temp(0), rows(N* M), cols(N* M)
{
    //zliczenie nnz i rezerwacja pamiêci row_offset colmnId i value
    nnz = (M * ((N >> 1) + 1) + (N >> 1) * (M / 3 + 2)) * 5 - 2 * N - 2 * M - N + 1;
    row_offset = new int[rows + 1];
    colmnId = new int[nnz];
    value = new float[nnz];
    coeff_v = new float[N * M];

    //wype³nienie row_offset, colmnId i value
    for (int l = 0; l >= 0 && l < M; l++)
    {
        for (int k = 0; k < N; k++)
        {
            if (k > ((N - 1) >> 1) && (l > (M / 6)) && l < (5 * M / 6)) {
                coeff_v[(long long)l * N + k] = 0;
                row_offset[l * N + k] = nnz_temp;
            }
            else if (l == M - 1 && k == 0) {   // wykrywanie k¹tów
                coeff_v[(long long)l * N + k] = 2 * BiF * Fo * Tamb;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -2 * Fo * (2 + BiF);

                colmnId[nnz_temp + 2] = l * N + k + 1;
                value[nnz_temp + 2] = 2 * Fo;

                nnz_temp += 3;
            }
            else if (l == M - 1 && k == N - 1) {
                coeff_v[(long long)l * N + k] = 2 * BiF * Fo * Tamb;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -2 * Fo * (2 + BiF);

                nnz_temp += 3;
            }
            else if (l == 5 * M / 6 && k == N - 1) {
                coeff_v[(long long)l * N + k] = 2 * BiH * Fo * TH;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k - 1;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -2 * Fo * (2 + BiH);

                colmnId[nnz_temp + 2] = (l + 1) * N + k;
                value[nnz_temp + 2] = 2 * Fo;

                nnz_temp += 3;
            }
            else if (l == 5 * M / 6 && k == (N >> 1)) { // 5i6
                coeff_v[(long long)l * N + k] = 4.0f / 3.0f * BiH * Fo * TH;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2.0f / 3.0f * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 4.0f / 3.0f * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -4 * Fo * (1 + 1.0f / 3.0f * BiH);

                colmnId[nnz_temp + 3] = l * N + k + 1;
                value[nnz_temp + 3] = 2.0f / 3.0f * Fo;

                colmnId[nnz_temp + 4] = (l + 1) * N + k;
                value[nnz_temp + 4] = 4.0f / 3.0f * Fo;

                nnz_temp += 5;
            }
            else if (l == M / 6 && k == (N >> 1)) { //6i7
                coeff_v[(long long)l * N + k] = 4.0f / 3.0f * BiH * Fo * TH;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 4.0f / 3.0f * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 4.0f / 3.0f * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -4 * Fo * (1 + 1.0f / 3.0f * BiH);

                colmnId[nnz_temp + 3] = l * N + k + 1;
                value[nnz_temp + 3] = 2.0f / 3.0f * Fo;

                colmnId[nnz_temp + 4] = (l + 1) * N + k;
                value[nnz_temp + 4] = 2.0f / 3.0f * Fo;

                nnz_temp += 5;
            }
            else if (l == M / 6 && k == N - 1) {
                coeff_v[(long long)l * N + k] = 2 * BiH * Fo * TH;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -2 * Fo * (2 + BiH);

                nnz_temp += 3;
            }
            else if (l == 0 && k == N - 1) {
                coeff_v[(long long)l * N + k] = 0;
                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k - 1;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -4 * Fo;

                colmnId[nnz_temp + 2] = (l + 1) * N + k;
                value[nnz_temp + 2] = 2 * Fo;

                nnz_temp += 3;
            }
            else if (l == 0 && k == 0) {
                coeff_v[(long long)l * N + k] = 0;
                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k;
                value[nnz_temp] = -4 * Fo;

                colmnId[nnz_temp + 1] = l * N + k + 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = (l + 1) * N + k;
                value[nnz_temp + 2] = 2 * Fo;

                nnz_temp += 3;
            }
            else if (k == 0) { // wykrywanie krawêdzi //2
                coeff_v[(long long)l * N + k] = 0;
                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -4 * Fo;

                colmnId[nnz_temp + 2] = l * N + k + 1;
                value[nnz_temp + 2] = 2 * Fo;

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else if (l == M - 1) { //3
                coeff_v[(long long)l * N + k] = 2 * BiF * Fo * Tamb;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -2 * Fo * (2 + BiF);

                colmnId[nnz_temp + 3] = l * N + k + 1;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else if (k == N - 1) { //4 lub 8
                coeff_v[(long long)l * N + k] = 0;
                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -4 * Fo;

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else if (k > (N >> 1) && l == (5 * M / 6)) { // 5
                coeff_v[(long long)l * N + k] = 2 * BiH * Fo * TH;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k - 1;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -2 * Fo * (2 + BiH);

                colmnId[nnz_temp + 2] = l * N + k + 1;
                value[nnz_temp + 2] = Fo;

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = 2 * Fo;

                nnz_temp += 4;
            }
            else if (k == (N >> 1) && (l > (M / 6)) && l < (5 * M / 6)) { // granica 6
                coeff_v[(long long)l * N + k] = 2 * BiH * Fo * TH;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -2 * Fo * (2 + BiH);

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else if (k > (N >> 1) && (l == (M / 6))) { // 7
                coeff_v[(long long)l * N + k] = 2 * BiH * Fo * TH;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -2 * Fo * (2 + BiH);

                colmnId[nnz_temp + 3] = l * N + k + 1;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else if (l == 0) { // 1
                coeff_v[(long long)l * N + k] = 0;
                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k - 1;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -4 * Fo;

                colmnId[nnz_temp + 2] = l * N + k + 1;
                value[nnz_temp + 2] = Fo;

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = 2 * Fo;

                nnz_temp += 4;
            }
            else                // niespe³nienie warunków znajdowania siê na brzegu siatki
            {
                coeff_v[(long long)l * N + k] = 0;
                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -4 * Fo;

                colmnId[nnz_temp + 3] = l * N + k + 1;
                value[nnz_temp + 3] = Fo;

                colmnId[nnz_temp + 4] = (l + 1) * N + k;
                value[nnz_temp + 4] = Fo;

                nnz_temp += 5;
            }
            //printf("\nnnz=%d       nnz_temp= %d, k= %d     l= %d\n", nnz, nnz_temp,k,l);
        }
    }
    //printf("nnz=%d       nnz_temp= %d\nPowinno byc nnz = nnz_temp \n\n", nnz, nnz_temp);

    // alokacja pamieci na GPU
    cudaError_t error = cudaMalloc((void**)&d_row_offset, (rows + 1) * sizeof(int));
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
    }
    error = cudaMalloc((void**)&d_colmnId, nnz * sizeof(float));
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
    }
    error = cudaMalloc((void**)&d_value, nnz * sizeof(float));
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
    }

    // kopiowanie danych do GPU
    cudaMemcpy(d_row_offset, row_offset, (rows + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colmnId, colmnId, nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_value, value, nnz * sizeof(float), cudaMemcpyHostToDevice);

    // tworzenie deskryptora wedlug metody CSR
    cusparseCreateCsr(&coeff_Matrix, rows, cols, nnz, d_row_offset, d_colmnId, d_value, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
}

// przeci¹¿enie konstruktora dla modelu szynoprzewodu
CsrM::CsrM(int M, int N, float a, float b, float Fo, float Bi, float Tamb, float I_sk, float resistivity, float diff, float conduct, float dt) : nnz_temp(0), rows(N* M), cols(N* M)
{
    //zliczenie nnz i rezerwacja pamiêci row_offset colmnId i value
    nnz = M * N * 5 - 2 * N - 2 * M;
    row_offset = new int[rows + 1];
    colmnId = new int[nnz];
    value = new float[nnz];
    coeff_v = new float[N * M];

    float S2 = 4 * a * a * 4 * b * b;
    float g = I_sk * I_sk * resistivity / S2;
    float g_1 = g / conduct * diff * dt;

    //wype³nienie row_offset, colmnId i value
    for (int l = 0; l >= 0 && l < M; l++)
    {
        for (int k = 0; k < N; k++)
        {
            if (l == 0 && k == N - 1) {   // 1 i 4
                coeff_v[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k - 1;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -2 * Fo * (2 + Bi);

                colmnId[nnz_temp + 2] = (l + 1) * N + k;
                value[nnz_temp + 2] = 2 * Fo;

                nnz_temp += 3;
            }
            else if (l == 0 && k == 0) { // 1 i 2
                coeff_v[(long long)l * N + k] = g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k;
                value[nnz_temp] = -4 * Fo;

                colmnId[nnz_temp + 1] = l * N + k + 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = (l + 1) * N + k;
                value[nnz_temp + 2] = 2 * Fo;

                nnz_temp += 3;
            }
            else if (l == M - 1 && k == 0) { // 2 i 3
                coeff_v[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -2 * Fo * (2 + Bi);

                colmnId[nnz_temp + 2] = l * N + k + 1;
                value[nnz_temp + 2] = 2 * Fo;

                nnz_temp += 3;
            }
            else if (l == M - 1 && k == N - 1) { // 3 i 4 
                coeff_v[(long long)l * N + k] = 4 * Bi * Fo * Tamb + g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -4 * Fo * (1 + Bi);

                nnz_temp += 3;
            }
            else if (l == 0) { // 1
                coeff_v[(long long)l * N + k] = g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = l * N + k - 1;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -4 * Fo;

                colmnId[nnz_temp + 2] = l * N + k + 1;
                value[nnz_temp + 2] = Fo;

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = 2 * Fo;

                nnz_temp += 4;
            }
            else if (k == 0) { // 2
                coeff_v[(long long)l * N + k] = g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k;
                value[nnz_temp + 1] = -4 * Fo;

                colmnId[nnz_temp + 2] = l * N + k + 1;
                value[nnz_temp + 2] = 2 * Fo;

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else if (l == M - 1) { // 3                                                       
                coeff_v[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = 2 * Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -2 * Fo * (2 + Bi);

                colmnId[nnz_temp + 3] = l * N + k + 1;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else if (k == N - 1) { // 4
                coeff_v[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = 2 * Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -2 * Fo * (2 + Bi);

                colmnId[nnz_temp + 3] = (l + 1) * N + k;
                value[nnz_temp + 3] = Fo;

                nnz_temp += 4;
            }
            else {
                coeff_v[(long long)l * N + k] = g_1;

                row_offset[l * N + k] = nnz_temp;

                colmnId[nnz_temp] = (l - 1) * N + k;
                value[nnz_temp] = Fo;

                colmnId[nnz_temp + 1] = l * N + k - 1;
                value[nnz_temp + 1] = Fo;

                colmnId[nnz_temp + 2] = l * N + k;
                value[nnz_temp + 2] = -4 * Fo;

                colmnId[nnz_temp + 3] = l * N + k + 1;
                value[nnz_temp + 3] = Fo;

                colmnId[nnz_temp + 4] = (l + 1) * N + k;
                value[nnz_temp + 4] = Fo;

                nnz_temp += 5;
            }

            //printf("\nnnz=%d       nnz_temp= %d, k= %d     l= %d\n", nnz, nnz_temp,k,l);
        }
    }
    //printf("nnz=%d       nnz_temp= %d\nPowinno byc nnz = nnz_temp \n\n", nnz, nnz_temp);

    //sprawdzenie CsrM
    /*
    printf("\n\nrow_off:    ");
    for (int i = 0; i < rows; i++)
        printf(" %d ", row_offset[i]);
    printf("\n\ncolmn_Id:    ");
    for (int i = 0; i < nnz; i++)
        printf(" %d ", colmnId[i]);
    printf("\n\nvalue:    ");
    for (int i = 0; i < nnz; i++)
        printf(" %g ", value[i]);*/

    // alokacja pamieci na GPU
    cudaError_t error = cudaMalloc((void**)&d_row_offset, (rows + 1) * sizeof(int));
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
    }
    error = cudaMalloc((void**)&d_colmnId, nnz * sizeof(float));
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
    }
    error = cudaMalloc((void**)&d_value, nnz * sizeof(float));
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
    }

    // kopiowanie danych do GPU
    cudaMemcpy(d_row_offset, row_offset, (rows + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colmnId, colmnId, nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_value, value, nnz * sizeof(float), cudaMemcpyHostToDevice);

    // tworzenie deskryptora wedlug metody CSR
    cusparseCreateCsr(&coeff_Matrix, rows, cols, nnz, d_row_offset, d_colmnId, d_value, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F);
}

CsrM::~CsrM()
{
    cusparseDestroySpMat(coeff_Matrix);
    delete[] row_offset;
    delete[] colmnId;
    delete[] value;
    cudaFree(d_row_offset);
    cudaFree(d_colmnId);
    cudaFree(d_value);
}




void Dense_coeff_M_podloga(int M, int N, float Fo, float BiF, float BiH, float Tamb, float TH, float** coeff_Matrix, float* coeff_vect)
{

    // tworzenie macierzy wspó³czynników
    for (int l = M - 1; l >= 0 && l <= M - 1; l--)
    {
        for (int k = 0; k < N; k++)
        {
            if (k > ((N - 1) >> 1) && (l > (M / 6)) && l < (5 * M / 6))
                break;
            else if (l == M - 1 && k == 0) {   // wykrywanie k¹tów
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiF);
                coeff_vect[(long long)l * N + k] = 2 * BiF * Fo * Tamb;
            }
            else if (l == M - 1 && k == N - 1) {
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiF);
                coeff_vect[(long long)l * N + k] = 2 * BiF * Fo * Tamb;
            }
            else if (l == 5 * M / 6 && k == N - 1) {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiH);
                coeff_vect[(long long)l * N + k] = 2 * BiH * Fo * TH;
            }
            else if (l == 5 * M / 6 && k == (N >> 1)) { // 5i6
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2.0f / 3.0f * Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 4.0f / 3.0f * Fo;
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 4.0f / 3.0f * Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2.0f / 3.0f * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo * (1 + 1.0f / 3.0f * BiH);
                coeff_vect[(long long)l * N + k] = 4.0f / 3.0f * BiH * Fo * TH;
            }
            else if (l == M / 6 && k == (N >> 1)) { //6i7
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2.0f / 3.0f * Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 4.0f / 3.0f * Fo;
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2.0f / 3.0f * Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 4.0f / 3.0f * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo * (1 + 1.0f / 3.0f * BiH);
                coeff_vect[(long long)l * N + k] = 4.0f / 3.0f * BiH * Fo * TH;
            }
            else if (l == M / 6 && k == N - 1) {
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiH);
                coeff_vect[(long long)l * N + k] = 2 * BiH * Fo * TH;
            }
            else if (l == 0 && k == N - 1) {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
            }
            else if (l == 0 && k == 0) {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
            }
            else if (k == 0) { // wykrywanie krawêdzi
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
            }
            else if (l == M - 1) {
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiF);
                coeff_vect[(long long)l * N + k] = 2 * BiF * Fo * Tamb;
            }
            else if (k == N - 1) {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
            }
            else if (k > (N >> 1) && l == (5 * M / 6)) {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiH);
                coeff_vect[(long long)l * N + k] = 2 * BiH * Fo * TH;
            }
            else if (k == (N >> 1) && (l > (M / 6)) && l < (5 * M / 6)) { // granica 6
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiH);
                coeff_vect[(long long)l * N + k] = 2 * BiH * Fo * TH;
            }
            else if (k > (N >> 1) && (l == (M / 6))) {
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + BiH);
                coeff_vect[(long long)l * N + k] = 2 * BiH * Fo * TH;
            }
            else if (l == 0) {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
            }
            else                // niespe³nienie warunków znajdowania siê na brzegu siatki
            {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
            }
        }
    }
}

void Dense_coeff_M_przewod(int M, int N, float a, float b, float Fo, float Bi, float Tamb, float I_sk, float resistivity, float diff, float conduct, float dt, float** coeff_Matrix, float* coeff_vect)
{
    float S2 = 4 * a * a * 4 * b * b;
    float g = I_sk * I_sk * resistivity / S2;
    float g_1 = g / conduct * diff * dt;
    // tworzenie macierzy wspó³czynników
    for (int l = M - 1; l >= 0 && l <= M - 1; l--)
    {
        for (int k = 0; k < N; k++)
        {
            if (l == 0 && k == N - 1) {   // 1 i 4
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + Bi);
                coeff_vect[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;
            }
            else if (l == 0 && k == 0) { // 1 i 2
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
                coeff_vect[(long long)l * N + k] = g_1;
            }
            else if (l == M - 1 && k == 0) { // 2 i 3
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + Bi);
                coeff_vect[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;
            }
            else if (l == M - 1 && k == N - 1) { // 3 i 4 
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo * (1 + Bi);
                coeff_vect[(long long)l * N + k] = 4 * Bi * Fo * Tamb + g_1;
            }
            else if (l == 0) { // 1
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
                coeff_vect[(long long)l * N + k] = g_1;
            }
            else if (k == 0) { // 2
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
                coeff_vect[(long long)l * N + k] = g_1;
            }
            else if (l == M - 1) { // 3                                                       
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + Bi);
                coeff_vect[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;
            }
            else if (k == N - 1) { // 4
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = 2 * Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -2 * Fo * (2 + Bi);
                coeff_vect[(long long)l * N + k] = 2 * Bi * Fo * Tamb + g_1;
            }
            else {
                coeff_Matrix[(long long)(l + 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)(l - 1) * N + k][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k + 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k - 1][(long long)l * N + k] = Fo;
                coeff_Matrix[(long long)l * N + k][(long long)l * N + k] = -4 * Fo;
                coeff_vect[(long long)l * N + k] = g_1;
            }
        }
    }
}

void mesh_przewod(float a, float b, int mesh_density, float* dx, int* M, int* N)
{
    float temp_a = a;
    float temp_b = b;

    long long scale = 1;

    long long int_a = temp_a;
    long long int_b = temp_b;

    while (int_a < temp_a || int_b < temp_b)
    {
        scale = scale * 10;

        temp_a = a * scale;
        temp_b = b * scale;

        int_a = a * scale;
        int_b = b * scale;
    }

    //NWD algorytm Euklidesa
    int temp;
    long long nwd_a = int_a;
    long long nwd_b = int_b;
    while (nwd_b != 0)
    {
        temp = nwd_a % nwd_b;
        nwd_a = nwd_b;
        nwd_b = temp;
    }

    *M = int_b / nwd_a * mesh_density + 1;
    *N = int_a / nwd_a * mesh_density + 1;
    *dx = (float)nwd_a / scale / mesh_density;
}
void mesh_przewod(float a, float b, int mesh_density, float* dx)
{
    float temp_a = a;
    float temp_b = b;

    long long scale = 1;

    long long int_a = temp_a;
    long long int_b = temp_b;

    while (int_a < temp_a || int_b < temp_b)
    {
        scale = scale * 10;

        temp_a = a * scale;
        temp_b = b * scale;

        int_a = a * scale;
        int_b = b * scale;
    }

    //NWD algorytm Euklidesa
    int temp;
    long long nwd_a = int_a;
    long long nwd_b = int_b;
    while (nwd_b != 0)
    {
        temp = nwd_a % nwd_b;
        nwd_a = nwd_b;
        nwd_b = temp;
    }

    *dx = (float)nwd_a / scale / mesh_density;
}

int samples(float dt, float t)
{
    long long int_dt = dt;
    long long int_t = t;

    float temp_dt = dt;
    float temp_t = t;

    long long scale = 1;

    while (int_dt < temp_dt || int_t < temp_t)
    {
        scale = scale * 10;

        temp_dt = dt * scale;
        temp_t = t * scale;

        int_dt = dt * scale;
        int_t = t * scale;
    }
    //wartoœci dt i t zosta³y zmienione na ca³kowite, aby unikn¹æ dodatkowyc zaokr¹gleñ
    //zaokr¹glenia float mog¹ jednak wp³ywaæ podczas zamiany z float na int 
    int samples = int_t / int_dt;
    samples++; //do uwzglêdnienia próbki zerowej
    samples = int_t % int_dt == 0 ? samples:++samples;
    
    return samples;
}


void Odczyt(int* n, float** A, int* M, int* N, float* dt)
{
    printf("probka: %d,  czas: %f [s]\n", *n, (*dt) * (*n));
    for (int l = (*M) - 1; l >= 0; l--)
    {
        for (int k = 0; k < (*N); k++)
        {
            if (k > (*N >> 1) && (l > (*M / 6)) && l < (5 * (*M) / 6))
                break;
            printf("%5.2f ", A[*n][l * (*N) + k]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void Odczyt_przewod(int* n, float** A, int* M, int* N, float* dt)
{
    printf("probka: %d,  czas: %f [s]\n", *n, (*dt) * (*n));
    for (int l = (*M) - 1; l >= 0; l--)
    {
        for (int k = 0; k < (*N); k++)
        {
            printf("%5.2f ", A[*n][l * (*N) + k]);
        }
        printf("\n");
    }
    printf("\n\n");
}
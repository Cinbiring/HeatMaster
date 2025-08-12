// HeatMaster.cpp
// IDE: Visual Studio 2022
#pragma once
#include <iostream>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include "CsrM&func.h"
#include <wx/wx.h>
//#pragma warning(disable : 6001)
#pragma warning(disable : 4996)
#define MAX_WATCHDOG_COUNT 20

union {
    char byte;
    struct {
        unsigned char repeat : 1;
        unsigned char dense_matrix : 1;
        unsigned char busbar_model : 1;
        unsigned char : 1;
        unsigned char : 1;
        unsigned char : 1;
        unsigned char : 1;
        unsigned char : 1;
    };
} Flags;



int calc(int busbar_model, bool dense_matrix, int mesh_density, float t, float L, float dt, float a_H, float a_F, float Tamb, float TH, float density_mat,
    float specific_heat_mat, float conduct_mat, float compute_accurancy, float a, float b, float I_sk, float resistivity, int* time_samples, float*** node, int N, int M, float dx)
{
    int k, l, n = 0;
    float* difference_Matrix;
    float Fo, BiF, BiH;
    float diff_mat;
    float** coeff_Matrix, * coeff_vect;
    const float alpha_plus = 1, alpha_minus = -1, betha = 1;
    int Watchdog;
    CsrM* CsrM1 = nullptr;
    

    //debug \/
    /*
    L = 0.5f; mesh_density = 6; specific_heat_mat = 800; conduct_mat = 1; density_mat = 19; a_H = 120; a_F = 9; dt = 0.1f;
    time_samples = 71; compute_accurancy = 0.01f; Tamb = 14; TH = 27; dx = L / (4 * mesh_density); t = dt * (time_samples - 1);
    Flags.dense_matrix = 0; Flags.busbar_model = 0;

    a = 0.04; b = 0.01; mesh_density = 5; specific_heat_mat = 910; conduct_mat = 229; density_mat = 2720; a_F = 9; dt = 0.1f;
    time_samples = 3600; compute_accurancy = 110.01f; Tamb = 20; I_sk = 200000; resistivity = 0.000000028264; t = dt * (time_samples - 1);
    Flags.dense_matrix = 0; Flags.busbar_model = 1;

    switch (Flags.busbar_model)
    {
    case 0:
        //dla podłogi
        M = abs(6 * mesh_density + 1);    // 0b0111 1111 1111 1111 = 32767
        N = abs((mesh_density << 2) + 1);
        printf("rozmiary M x N: %d x %d\n", M, N);
    case 1:
        mesh_przewod(a, b, mesh_density, &dx, &M, &N);
    }
    */
    //      /\
    
    //poprawne ustawienie flagi
    if (busbar_model == 1)
        Flags.busbar_model = 1;
    else
        Flags.busbar_model = 0;
    
    
    //wybór sposobu obliczeń macierz wpółczynników gęsta/rzadka
    if (dense_matrix == true)
    {
        Flags.dense_matrix = 1;
        printf("\nrealizacja macierzy gestej\n");
    }
    else
    {
        Flags.dense_matrix = 0;
        printf("\nrealizacja macierzy rzadkiej\n");
    }

    //wxLogError(wxString::Format("busbar: %d, gęsta: %d", Flags.busbar_model, Flags.dense_matrix));

    int NM = N * M; // w celu uniknięcia tego samego mnożenia wielokrotnie
    if (M <= 0 || N <= 0)
    {
        printf("\nBlad tworzenia siatki!\n\n");
        return -10;
    }


    //rezerwacja przestrzeni dla wektora wspolczynnikow
    // może lepiej byłoby użyć cudaMallocHost()
    coeff_vect = (float*)calloc((long long)NM + 1, sizeof(float));
    if (coeff_vect == NULL) {
        printf("Blad alokacji pamieci dla coeff_vect\n");
        return 1;
    }

    diff_mat = (float)conduct_mat / (specific_heat_mat * density_mat);
    // liczby podobienstwa
    Fo = (diff_mat * dt) / (dx * dx);
    BiF = a_F * dx / conduct_mat;
    BiH = a_H * dx / conduct_mat;

    // alokowanie pamięci dla węzłów siatki
    (*node) = (float**)calloc(*time_samples, sizeof(float*)); // n

    for (int i = 0; i < *time_samples; i++)
    {
        (*node)[i] = (float*)calloc(NM, sizeof(float)); // l i k
    }
    printf("\n");

    // alokacja pamieci dla macierzy roznic
    cudaMallocHost(&difference_Matrix, NM * sizeof(float));

    // Inicjalizacja i wypisanie wartości w chwili n = 0
    switch(Flags.busbar_model)
    { 
        case 0:
            for (l = M - 1; l >= 0 && l <= M - 1; l--)
            {
                for (k = 0; k < N; k++)
                {
                    if (k >(N >> 1) && (l > (M / 6)) && l < (5 * (M) / 6))
                        (*node)[0][l * N + k] = TH;
                    else
                    (*node)[0][l * N + k] = Tamb;
                }
            }
            break;
        case 1:
            for (l = M - 1; l >= 0 && l <= M - 1; l--)
            {
                for (k = 0; k < N; k++)
                {
                    (*node)[0][l * N + k] = Tamb;
                }
            }
            break;
    }

    switch (Flags.busbar_model)
    {
    case 0:
        Odczyt(&n, (*node), &M, &N, &dt); break;
    case 1:
        Odczyt_przewod(&n, (*node), &M, &N, &dt); break;
    }

    float* d_present_value, * d_previous_value, * d_prev_i_value, * d_coeff_Matrix, * d_coeff_vect, * d_difference_Matrix;
    int diff_index; // jako indeks największej róznicy
    cublasHandle_t handle_cublas;
    cublasCreate(&handle_cublas);
    cublasStatus_t error_dn;
    cusparseDnVecDescr_t X, Y;
    cusparseHandle_t handle_sparse;
    cusparseCreate(&handle_sparse);


    //alokacja pamięci na GPU
    cudaError_t error = cudaMalloc(&d_present_value, NM * sizeof(float));            // X
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
        cudaFree(d_present_value);
        return 1;
    }

    error = cudaMalloc(&d_previous_value, NM * sizeof(float));                      // Y
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
        cudaFree(d_previous_value);
        return 1;
    }

    error = cudaMalloc(&d_prev_i_value, NM * sizeof(float));                      // Xprev
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
        cudaFree(d_prev_i_value);
        return 1;
    }

    error = cudaMalloc(&d_coeff_vect, NM * sizeof(float));                          // coeff_v
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
        cudaFree(d_coeff_vect);
        return 1;
    }

    error = cudaMalloc(&d_difference_Matrix, NM * sizeof(float));                   // diff_M
    if (error != cudaSuccess) {
        printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
        cudaFree(d_difference_Matrix);
        return 1;
    }

    //tworzenie macierzy współczynników i rozwiazanie obliczen podlug spospbu tworzenia macierzy wspolczynnikow
    switch (Flags.dense_matrix)
    {
    case 1: // DLA DENSE COEFFICIENT MATRIX
        //rezerwacja przestrzeni dla macierzy wspolczynnikow
        error = cudaMalloc(&d_coeff_Matrix, NM * NM * sizeof(float));                   // coeff_M
        if (error != cudaSuccess) {
            printf("Blad alokacji pamieci na GPU: %s", cudaGetErrorString(error));
            cudaFree(d_coeff_Matrix);
            return 1;
        }
        
        coeff_Matrix = (float**)calloc((long long)NM, sizeof(float*));
        if (coeff_Matrix == NULL) {
            printf("Blad alokacji pamieci dla coeff_Matrix\n");
            return 1;
        }
        for (int i = 0; i < NM; i++)
        {
            coeff_Matrix[i] = (float*)calloc((long long)NM, sizeof(float));
            if (coeff_Matrix[i] == NULL) {
                printf("Blad alokacji pamieci dla coeff_Matrix[%d]\n", i);
                return 1;
            }
        }

        switch (Flags.busbar_model)
        {
        case 0:
            Dense_coeff_M_podloga(M, N, Fo, BiF, BiH, Tamb, TH, coeff_Matrix, coeff_vect);
            break;
        case 1:
            Dense_coeff_M_przewod(M, N, a, b, Fo, BiF, Tamb, I_sk, resistivity, diff_mat, conduct_mat, dt, coeff_Matrix, coeff_vect);
            break;
        }


        //przeniesienie danych do GPU
        cudaMemcpy(d_present_value, (*node)[0], NM * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_previous_value, d_present_value, NM * sizeof(float), cudaMemcpyDeviceToDevice);
        for (int i = 0;i < NM;i++)
            cudaMemcpy(d_coeff_Matrix + i * NM, coeff_Matrix[i], NM * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_coeff_vect, coeff_vect, NM * sizeof(float), cudaMemcpyHostToDevice);

        for (n = 1;n < *time_samples;n++)
        {
            Watchdog = 0;
            Flags.repeat = true;
            while (Flags.repeat)
            {
                cudaMemcpy(d_prev_i_value, d_present_value, NM * sizeof(float), cudaMemcpyDeviceToDevice); // Xprev <- X
                cublasSgemv_v2(handle_cublas, CUBLAS_OP_N, NM, NM, &alpha_plus, d_coeff_Matrix, NM, d_present_value, 1, &betha, d_previous_value, 1); // Y <- nowy wynik
                cublasSaxpy_v2(handle_cublas, NM, &alpha_plus, d_coeff_vect, 1, d_previous_value, 1); // Y=Y+coeff_vect
                cudaMemcpy(d_present_value, d_previous_value, NM * sizeof(float), cudaMemcpyDeviceToDevice); // X <- Y
                cudaMemcpy(d_difference_Matrix, d_previous_value, NM * sizeof(float), cudaMemcpyDeviceToDevice); //diff_M <- Y
                cudaMemcpy(d_previous_value, (*node)[n - 1], NM * sizeof(float), cudaMemcpyHostToDevice); // Y <- wynik z wczesniejszej chwili
                cublasSaxpy_v2(handle_cublas, NM, &alpha_minus, d_prev_i_value, 1, d_difference_Matrix, 1); // diff_M = diff_M - Xprev
                cublasIsamax_v2(handle_cublas, NM, d_difference_Matrix, 1, &diff_index);
                cudaMemcpy(difference_Matrix, d_difference_Matrix, NM * sizeof(float), cudaMemcpyDeviceToHost); // diff_M -> do CPU

                Flags.repeat = (abs(difference_Matrix[diff_index - 1]) < compute_accurancy) ? false : true;
                Watchdog++;
                Flags.repeat = Watchdog < MAX_WATCHDOG_COUNT ? Flags.repeat : false;
            }
            cudaMemcpy(d_previous_value, d_present_value, NM * sizeof(float), cudaMemcpyDeviceToDevice); // Y <- wynik z aktualnej chwili [X]
            cudaMemcpy((*node)[n], d_present_value, NM * sizeof(float), cudaMemcpyDeviceToHost); // Xcpu <- X
            //wxLogError(wxString::Format("node: %f", (*node)[n][0]));
            /*
            switch (Flags.busbar_model)
            {
            case 0:
                Odczyt(&n, (*node), &M, &N, &dt); break;
            case 1:
                Odczyt_przewod(&n, (*node), &M, &N, &dt);  break;
            }
            */
            //printf("prob: %d        ostatnia roznica: %f\n\n", probka, difference_Matrix[diff_index - 1]);
        }
        break;
    case 0: // DLA SPARSE COEFFICIENT MATRIX


        switch (Flags.busbar_model)
        {
        case 0:
            CsrM1 = new CsrM(M, N, Fo, BiF, BiH, Tamb, TH); break;
        case 1:
            CsrM1 = new CsrM(M, N, a, b, Fo, BiF, Tamb, I_sk, resistivity, diff_mat, conduct_mat, dt); break;
        }

        cusparseCreateDnVec(&X, NM, d_present_value, CUDA_R_32F);
        cusparseCreateDnVec(&Y, NM, d_previous_value, CUDA_R_32F);
        cusparseStatus_t error_sp;
        void* d_buffer = NULL;
        size_t bufferSize = 0;
        //przeniesienie danych do GPU
        cudaMemcpy(d_present_value, (*node)[0], NM * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_previous_value, d_present_value, NM * sizeof(float), cudaMemcpyDeviceToDevice);
        cudaMemcpy(d_coeff_vect, CsrM1->coeff_v, NM * sizeof(float), cudaMemcpyHostToDevice);

        cusparseSpMV_bufferSize(handle_sparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha_plus, CsrM1->coeff_Matrix, X, &betha, Y, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, &bufferSize);
        error = cudaMalloc(&d_buffer, bufferSize);
        if (error != cudaSuccess) {
            printf("Blad alokacji pamieci na GPU (bufor): %s", cudaGetErrorString(error));
            cudaFree(d_buffer);
            return 1;
        }

        for (n = 1;n < *time_samples;n++)
        {
            Watchdog = 0;
            Flags.repeat = true;
            while (Flags.repeat)
            {
                error = cudaMemcpy(d_prev_i_value, d_present_value, NM * sizeof(float), cudaMemcpyDeviceToDevice); // Xprev <- X
                if (error != cudaSuccess) {
                    printf("Blad kopiowania danych w GPU [Xprev <- X]: %s", cudaGetErrorString(error)); return 1;
                }

                // Y <- nowy wynik
                error_sp = cusparseSpMV(handle_sparse, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha_plus, CsrM1->coeff_Matrix, X, &betha, Y, CUDA_R_32F, CUSPARSE_SPMV_ALG_DEFAULT, d_buffer);
                if (error_sp != CUSPARSE_STATUS_SUCCESS) {
                    printf("Blad obliczen MV: %d", error_sp); return 1;
                }

                error = cudaDeviceSynchronize();
                if (error != cudaSuccess) {
                    printf("Blad synchronizacji urzadzenia: %s\n", cudaGetErrorString(error));
                    return 1;
                }

                // Y=Y+coeff_vect
                error_dn = cublasSaxpy_v2(handle_cublas, NM, &alpha_plus, d_coeff_vect, 1, d_previous_value, 1);
                if (error_dn != CUBLAS_STATUS_SUCCESS) {
                    printf("Blad dodawania wektorow: %d\n", error_dn); return 1;
                }

                // X <- Y
                error = cudaMemcpy(d_present_value, d_previous_value, NM * sizeof(float), cudaMemcpyDeviceToDevice);
                if (error != cudaSuccess) {
                    printf("Blad kopiowania danych w GPU [X <- Y]: %s", cudaGetErrorString(error)); return 1;
                }

                //diff_M <- Y
                error = cudaMemcpy(d_difference_Matrix, d_previous_value, NM * sizeof(float), cudaMemcpyDeviceToDevice);
                if (error != cudaSuccess) {
                    printf("Blad kopiowania danych w GPU [Macierz dokladnosi <- Y]: %s\n", cudaGetErrorString(error)); return 1;
                }

                // Y <- wynik z wczesniejszej chwili
                error = cudaMemcpy(d_previous_value, (*node)[n - 1], NM * sizeof(float), cudaMemcpyHostToDevice);
                if (error != cudaSuccess) {
                    printf("Blad kopiowania danych do GPU [Y <- node[n-1] ]: %s\n", cudaGetErrorString(error)); return 1;
                }

                // diff_M = diff_M - Xprev
                cublasSaxpy_v2(handle_cublas, NM, &alpha_minus, d_prev_i_value, 1, d_difference_Matrix, 1);
                cublasIsamax_v2(handle_cublas, NM, d_difference_Matrix, 1, &diff_index);

                // diff_M -> do CPU
                error = cudaMemcpy(difference_Matrix, d_difference_Matrix, NM * sizeof(float), cudaMemcpyDeviceToHost);
                if (error != cudaSuccess) {
                    printf("Blad kopiowania danych do Hosta [Macierz dokladnosci]: %s\n", cudaGetErrorString(error)); return 1;
                }
                Flags.repeat = abs(difference_Matrix[diff_index - 1]) < compute_accurancy ? false : true;
                Watchdog++;
                Flags.repeat = Watchdog < MAX_WATCHDOG_COUNT ? Flags.repeat : false;
            }
            error = cudaMemcpy(d_previous_value, d_present_value, NM * sizeof(float), cudaMemcpyDeviceToDevice); // Y <- wynik z aktualnej chwili [X]
            if (error != cudaSuccess) {
                printf("Blad kopiowania danych w GPU [Y <- X]: %s", cudaGetErrorString(error)); return 1;
            }
            error = cudaMemcpy((*node)[n], d_present_value, NM * sizeof(float), cudaMemcpyDeviceToHost); // Xcpu <- X
            if (error != cudaSuccess) {
                printf("Blad kopiowania danych do Hosta: %s", cudaGetErrorString(error)); return 1;
            }
            /* Odczyt danych tutaj niepotrzebny
            switch (Flags.busbar_model)
            {
            case 0:
                Odczyt(&n, (*node), &M, &N, &dt); break;
            case 1:
                Odczyt_przewod(&n, (*node), &M, &N, &dt);  break;
            }
            */
            //printf("prob: %d        ostatnia roznica: %f\n\n", proba, difference_Matrix[diff_index - 1]);
        }
        cudaFree(d_buffer);

        // podgląd macierzy rzadkiej
        /*
        printf("\n\nrow_offset:");
        for (int i = 0;i < CsrM1.rows;i++)
            printf(" %d ",CsrM1.row_offset[i]);
        printf("\n\ncolumnID:");
        for (int i = 0;i < CsrM1.nnz;i++)
            printf(" %d ", CsrM1.colmnId[i]);
        printf("\n\nValues:");
        for (int i = 0;i < CsrM1.nnz;i++)
            printf(" %5.3f ", CsrM1.value[i]);
        printf("\n\n");
        */
        break;
    }

    // Dealokowanie pamięci      
    cublasDestroy(handle_cublas);
    cusparseDestroy(handle_sparse);
    ///*
    if (!Flags.dense_matrix)
    {
        cusparseDestroyDnVec(X);
        cusparseDestroyDnVec(Y);
    }
    //*/

    cudaFree(d_difference_Matrix);
    cudaFree(d_coeff_vect);
    cudaFree(d_prev_i_value);
    cudaFree(d_previous_value);
    cudaFree(d_present_value);
    cudaFreeHost(difference_Matrix);

    if (Flags.dense_matrix)
    {
        for (int i = 0; i < NM; i++) {              // coeff_Matrix
            free(coeff_Matrix[i]);
        }
        free(coeff_Matrix);
        cudaFree(d_coeff_Matrix);
    }
    free(coeff_vect);                           // coeff_vect

    delete CsrM1;
    CsrM1 = nullptr;

    return 0;
}
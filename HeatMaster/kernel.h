#pragma once

int calc(int busbar_model, bool dense_matrix, int mesh_density, float t, float L, float dt, float a_H, float a_F, float Tamb, float TH, float density_mat,
    float specific_heat_mat, float conduct_mat, float compute_accurancy, float a, float b, float I_sk, float resistivity, int* time_samples, float*** node, int N, int M, float dx);
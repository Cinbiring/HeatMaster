#pragma once
#include <wx/wx.h>
#include "IsoLines.h"

class MainFrame : public wxFrame {
public:
    MainFrame(const wxString& title);
    ~MainFrame();

    wxChoice* Get_Model();
    float Get_L(), Get_a(), Get_b();
    int Get_mesh_density();
    float Get_specific_heat(), Get_conduct_mat(), Get_density_mat();
    float Get_a_H(), Get_a_F();
    float Get_dt(), Get_sim_time();
    float Get_compute_accurancy();
    float Get_T_amb(), Get_T_H();
    float Get_I_sk(), Get_resistivity();

    IsoLines* IsoLines1;
    wxSlider* time_Slider;
    int time_samples;
    int N, M;
    float** node;
private:
    wxChoice* Model;
    wxCheckBox* Dense_M;
    wxStaticText* current_time;
    
    float L, a, b, dx;
    int mesh_density;
    float specific_heat, conduct_mat, density_mat;
    float a_H, a_F;
    float dt, sim_time;
    float compute_accurancy;
    float T_amb, T_H;
    float I_sk, resistivity;

    void OnOpen(wxCommandEvent& event);
    void OnInput(wxCommandEvent& event);
    void OnExit(wxCommandEvent& event);
    void OnKeyPressed(wxKeyEvent& event);
    void OnChoice(wxCommandEvent& event);
    void OnRunButton(wxCommandEvent& event);
    void OnTimeSlider(wxCommandEvent& event);
};

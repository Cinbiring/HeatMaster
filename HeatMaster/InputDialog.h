#pragma once
#include <wx/wx.h>
#include <wx/spinctrl.h>

class InputDialog : public wxDialog {
public:
    InputDialog(wxFrame* parent, const wxString& title);

    float GetValue_L() const;
    float GetValue_a() const;
    float GetValue_b() const;
    int GetValue_mesh_density() const;
    float GetValue_specific_heat() const;
    float GetValue_conduct_mat() const;
    float GetValue_density_mat() const;
    float GetValue_a_H() const;
    float GetValue_a_F() const;
    float GetValue_dt() const;
    float GetValue_sim_time() const;
    float GetValue_compute_accurancy() const;
    float GetValue_T_amb() const;
    float GetValue_T_H() const;
    float GetValue_I_sk() const;
    float GetValue_resistivity() const;
private:
    int Model();
    wxString L(), a(), b(), mesh_density(), specific_heat(), conduct_mat(), density_mat(), a_H(), a_F(), dt(), sim_time(), compute_accurancy(), T_amb(), T_H(), I_sk(), resistivity();
    float dx();
    wxFrame* parent_ptr;
    wxTextCtrl* L_field;
    wxTextCtrl* a_field;
    wxTextCtrl* b_field;
    wxSpinCtrl* mesh_density_field;
    wxStaticText* dx_label;
    wxTextCtrl* specific_heat_field;
    wxTextCtrl* conduct_mat_field;
    wxTextCtrl* density_mat_field;
    wxTextCtrl* a_H_field;
    wxTextCtrl* a_F_field;
    wxTextCtrl* dt_field;
    wxTextCtrl* time_field;
    wxSpinCtrl* time_samples_field;
    wxTextCtrl* compute_accurancy_field;
    wxTextCtrl* T_amb_field;
    wxTextCtrl* T_H_field;
    wxTextCtrl* I_sk_field;
    wxTextCtrl* resistivity_field;

    void Update_dx_Label(wxCommandEvent& event);
};
#include "InputDialog.h"
#include <wx\wx.h>
#include <wx/spinctrl.h>
#include <wx/valnum.h>
#include "MainFrame.h"
#include "CsrM&func.h"

InputDialog::InputDialog(wxFrame* parent, const wxString& title)
    : wxDialog(parent, wxID_ANY, title, wxDefaultPosition, wxSize(400, 550)) {
    wxBoxSizer* mainSizer = new wxBoxSizer(wxVERTICAL);
    

    // Walidacja wprowadzonych danych
    wxFloatingPointValidator<double> valid_input_double(20, nullptr, wxNUM_VAL_ZERO_AS_BLANK | wxNUM_VAL_NO_TRAILING_ZEROES); // bez liter i 20 cyfr po przeciku
    valid_input_double.SetRange(0.00000000000000001, std::numeric_limits<double>::max()); // Zakres: od 0 do maksimum
    
    wxBoxSizer* line1Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line2Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line3Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line4Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line5Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line6Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line7Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line8Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line9Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line10Sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer* line11Sizer = new wxBoxSizer(wxHORIZONTAL);
    mainSizer->Add(line1Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line2Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line3Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line4Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line5Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line6Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line7Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line8Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line9Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line10Sizer, 0, wxEXPAND | wxALL, 5);
    mainSizer->Add(line11Sizer, 0, wxEXPAND | wxALL, 5);

    switch (Model())
    {
    default:
        break;
    case 0:                                 //Dane wejœciowe dla modelu pod³ogi
        // d³ugoœæ 'L' dla modelu pod³ogi
        
        line1Sizer->Add(new wxStaticText(this, wxID_ANY, "D³ugoœæ L [m]:"), 0, wxALL, 5);
        L_field = new wxTextCtrl(this, wxID_ANY, L(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
        line1Sizer->Add(L_field, 1, wxALL, 5);
        
        //wspó³czynnik przenikalnoœci cieplnej
        line7Sizer->Add(new wxStaticText(this, wxID_ANY, "Przenikalnoœæ cieplna od kana³ów [W/(m^2*K)]:"), 0, wxALL, 5);
        a_H_field = new wxTextCtrl(this, wxID_ANY, a_H(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
        line7Sizer->Add(a_H_field, 1, wxALL, 5);
        
        //Temperatura w kana³ach
        line11Sizer->Add(new wxStaticText(this, wxID_ANY, "temperatura w kana³ach [^C]:"), 0, wxALL, 5);
        T_H_field = new wxTextCtrl(this, wxID_ANY, T_H(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
        line11Sizer->Add(T_H_field, 1, wxALL, 5);
        break;
    case 1:                                 //Dane wejœciowe dla modelu przewodu p³askiego
        //d³ugoœci 'a' i 'b'
        line1Sizer->Add(new wxStaticText(this, wxID_ANY, "D³ugoœæ a [m]:"), 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
        a_field = new wxTextCtrl(this, wxID_ANY, a(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
        line1Sizer->Add(a_field, 1, wxALL, 5);
        a_field->Bind(wxEVT_TEXT, &InputDialog::Update_dx_Label, this);

        line1Sizer->Add(new wxStaticText(this, wxID_ANY, "b [m]:"), 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
        b_field = new wxTextCtrl(this, wxID_ANY, b(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
        line1Sizer->Add(b_field, 1, wxALL, 5);
        b_field->Bind(wxEVT_TEXT, &InputDialog::Update_dx_Label, this);

        //elektryczne parametry
        line11Sizer->Add(new wxStaticText(this, wxID_ANY, "I_sk [A]:"), 0, wxALL, 5);
        I_sk_field = new wxTextCtrl(this, wxID_ANY, I_sk(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
        line11Sizer->Add(I_sk_field, 1, wxALL, 5);
 
        line11Sizer->Add(new wxStaticText(this, wxID_ANY, "rezystywnoœæ [Ohm*m]:"), 0, wxALL, 5);
        resistivity_field = new wxTextCtrl(this, wxID_ANY, resistivity(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
        line11Sizer->Add(resistivity_field, 1, wxALL, 5);

        break;
    }
    
    //gêstoœæ siatki
    line2Sizer->Add(new wxStaticText(this, wxID_ANY, "Gêstoœæ siatki:"), 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
    mesh_density_field = new wxSpinCtrl(this, wxID_ANY, mesh_density(), wxDefaultPosition, wxSize(50, -1));
    mesh_density_field->SetRange(1,100);
    line2Sizer->Add(mesh_density_field, 1, wxALL, 5);
    mesh_density_field->Bind(wxEVT_SPINCTRL, &InputDialog::Update_dx_Label, this);
    
    dx_label = new wxStaticText(this, wxID_ANY, wxString::Format("dx [m]: %g", dx()));
    line2Sizer->Add(dx_label, 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
    line2Sizer->AddSpacer(60);

    //w³aœciwoœci materia³u
    line3Sizer->Add(new wxStaticText(this, wxID_ANY, "Ciep³o w³aœciwe materia³u [J/(kg*K)]:"), 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
    specific_heat_field = new wxTextCtrl(this, wxID_ANY, specific_heat(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line3Sizer->Add(specific_heat_field, 1, wxALL, 5);

    line4Sizer->Add(new wxStaticText(this, wxID_ANY, "przewodnoœæ cieplna materia³u  [W/(m*K)]:"), 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
    conduct_mat_field = new wxTextCtrl(this, wxID_ANY, conduct_mat(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line4Sizer->Add(conduct_mat_field, 1, wxALL, 5);

    line5Sizer->Add(new wxStaticText(this, wxID_ANY, "Gêstoœæ materia³u [kg/m^3]:"), 0, wxALL | wxALIGN_CENTER_VERTICAL, 5);
    density_mat_field = new wxTextCtrl(this, wxID_ANY, density_mat(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line5Sizer->Add(density_mat_field, 1, wxALL, 5);

    //wspó³czynnik przenikalnoœci cieplnej
    line6Sizer->Add(new wxStaticText(this, wxID_ANY, "Przenikalnoœæ cieplna na powierzchni [W/(m^2*K)]:"), 0, wxALL, 5);
    a_F_field = new wxTextCtrl(this, wxID_ANY, a_F(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line6Sizer->Add(a_F_field, 1, wxALL, 5);

    //czas
    line8Sizer->Add(new wxStaticText(this, wxID_ANY, "okres próbkowania [s]:"), 0, wxALL, 5);
    dt_field = new wxTextCtrl(this, wxID_ANY, dt(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line8Sizer->Add(dt_field, 1, wxALL, 5);

    line8Sizer->Add(new wxStaticText(this, wxID_ANY, "czas [s]:"), 0, wxALL, 5);
    time_field = new wxTextCtrl(this, wxID_ANY, sim_time(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line8Sizer->Add(time_field, 1, wxALL, 5);

    //dok³adnoœæ obliczeñ
    line9Sizer->Add(new wxStaticText(this, wxID_ANY, "dok³adnoœæ obliczeñ [^C]:"), 0, wxALL, 5);
    compute_accurancy_field = new wxTextCtrl(this, wxID_ANY, compute_accurancy(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line9Sizer->Add(compute_accurancy_field, 1, wxALL, 5);

    //Temperatura otoczenia
    line10Sizer->Add(new wxStaticText(this, wxID_ANY, "temperatura otoczenia [^C]:"), 0, wxALL, 5);
    T_amb_field = new wxTextCtrl(this, wxID_ANY, T_amb(), wxDefaultPosition, wxSize(50, -1), 0, valid_input_double);
    line10Sizer->Add(T_amb_field, 1, wxALL, 5);

    // Przyciski OK/Anuluj
    mainSizer->AddStretchSpacer();
    mainSizer->Add(CreateButtonSizer(wxOK | wxCANCEL), 0, wxALL | wxALIGN_CENTER, 5);

    SetSizer(mainSizer);
}

//funkcje zwracaj¹ce dane wprowadzone do okienek
float InputDialog::GetValue_L() const {
    double value;
    if (L_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value); //rzutowanie double na flaot
    }
    return 0.0f; // Domyœlna wartoœæ, jeœli nie uda siê konwersja
}

float InputDialog::GetValue_a() const {
    double value;
    if (a_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_b() const {
    double value;
    if (b_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}


int InputDialog::GetValue_mesh_density() const {
    if (mesh_density_field) {
        return mesh_density_field->GetValue();
    }
    return 0.0f;
}

float InputDialog::GetValue_specific_heat() const {
    double value;
    if (specific_heat_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_conduct_mat() const {
    double value;
    if (conduct_mat_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_density_mat() const {
    double value;
    if (density_mat_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_a_H() const {
    double value;
    if (a_H_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_a_F() const {
    double value;
    if (a_F_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_dt() const {
    double value;
    if (dt_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_sim_time() const {
    double value;
    if (time_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_compute_accurancy() const {
    double value;
    if (compute_accurancy_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_T_amb() const {
    double value;
    if (T_amb_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_T_H() const {
    double value;
    if (T_H_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_I_sk() const {
    double value;
    if (I_sk_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}
float InputDialog::GetValue_resistivity() const {
    double value;
    if (resistivity_field->GetValue().ToDouble(&value)) {
        return static_cast<float>(value);
    }
    return 0.0f;
}

// uzyskanie informacji o wyborze modelu z MainFrame
int InputDialog::Model() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    wxChoice* model = parent_ptr->Get_Model();
    if (model) {
        int selectedValue = model->GetSelection();
        return selectedValue;
    }
}

//wpisanie wartoœci w okienka
wxString InputDialog::L() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float L = parent_ptr->Get_L();
    wxString text_L = wxString::Format("%g", L);
    return text_L;
}

wxString InputDialog::a() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float a = parent_ptr->Get_a();
    wxString text_a = wxString::Format("%g", a);
    return text_a;
}

wxString InputDialog::b() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float b = parent_ptr->Get_b();
    wxString text_b = wxString::Format("%g", b);
    return text_b;
}

wxString InputDialog::mesh_density() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    int mesh_density = parent_ptr->Get_mesh_density();
    wxString text_mesh_density = wxString::Format("%d", mesh_density);
    return text_mesh_density;
}

wxString InputDialog::specific_heat() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float specific_heat = parent_ptr->Get_specific_heat();
    wxString text_specific_heat = wxString::Format("%g", specific_heat);
    return text_specific_heat;
}

wxString InputDialog::conduct_mat() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_conduct_mat();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::density_mat() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_density_mat();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::a_H() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_a_H();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::a_F() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_a_F();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::dt() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_dt();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::sim_time() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_sim_time();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::compute_accurancy() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_compute_accurancy();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::T_amb() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_T_amb();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::T_H() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_T_H();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::I_sk() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_I_sk();
    wxString text = wxString::Format("%g", value);
    return text;
}

wxString InputDialog::resistivity() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float value = parent_ptr->Get_resistivity();
    wxString text = wxString::Format("%g", value);
    return text;
}

// do inicjalizacji wartoœci dx
float InputDialog::dx() {
    MainFrame* parent_ptr = dynamic_cast<MainFrame*>(GetParent());
    float dx;
    if (parent_ptr->Get_mesh_density()) {
        switch (Model())
        {
        case 0:
            //dla pod³ogi
            dx = GetValue_L() / (4 * parent_ptr->Get_mesh_density());
            break;
        case 1:
            // dla szynoprzewodu
            mesh_przewod(GetValue_a(), GetValue_b(), parent_ptr->Get_mesh_density(), &dx);
            break;
        }
        return dx;
    }
    wxMessageBox("B³¹d odczytu danych");
}
// aktualizacja wartoœci dx
void InputDialog::Update_dx_Label(wxCommandEvent& event)
{
    float dx = mesh_density_field->GetValue();
    switch (Model())
    {
    case 0:
        //dla pod³ogi
        dx = GetValue_L() / (4 * dx);
        break;
    case 1:
        // dla szynoprzewodu
        mesh_przewod(GetValue_a(), GetValue_b(), mesh_density_field->GetValue(), &dx);
        break;
    }
    dx_label->SetLabel(wxString::Format("dx [m]: %g", dx));
}
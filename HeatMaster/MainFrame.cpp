#include <wx/wx.h>
#include <wx/slider.h>
#include "MainFrame.h"
#include "InputDialog.h"
#include "kernel.h"
#include <wx/glcanvas.h>
#include "IsoLines.h"
#include "CsrM&func.h"

enum IDs {
	Input_ID = 2
};

MainFrame::MainFrame(const wxString& title): wxFrame(nullptr, wxID_ANY, title), 
L(0.2), a(0.04), b(0.01), mesh_density(5), specific_heat(910), conduct_mat(229), density_mat(2720), a_H(18), a_F(9), dt(0.01), sim_time(36),
compute_accurancy(1), T_amb(14), T_H(27), I_sk(20000), resistivity(0.000000028264), time_samples(0), N(1), M(1), dx(1)
{	
	//inicjalizacja tablicy wartoœci, aby unikn¹æ b³êdów i komplikacji algorytmów
	node = (float**)calloc(1, sizeof(float*));

	//górny pasek
	wxMenuBar* menuBar = new wxMenuBar();
	SetMenuBar(menuBar);
	//menu "Plik"
	wxMenu* menuFile = new wxMenu();
	menuFile->Append(wxID_OPEN, "&Otwórz plik\tCtrl-O", "Otwórz plik");
	menuBar->Append(menuFile, "&Plik");
	menuFile->AppendSeparator();
	menuFile->Append(wxID_EXIT, "Wyjœcie", "Zamknij program");
	//menu "Parametry"
	wxMenu* menuInput = new wxMenu();
	menuInput->Append(Input_ID, "&Wprowad¿ dane\tCtrl-P", "WprowadŸ parametry uk³adu");
	menuBar->Append(menuInput, "&Parametry");
	
	//obs³uga zdarzeñ górnego paska
	Bind(wxEVT_MENU, &MainFrame::OnOpen, this, wxID_OPEN);
	Bind(wxEVT_MENU, &MainFrame::OnExit, this, wxID_EXIT);
	Bind(wxEVT_MENU, &MainFrame::OnInput, this, Input_ID);

	wxPanel* panel = new wxPanel(this);

	wxBoxSizer* mainSizer = new wxBoxSizer(wxVERTICAL);
	wxBoxSizer* widthSizer1 = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer* widthSizer2 = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer* widthSizer3 = new wxBoxSizer(wxHORIZONTAL);
	wxBoxSizer* widthSizer4 = new wxBoxSizer(wxHORIZONTAL);
	mainSizer->Add(widthSizer1, 0, wxEXPAND | wxALL, 0);
	mainSizer->Add(widthSizer2, 0, wxEXPAND | wxALL, 0);
	mainSizer->Add(widthSizer3, 0, wxEXPAND | wxALL, 0);
	mainSizer->Add(widthSizer4, 1, wxEXPAND | wxALL, 0);

	// wybór modelu
	wxArrayString models;
	models.Add("Ogrzewana pod³oga");
	models.Add("P³aski przewód DC");

	//okno wyboru modelu
	Model = new wxChoice(panel, wxID_ANY, wxDefaultPosition, wxSize(200, -1), models);
	Model->Select(1);
	widthSizer1->Add(Model, 0, wxALL, 10);
	Model->Bind(wxEVT_CHOICE, &MainFrame::OnChoice, this);

	
	//wybór sposobu obliczeñ ( macierz gêsta/rzadka )
	Dense_M = new wxCheckBox(panel, wxID_ANY, "wymuœ macierz gêst¹", wxDefaultPosition);
	widthSizer1->Add(Dense_M, 0, wxALL, 14);

	//suwak czasu
	time_Slider = new wxSlider(panel, wxID_ANY,0,0,1,wxDefaultPosition,wxSize(400,-1), wxSL_VALUE_LABEL);
	widthSizer2->Add(time_Slider, 1, wxALL, 5);
	time_Slider->Bind(wxEVT_SLIDER, &MainFrame::OnTimeSlider, this);
	//time_Slider->Show(false);

	current_time = new wxStaticText(panel, wxID_ANY, "t = ?", wxDefaultPosition, wxDefaultSize);
	widthSizer3->AddStretchSpacer();
	widthSizer3->Add(current_time, 0, wxALL, 5);
	widthSizer3->AddStretchSpacer();
	//current_time->Show(false);

	//Izolinie
	IsoLines1 = new IsoLines(panel);
	widthSizer4->Add(IsoLines1, 1, wxEXPAND | wxALL, 5);
	IsoLines1->SetMinSize(wxSize(200,200));
	//IsoLines1->Show(false);

	//przycisk rozpoczynaj¹cy obliczenia
	wxButton* button_RUN = new wxButton(panel, wxID_ANY, "Uruchom", wxDefaultPosition, wxSize(100, 50));
	button_RUN->Bind(wxEVT_BUTTON, &MainFrame::OnRunButton, this);
	widthSizer1->AddStretchSpacer();
	widthSizer1->Add(button_RUN, 0, wxALL, 10);

	//najmnieszy rozmiar okna
	panel->SetSizer(mainSizer);
	mainSizer->SetSizeHints(this);

	//wykrycie wciœniêcia klawiszy
	panel->Bind(wxEVT_CHAR_HOOK, &MainFrame::OnKeyPressed, this);


}

// destruktor

MainFrame::~MainFrame()
{
	
	for (int i = 0; i < time_samples; i++) {
		free(node[i]);
	}
	free(node);
}

// inne metody
void MainFrame::OnExit(wxCommandEvent& event) {
	Close(true);
}

void MainFrame::OnOpen(wxCommandEvent& event) {
	wxMessageBox("Opcja jeszcze nie dostêpna", "Brak implementacji", wxOK | wxICON_INFORMATION);
}

void MainFrame::OnInput(wxCommandEvent& event) {
	InputDialog InputDialog(this, "WprowadŸ wartoœci");
	if (InputDialog.ShowModal() == wxID_OK)
	{
		int choice = Model->GetSelection();

		wxString msg;
		
		switch (choice)
		{
		case 0: 
			
			L = InputDialog.GetValue_L();
			a_H = InputDialog.GetValue_a_H();
			T_H = InputDialog.GetValue_T_H();
			msg.Printf("Wartoœæ 1: %.2f\n", L);
			break;
		case 1:
			a = InputDialog.GetValue_a();
			b = InputDialog.GetValue_b();
			I_sk = InputDialog.GetValue_I_sk();
			resistivity = InputDialog.GetValue_resistivity();
			break;
		}
		mesh_density = InputDialog.GetValue_mesh_density();
		specific_heat = InputDialog.GetValue_specific_heat();
		conduct_mat = InputDialog.GetValue_conduct_mat();
		density_mat = InputDialog.GetValue_density_mat();
		a_F = InputDialog.GetValue_a_F();
		dt = InputDialog.GetValue_dt();
		sim_time = InputDialog.GetValue_sim_time();
		compute_accurancy = InputDialog.GetValue_compute_accurancy();
		T_amb = InputDialog.GetValue_T_amb();
	}
}

void MainFrame::OnKeyPressed(wxKeyEvent& event) {
	if (event.ControlDown() && event.GetKeyCode() == 'P') { // klawisz ctrl+'P'
		OnInput(wxCommandEvent());
	}
	else {
		event.Skip();
	}
}

void MainFrame::OnChoice(wxCommandEvent& event) {
	switch (Model->GetSelection())
	{
	case 0:
		L = 0.15; 
		specific_heat = 840;
		conduct_mat = 1; 
		density_mat = 2000; 
		a_H = 60; 
		a_F = 9;
		break;
	case 1:
		a = 0.04f;
		b = 0.01;
		specific_heat = 910;
		conduct_mat = 229;
		density_mat = 2720;
		a_F = 9;
		resistivity = 0.000000028264;
		break;
	}
}

void MainFrame::OnRunButton(wxCommandEvent& event)
{
	// gêstoœæ siatki
	switch (Model->GetSelection())
	{
	case 0:
		//dla pod³ogi
		M = abs(6 * mesh_density + 1);
		N = abs((mesh_density << 2) + 1);
		dx = L / (4 * mesh_density);
		//printf("rozmiary M x N: %d x %d\n", M, N);
		break; 
	case 1:
		//dla szynoprzewodu
		mesh_przewod(a, b, mesh_density, &dx, &M, &N);
		break;
	}

	//----------------------------------------------------------
	//Zwolnienie pamiêci przed przed zaalokowaniem nowej
	//----------------------------------------------------------

	for (int i = 0; i < time_samples; i++) {    //node
		free(node[i]);
	}
	free(node);

	time_samples = samples(dt, sim_time);
	
	// alokowanie pamiêci dla wêz³ów siatki
	node = (float**)calloc(time_samples, sizeof(float*)); // n

	for (int i = 0; i < time_samples; i++)
	{
		node[i] = (float*)calloc(N * M, sizeof(float)); // l i k
	}
	
	calc(Model->GetSelection(), Dense_M->GetValue(), mesh_density, sim_time, L, dt, a_H,  a_F, T_amb, T_H, density_mat, specific_heat, conduct_mat, compute_accurancy, a, b, I_sk, resistivity, &time_samples, &node, N, M, dx);
	time_Slider->SetRange(0, time_samples-1);
	time_Slider->SetValue(time_samples-1);
	IsoLines1->Refresh();
}

void MainFrame::OnTimeSlider(wxCommandEvent& event)
{
	current_time->SetLabel(wxString::Format("t = %g [s]", dt * time_Slider->GetValue()));
	IsoLines1->Refresh();
}

//transfer danych
wxChoice* MainFrame::Get_Model() 
{
	return Model;
}
 
float MainFrame::Get_L()
{
	return L;
}

float MainFrame::Get_a()
{
	return a;
}
float MainFrame::Get_b()
{
	return b;
}
int MainFrame::Get_mesh_density()
{
	return mesh_density;
}
float MainFrame::Get_specific_heat()
{
	return specific_heat;
}
float MainFrame::Get_conduct_mat()
{
	return conduct_mat;
}
float MainFrame::Get_density_mat()
{
	return density_mat;
}
float MainFrame::Get_a_H()
{
	return a_H;
}
float MainFrame::Get_a_F()
{
	return a_F;
}
float MainFrame::Get_dt()
{
	return dt;
}
float MainFrame::Get_sim_time()
{
	return sim_time;
}
float MainFrame::Get_compute_accurancy()
{
	return compute_accurancy;
}
float MainFrame::Get_T_amb()
{
	return T_amb;
}
float MainFrame::Get_T_H()
{
	return T_H;
}
float MainFrame::Get_I_sk()
{
	return I_sk;
}
float MainFrame::Get_resistivity()
{
	return resistivity;
}
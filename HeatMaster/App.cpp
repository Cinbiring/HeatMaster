#include <wx/wx.h>
#include "MainFrame.h"
#include "App.h"

wxIMPLEMENT_APP(App);

bool App::OnInit() {
	MainFrame* mainFrame = new MainFrame("electroHeat calc");
	mainFrame->SetClientSize(1000, 800);
	mainFrame->Center();
	mainFrame->Show();
	return true;
}
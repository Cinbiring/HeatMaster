#pragma once
#include <wx/wx.h>
#include <wx/glcanvas.h>

class MainFrame;

class IsoLines : public wxGLCanvas
{
public:
	IsoLines(wxWindow* parent);
	~IsoLines();
private:
	int N, M;
	
	MainFrame* parentFrame;
	wxGLContext* context;

	void OnPaint(wxPaintEvent&);

	void DrawIsoLines();
};


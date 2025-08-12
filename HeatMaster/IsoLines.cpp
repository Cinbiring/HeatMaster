#include "IsoLines.h"
#include <wx/wx.h>
#include<wx/glcanvas.h>
#include "MainFrame.h"

IsoLines::IsoLines(wxWindow* parent) : wxGLCanvas(parent, wxID_ANY, nullptr)
{
	// utworzenie kontekstu OpenGL
	context = new wxGLContext(this);

	Bind(wxEVT_PAINT, &IsoLines::OnPaint, this);
}

IsoLines::~IsoLines()
{
	delete context;
}

void IsoLines::OnPaint(wxPaintEvent&)
{
	//Ustawienie kontekstu OpenGL
	SetCurrent(*context);

	//Pobranie rozmiaru izolinii
	int width, height;
	GetClientSize(&width, &height);

	//Ustawienie propocji widoku
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//Czyszczenie ekarnu
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

	//Rysowanie linii izometrycznych
	DrawIsoLines();

	//Wyœwietlenie zawaratoœci
	SwapBuffers();
}

void IsoLines::DrawIsoLines()
{
	MainFrame* parentFrame = dynamic_cast<MainFrame*>(GetParent()->GetParent());
	// Pobranie wartoœci z klasy nadrzêdnej
	float** node = parentFrame->node;
	N = parentFrame->N;
	M = parentFrame->M;
	int sample = parentFrame->time_Slider->GetValue();

	//liczba izolinii
	const int nbiso = 20;

	// Zakres temperatur i kolorów
	const float minTemp = 10.0f;
	const float maxTemp = 100.0f;
	const int numColors = 256;
	
	// Gradient kolorów niebieski-czerwony
	float colorMap[numColors][3]; // Tablica kolorów
	for (int i = 0; i < numColors; i++)
	{
		float t = (float)i / (numColors - 1); // skalowanie do jednoœci
		colorMap[i][0] = t;//(1 - t); // gradientowanie czerwonego
		colorMap[i][1] = 0; // wy³¹czenie zielonego
		colorMap[i][2] = 1-t; // gradientowanie niebieskiego
	}

	// Rysowanie
	glBegin(GL_QUADS);
	for (int l = 0; l < M - 1; l++)
	{
		for (int k = 0; k < N - 1; k++)
		{
			// wartoœci temperatury dla danego wêz³a
			float t1 = node[sample][l * N + k];
			float t2 = node[sample][l * N + k + 1];
			float t3 = node[sample][(l + 1) * N + k + 1];
			float t4 = node[sample][(l + 1) * N + k];

			// normalizacja temperatur do zakresu <0,1>
			float norm_t1 = (t1 - minTemp) / (maxTemp - minTemp);
			float norm_t2 = (t2 - minTemp) / (maxTemp - minTemp);
			float norm_t3 = (t3 - minTemp) / (maxTemp - minTemp);
			float norm_t4 = (t4 - minTemp) / (maxTemp - minTemp);
			//wxLogError(wxString::Format("%g", node[0][0]));
			
			norm_t1 = norm_t1 > 1.0f ? 1.0f : norm_t1;
			norm_t2 = norm_t2 > 1.0f ? 1.0f : norm_t2;
			norm_t3 = norm_t3 > 1.0f ? 1.0f : norm_t3;
			norm_t4 = norm_t4 > 1.0f ? 1.0f : norm_t4;

			norm_t1 = norm_t1 < 0.0f ? 0.0f : norm_t1;
			norm_t2 = norm_t2 < 0.0f ? 0.0f : norm_t2;
			norm_t3 = norm_t3 < 0.0f ? 0.0f : norm_t3;
			norm_t4 = norm_t4 < 0.0f ? 0.0f : norm_t4;
			
			// indeksowanie znormalizowanych kolorów
			int colorId1 = (int)(norm_t1 * (numColors - 1));
			int colorId2 = (int)(norm_t2 * (numColors - 1));
			int colorId3 = (int)(norm_t3 * (numColors - 1));
			int colorId4 = (int)(norm_t4 * (numColors - 1));

			// Rysowanie kwadratami GL_QUADS
			glColor3f(colorMap[colorId1][0], colorMap[colorId1][1], colorMap[colorId1][2]);
			glVertex2f(-1.0f + 2.0f * (float)k / (N - 1), -1.0f + 2.0f * (float)l / (M - 1));

			glColor3f(colorMap[colorId2][0], colorMap[colorId2][1], colorMap[colorId2][2]);
			glVertex2f(-1.0f + 2.0f * (float)(k + 1) / (N - 1), -1.0f + 2.0f * (float)l / (M - 1));

			glColor3f(colorMap[colorId3][0], colorMap[colorId3][1], colorMap[colorId3][2]);
			glVertex2f(-1.0f + 2.0f * (float)(k + 1) / (N - 1), -1.0f + 2.0f * (float)(l + 1) / (M - 1));

			glColor3f(colorMap[colorId4][0], colorMap[colorId4][1], colorMap[colorId4][2]);
			glVertex2f(-1.0f + 2.0f * (float)k / (N - 1), -1.0f + 2.0f * (float)(l + 1) / (M - 1));
		}
	}

	glEnd();
}
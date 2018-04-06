#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#ifndef CustomCanvas_h
#define CustomCanvas_h

/*
  HOW TO USE
  Declare your TCanvas as a CustomCanvas instead.  It'll keep track of the file names you save with that canvas.  To retrieve them call GetPointer().
  The output of GetPointer() easily usable for the input of TexSlides.C
 */

class CustomCanvas : public TCanvas{
 private:
  std::vector<std::string>* histNames=new std::vector<std::string>();

 public:
  std::vector<std::string>* GetPointer()
    {
      return histNames;
    }
  void SaveAs(const char* filename="", Option_t* option="")
  {
    std::string name=(std::string)filename;
    histNames->push_back(name);

    if (name.length()>0) TCanvas::SaveAs(filename, option);
  }
  CustomCanvas(const char* name, const char* title, Int_t ww, Int_t wh)
    : TCanvas(name, title, ww, wh)
    {
    }
};

#endif

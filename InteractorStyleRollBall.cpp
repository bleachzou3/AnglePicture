#include "InteractorStyleRollBall.h"
#include <vtkObjectFactory.h>
#include <vtkRenderWindowInteractor.h>
#include <iostream>
#include <string>
vtkStandardNewMacro(InteractorStyleRollBall);
InteractorStyleRollBall::InteractorStyleRollBall()
{
	ctrlPressed = false;
	centerLineOnlyRenderer = 0;
	vascularRenderer = 0;
}
void InteractorStyleRollBall::OnLeftButtonUp()
{
	
	vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
}
void InteractorStyleRollBall::OnLeftButtonDown()
{

	vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void InteractorStyleRollBall::OnKeyPress()
{
	  vtkRenderWindowInteractor *rwi = this->Interactor;
      std::string key = rwi->GetKeySym();
      if(key == "Control_L")
	  {
		  ctrlPressed = true;
	  }
      // Output the key that was pressed
     vtkInteractorStyleTrackballCamera::OnKeyPress();
}

void InteractorStyleRollBall::OnKeyRelease()
{
	  vtkRenderWindowInteractor *rwi = this->Interactor;
      std::string key = rwi->GetKeySym();
	  if(key == "Control_L")
	  {
		  ctrlPressed = false;
	  }
	   vtkInteractorStyleTrackballCamera::OnKeyRelease();
}

void InteractorStyleRollBall::OnMouseMove()
{
	if(this->GetCurrentRenderer() == centerLineOnlyRenderer && ctrlPressed == true)
	{
		cout << "h" << endl;
	}
	vtkInteractorStyleTrackballCamera::OnMouseMove();
}

void InteractorStyleRollBall::setCenterLineOnlyRenderer(vtkRenderer* _centerLineOnlyRenderer)
{
	centerLineOnlyRenderer = _centerLineOnlyRenderer;
}

void InteractorStyleRollBall::setVascularRenderer(vtkRenderer* _vascularRenderer)
{
	vascularRenderer = _vascularRenderer;
}
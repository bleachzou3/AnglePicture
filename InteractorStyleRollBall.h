#ifndef INTERACTOR_STYLE_ROLL_BALL_H_
#define INTERACTOR_STYLE_ROLL_BALL_H_
#include<vtkInteractorStyleTrackballCamera.h>
#include<vtkRenderer.h>
#include<vtkSetGet.h>
class InteractorStyleRollBall:public vtkInteractorStyleTrackballCamera
{
public:
	static InteractorStyleRollBall* New();
	vtkTypeMacro(InteractorStyleRollBall,vtkInteractorStyleTrackballCamera);

	InteractorStyleRollBall();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnKeyPress();
	virtual void OnKeyRelease();
	virtual void OnMouseMove();
	void setCenterLineOnlyRenderer(vtkRenderer* _centerLineOnlyRenderer);
	void setVascularRenderer(vtkRenderer* _vascularRenderer);
private:
    bool ctrlPressed;
	vtkRenderer* centerLineOnlyRenderer;
	vtkRenderer* vascularRenderer;
};

#endif
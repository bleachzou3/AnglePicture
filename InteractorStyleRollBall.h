#ifndef INTERACTOR_STYLE_ROLL_BALL_H_
#define INTERACTOR_STYLE_ROLL_BALL_H_
#include<vtkInteractorStyleTrackballCamera.h>
#include<vtkRenderer.h>
#include<vtkSetGet.h>
#include<vtkSmartPointer.h>
#include<vtkPointPicker.h>
#include<vtkPolyData.h>
class InteractorStyleRollBall:public vtkInteractorStyleTrackballCamera
{
public:
	static InteractorStyleRollBall* New();
	vtkTypeMacro(InteractorStyleRollBall,vtkInteractorStyleTrackballCamera);

	InteractorStyleRollBall();
	~InteractorStyleRollBall();
	virtual void OnLeftButtonDown();
	virtual void OnLeftButtonUp();
	virtual void OnKeyPress();
	virtual void OnKeyRelease();
	virtual void OnMouseMove();
	void setCenterLineOnlyRenderer(vtkRenderer* _centerLineOnlyRenderer);
	void setVascularRenderer(vtkRenderer* _vascularRenderer);
	void setCenterLineData(vtkPolyData* _centerLineData);
private:
    bool ctrlPressed;
	//从外部传进来的
	vtkRenderer* centerLineOnlyRenderer;

	//从外部传进来的
	vtkRenderer* vascularRenderer;


	//从外表传入,这个是已经计算出了normal的中心线
	vtkPolyData* centerLineData;



	//构造函数里生成
	vtkSmartPointer<vtkPointPicker> PointPicker;
	

	//显示的小球
	vtkActor* indicatorBall;

	vtkIdType curId;
private:
	void showBall(vtkIdType _id);
};

#endif
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
	//���ⲿ��������
	vtkRenderer* centerLineOnlyRenderer;

	//���ⲿ��������
	vtkRenderer* vascularRenderer;


	//�������,������Ѿ��������normal��������
	vtkPolyData* centerLineData;



	//���캯��������
	vtkSmartPointer<vtkPointPicker> PointPicker;
	

	//��ʾ��С��
	vtkActor* indicatorBall;

	vtkIdType curId;
private:
	void showBall(vtkIdType _id);
};

#endif
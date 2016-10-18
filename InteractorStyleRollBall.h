#ifndef INTERACTOR_STYLE_ROLL_BALL_H_
#define INTERACTOR_STYLE_ROLL_BALL_H_
#include<vtkInteractorStyleTrackballCamera.h>
#include<vtkRenderer.h>
#include<vtkSetGet.h>
#include<vtkSmartPointer.h>
#include<vtkPointPicker.h>
#include<vtkPolyData.h>
#include<vtkImageData.h>
#include<vtkImageActor.h>
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
	void setImageActor(vtkImageActor*_imageActor);
	void setCurrentOblique(vtkImageData* _currentOblique);
	void setOriginalImage(vtkImageData* _originalImage);
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

	//ԭʼ�����ݼ��������洫����
	vtkImageData* originalImage;

	//��ǰ��Ƭ�����ݣ�����Ƚ�������Ȼ�����洫����������Ҫ����������
	vtkImageData* currentOblique;

	//�����洫����
	vtkImageActor* imageActor;
private:
	void showBall(vtkIdType _id);
};

#endif
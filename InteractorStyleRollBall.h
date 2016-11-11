#ifndef INTERACTOR_STYLE_ROLL_BALL_H_
#define INTERACTOR_STYLE_ROLL_BALL_H_
using namespace std;
#include<vtkInteractorStyleTrackballCamera.h>
#include<vtkRenderer.h>
#include<vtkSetGet.h>
#include<vtkSmartPointer.h>
#include<vtkPointPicker.h>
#include<vtkPolyData.h>
#include<vtkImageData.h>
#include<vtkImageActor.h>
#include<unordered_map>
#include<vtkCursor2D.h>
#include<vtkCamera.h>
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
	void setImageRenderer(vtkRenderer* _imageRenderer);
	void setAllAnglesImages(unordered_map<vtkIdType,vtkImageData*>& _allImages);
private:
    bool ctrlPressed;

	bool leftButtonDown;
	//���ⲿ��������
	vtkRenderer* centerLineOnlyRenderer;

	//���ⲿ��������
	vtkRenderer* vascularRenderer;

	//��ʾ��Ƭ���ݵ�renderer
	vtkRenderer* imageRenderer;

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

	//������нǶȵ�ͼƬ,�������ʵʱ����İ汾�������ע�͵�
	std::unordered_map<vtkIdType,vtkImageData*> allAnglesImages;


	bool calculateAllAnglesImages;

	//��ʾcursor
	vtkSmartPointer<vtkCursor2D> cursor;

	vtkSmartPointer<vtkPolyDataMapper> cursorMapper;

	vtkSmartPointer<vtkActor> cursorActor;

	bool showCursor;
private:
	void showBall(vtkIdType _id);

	//��showBall�������changePicture
	void changePicture(vtkIdType _id,double* fixedPoint);

	//��image renderer����������Ѫ��λ��
	void initLesionCursor();

	//ʵ������ͬ������
	void simultaneousOrientation();

	//��ֻ��ʾ�����ߵ�renderer�����ڵ�������point��id,���û�е����ã��ͷ���-1;
	vtkIdType getCenterLinePointId();
};

#endif
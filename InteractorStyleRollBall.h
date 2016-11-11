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
	//从外部传进来的
	vtkRenderer* centerLineOnlyRenderer;

	//从外部传进来的
	vtkRenderer* vascularRenderer;

	//显示切片数据的renderer
	vtkRenderer* imageRenderer;

	//从外表传入,这个是已经计算出了normal的中心线
	vtkPolyData* centerLineData;



	//构造函数里生成
	vtkSmartPointer<vtkPointPicker> PointPicker;
	

	//显示的小球
	vtkActor* indicatorBall;

	vtkIdType curId;

	//原始的数据集，从外面传进来
	vtkImageData* originalImage;

	//当前切片的数据，这个比较特殊虽然从外面传进来，但是要在这里销毁
	vtkImageData* currentOblique;

	//从外面传进来
	vtkImageActor* imageActor;

	//存放所有角度的图片,如果不是实时计算的版本这个可以注释掉
	std::unordered_map<vtkIdType,vtkImageData*> allAnglesImages;


	bool calculateAllAnglesImages;

	//显示cursor
	vtkSmartPointer<vtkCursor2D> cursor;

	vtkSmartPointer<vtkPolyDataMapper> cursorMapper;

	vtkSmartPointer<vtkActor> cursorActor;

	bool showCursor;
private:
	void showBall(vtkIdType _id);

	//在showBall里面调用changePicture
	void changePicture(vtkIdType _id,double* fixedPoint);

	//在image renderer里面标记现在血管位置
	void initLesionCursor();

	//实现两边同步操作
	void simultaneousOrientation();

	//在只显示中心线的renderer区域内点击，获得point的id,如果没有点击获得，就返回-1;
	vtkIdType getCenterLinePointId();
};

#endif
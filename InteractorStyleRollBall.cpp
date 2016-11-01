#include "InteractorStyleRollBall.h"
#include <vtkObjectFactory.h>
#include <vtkRenderWindowInteractor.h>
#include <iostream>
#include <string>
#include <vtkRenderWindow.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageMapper3D.h>
#include "AnglePictureUtility.h"
#include <vtkJPEGWriter.h>
#include <vtkImageWriter.h>
#include <vtkDoubleArray.h>
#include <vtkCursor2D.h>
#include <vtkProperty.h>
vtkStandardNewMacro(InteractorStyleRollBall);
InteractorStyleRollBall::InteractorStyleRollBall()
{
	ctrlPressed = false;
	centerLineOnlyRenderer = 0;
	vascularRenderer = 0;
	indicatorBall = 0;
	curId = -1;
	PointPicker = vtkSmartPointer<vtkPointPicker>::New();
	imageActor = 0;
	currentOblique = 0;
	originalImage = 0;
	initLesionCursor();

}


void InteractorStyleRollBall::setImageRenderer(vtkRenderer* _imageRenderer)
{
	imageRenderer = _imageRenderer;

}
void InteractorStyleRollBall::setCurrentOblique(vtkImageData*_currentOblique)
{
	currentOblique = _currentOblique;
}

void InteractorStyleRollBall::setImageActor(vtkImageActor* _imageActor)
{
	imageActor = _imageActor;
}

void InteractorStyleRollBall::setOriginalImage(vtkImageData* _originalImage)
{
	originalImage = _originalImage;
}
void InteractorStyleRollBall::OnLeftButtonUp()
{
	
	vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
}
void InteractorStyleRollBall::OnLeftButtonDown()
{

	int x = this->Interactor->GetEventPosition()[0];
    int y = this->Interactor->GetEventPosition()[1];
	
	//switch current renderer
    this->FindPokedRenderer(x, y);
	if(this->GetCurrentRenderer() == imageRenderer)
	{
		if(showCursor)
		{
		   imageRenderer->AddActor(cursorActor);
		   showCursor = false;
		}else
		{  
			imageRenderer->RemoveActor(cursorActor);
			showCursor = true;
		}
	}

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
		  if(indicatorBall != 0)
		  {
			  vascularRenderer->RemoveActor(indicatorBall);
			  vascularRenderer->Render();
			  vascularRenderer->GetRenderWindow()->Render();
			  if(indicatorBall != 0 && indicatorBall->GetReferenceCount() == 1)
			  {
				indicatorBall->Delete();
				curId = -1;
			  }

		  }
	  }
	   vtkInteractorStyleTrackballCamera::OnKeyRelease();
}

void InteractorStyleRollBall::OnMouseMove()
{
    int x = this->Interactor->GetEventPosition()[0];
    int y = this->Interactor->GetEventPosition()[1];
	
	//switch current renderer
    this->FindPokedRenderer(x, y);
	if(this->GetCurrentRenderer() == centerLineOnlyRenderer && ctrlPressed == true)
	{
	   this->PointPicker->Pick(this->Interactor->GetEventPosition()[0],
                 this->Interactor->GetEventPosition()[1],
                 0,  // always zero.
				 this->GetCurrentRenderer());
	    if(this->PointPicker->GetPointId() >= 0)
		{
			vtkIdType id = this->PointPicker->GetPointId();
			showBall(id);
			cout << "x:" << x <<"y:" << y << endl;
			
		}

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

void InteractorStyleRollBall::setCenterLineData(vtkPolyData* _centerlineData)
{
	centerLineData = _centerlineData;
}
/**
*这个点的切片实时计算的版本
*/
void InteractorStyleRollBall::changePicture(vtkIdType _id,double* fixedPoint)
{
	
	vtkImageData *oldOblique = currentOblique;
    currentOblique = vtkImageData::New();	
	double* normal = centerLineData->GetPointData()->GetArray(AnglePictureUtility::OBLIQUE_NORMAL.c_str())->GetTuple(_id);
	Vector3 normalV(normal[0],normal[1],normal[2]);
	Vector3 fixedPointV(fixedPoint[0],fixedPoint[1],fixedPoint[2]);
	AnglePictureUtility::computeOblique(originalImage,normalV,fixedPointV,currentOblique);	
	imageActor->GetMapper()->SetInputData(currentOblique);
	imageActor->Update();
	imageRenderer->AddActor(imageActor);
    int extent[6];
	int dimension[3];
	double origin[3];
	double spacing[3];
	currentOblique->GetExtent(extent);
	currentOblique->GetDimensions(dimension);
	currentOblique->GetOrigin(origin);
	currentOblique->GetSpacing(spacing);
	cout << "currentOblique:" << extent[0] << " " <<extent[1]<< " " << extent[2] <<" " << extent[3] <<"  " << extent[4] << "  " <<extent[5] << endl;
	cout << "currentOblique:dimension" << dimension[0] << " " << dimension[1] << " " << dimension[2] << endl;
	cout << "currentOblique: origin:" << origin[0] << " " << origin[1] << " " << origin[2] << endl;
	cout << "currentOblique:spacing" << spacing[0] << " " << spacing[1] << " " << spacing[2] << endl;


	if(oldOblique != 0 && oldOblique->GetReferenceCount() == 1)
	{
		oldOblique->Delete();
	}

	imageRenderer->Render();
	imageRenderer->GetRenderWindow()->Render();

}

void InteractorStyleRollBall::showBall(vtkIdType _id)
{
			if(_id == curId)
			{
				return;
			}
			cout << "id:" << _id << endl;
			
			curId = _id;
			//有小球存在，保留着
			if(indicatorBall != 0)
			{
				vascularRenderer->RemoveActor(indicatorBall);
			}
			vascularRenderer->Render();
			if(indicatorBall!= 0 && indicatorBall->GetReferenceCount() == 1)
			{
				indicatorBall->Delete();
			}

			indicatorBall = vtkActor::New();
            //获取球心坐标
			double center[3];
			centerLineData->GetPoint(_id,center);
			cout << "center: " <<center[0] <<" " << center[1] << "  " << center[2] << endl;

			//获取半径
			double* radius = centerLineData->GetPointData()->GetArray(AnglePictureUtility::MAXIMUM_INSCRIBED_SPHERE_RADIUS.c_str())->GetTuple(curId);
			cout << "radius:" << radius[0] << endl;

			//生成小球
			vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
			sphere->SetCenter(center[0],center[1],center[2]);
			sphere->SetRadius(radius[0]);
			sphere->Update();
			vtkSmartPointer<vtkPolyDataMapper> sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			//sphereMapper->SetInputConnection(sphereMapper->GetOutputPort());
			sphereMapper->SetInputData(sphere->GetOutput());
			indicatorBall->SetMapper(sphereMapper);
			vascularRenderer->AddActor(indicatorBall);

			changePicture(_id,center);

			imageRenderer->Render();
			vascularRenderer->Render();
			vascularRenderer->GetRenderWindow()->Render();
			imageRenderer->GetRenderWindow()->Render();
			
}

InteractorStyleRollBall::~InteractorStyleRollBall()
{

	if(indicatorBall != 0 && indicatorBall->GetReferenceCount() == 1)
	{
		indicatorBall->Delete();
	}
	if(imageActor != 0 && imageActor->GetReferenceCount() == 1)
	{
		imageActor->Delete();
	}
	if(currentOblique != 0 && currentOblique->GetReferenceCount() == 1)
	{
		currentOblique->Delete();
	}
	for(pair<vtkIdType,vtkImageData*> cur:allAnglesImages)
	{
		if(cur.second != 0 && cur.second->GetReferenceCount() == 1)
		{
			cur.second->Delete();
		}
	}
	
}

void InteractorStyleRollBall::setAllAnglesImages(unordered_map<vtkIdType,vtkImageData*> & _allImage)
{
	allAnglesImages = _allImage;
}


void InteractorStyleRollBall::initLesionCursor()
{
		cursor = vtkSmartPointer<vtkCursor2D>::New();
		cursor->SetModelBounds(-10,10,-10,10,0,0);
		cursor->AllOn();
		cursor->OutlineOff();
		cursor->Update();
 
	    cursorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		cursorMapper->SetInputConnection(cursor->GetOutputPort());
		cursorActor = vtkSmartPointer<vtkActor>::New();
		cursorActor->GetProperty()->SetColor(1,0,0);
		cursorActor->SetMapper(cursorMapper);

		showCursor = true;
}
/**
*这个是把所有点切片都都计算出来版本(实时显示图像)
*
void InteractorStyleRollBall::changePicture(vtkIdType _id,double* fixedPoint)
{
	imageActor->GetMapper()->SetInputData(allAnglesImages[_id]);
	imageActor->Update();
	//imageRenderer->AddActor(imageActor);
	imageRenderer->Render();
	imageRenderer->GetRenderWindow()->Render();

}
*/
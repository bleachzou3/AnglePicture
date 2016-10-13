#include "InteractorStyleRollBall.h"
#include <vtkObjectFactory.h>
#include <vtkRenderWindowInteractor.h>
#include <iostream>
#include <string>
#include <vtkRenderWindow.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
vtkStandardNewMacro(InteractorStyleRollBall);
InteractorStyleRollBall::InteractorStyleRollBall()
{
	ctrlPressed = false;
	centerLineOnlyRenderer = 0;
	vascularRenderer = 0;
	indicatorBall = 0;
	curId = -1;
	PointPicker = vtkSmartPointer<vtkPointPicker>::New();
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
			double* radius = centerLineData->GetPointData()->GetArray("MaximumInscribedSphereRadius")->GetTuple(curId);
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
			vascularRenderer->Render();
			vascularRenderer->GetRenderWindow()->Render();
			
}

InteractorStyleRollBall::~InteractorStyleRollBall()
{
	if(indicatorBall != 0 && indicatorBall->GetReferenceCount() == 1)
	{
		indicatorBall->Delete();
	}
}
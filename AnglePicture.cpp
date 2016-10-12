#include <vtkJPEGReader.h>
#include <vtkImageMapper3D.h>
#include <vtkImageActor.h> // Note: this is a 3D actor (c.f. vtkImageMapper which is 2D)
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkSmartPointer.h>
#include <vtkMatrix4x4.h>
#include <vtkImageReslice.h>
#include <vtkDataObject.h>
#include <vtkImageData.h>
#include "Vector3.h"
#include "ONB.h"
#include <vtkDICOMImageReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include "InteractorStyleRollBall.h"
/**
*@param original 原始的图像切片, direction 中心线的方向, fixedPoint图像切片中心经过的点
*
*/
void computeOblique(vtkImageData* original,Vector3 direction,Vector3 fixedPoint,vtkImageData*res)
{
	
	if(res == 0)
	{
		return;
	}
	ONB onb;
	onb.initFromU(direction);
	Vector3 U = onb.u();
	Vector3 V = onb.v();
	Vector3 W = onb.w();

	static double axialElements[16] = {
		V.x(),W.x(),U.x(),fixedPoint.x(),
		V.y(),W.y(),U.y(),fixedPoint.y(),
		V.z(),W.z(),U.z(),fixedPoint.z(),
		0,    0,    0,    1
	};

	vtkSmartPointer<vtkMatrix4x4> resliceAxes = vtkSmartPointer<vtkMatrix4x4>::New();
	resliceAxes->DeepCopy(axialElements);

	vtkSmartPointer<vtkImageReslice> reslice = vtkSmartPointer<vtkImageReslice>::New();
	reslice->SetInputData(original);
	reslice->SetOutputDimensionality(2);
	reslice->SetResliceAxes(resliceAxes);
	reslice->SetInterpolationModeToLinear();
	reslice->Update();

	res->DeepCopy(reslice->GetOutput());
}

/**
*计算每个切面的法向量obliqueNormal
* 计算左右相邻的两个点，哪个点近就用哪个
*/


void computeNormalByPoints(vtkPolyData * data)
{
	vtkDataArray* normals = data->GetPointData()->GetArray("obliqueNormal");

	if(normals != 0)
	{
		return;
	}
	int N = data->GetNumberOfPoints();
	int line = data->GetNumberOfLines();
	if(N < 3)
	{
		return;
	}
	vtkSmartPointer<vtkFloatArray> newNormals = vtkSmartPointer<vtkFloatArray>::New();
	newNormals->SetName("obliqueNormal");

	newNormals->SetNumberOfComponents(3);
	newNormals->SetNumberOfTuples(N);

	double firstNode[3]; 
	data->GetPoint(0,firstNode);
	double secondNode[3];
	data->GetPoint(1,secondNode);
	newNormals->SetComponent(0,0,firstNode[0] - secondNode[0]);
	newNormals->SetComponent(0,1,firstNode[1] - secondNode[1]);
	newNormals->SetComponent(0,2,firstNode[2] - secondNode[2]);

	for(int i = 1; i < N - 1; i++)
	{
		double currentNode[3];
		data->GetPoint(i,currentNode);
		double previousNode[3];
		double nextNode[3] ;
		data->GetPoint(i - 1,previousNode);
		data->GetPoint(i + 1,nextNode);

		double prevCurr = sqrt(vtkMath::Distance2BetweenPoints(currentNode,previousNode));
		double nextCurr = sqrt(vtkMath::Distance2BetweenPoints(currentNode,nextNode));

		double curNormal[3];
		if(prevCurr < nextCurr)
		{
			curNormal[0] = previousNode[0] - currentNode[0];
			curNormal[1] = previousNode[1] - currentNode[1];
			curNormal[2] = previousNode[2] - currentNode[2];
		}else
		{
			curNormal[0] = currentNode[0] - nextNode[0];
			curNormal[1] = currentNode[1] - nextNode[1];
			curNormal[2] = currentNode[2] - nextNode[2];
		}
		newNormals->SetComponent(i,0,curNormal[0]);
		newNormals->SetComponent(i,1,curNormal[1]);
	    newNormals->SetComponent(i,2,curNormal[2]);
	}

	double* rFirstNode = data->GetPoint(N - 1);
	double* rSecondNode = data->GetPoint(N - 2);

    newNormals->SetComponent(N - 1,0,firstNode[0] - secondNode[0]);
	newNormals->SetComponent(N - 1,1,firstNode[1] - secondNode[1]);
	newNormals->SetComponent(N - 1,2,firstNode[2] - secondNode[2]);

	data->GetPointData()->AddArray(newNormals);

	
}

/**
*根据图元进行计算
*/


int main(int argc, char* argv[])
{

   // Setup render window
  vtkSmartPointer<vtkRenderWindow> window = 
    vtkSmartPointer<vtkRenderWindow>::New();
    // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(window);
  vtkSmartPointer<InteractorStyleRollBall> interactorStyle = vtkSmartPointer<InteractorStyleRollBall>::New();
  interactor->SetInteractorStyle(interactorStyle);

  
  double tangentPicture[4] = {0.0,0.0,0.33,0.5};
  double modelRegin[4]= {0.33,0.5,1.0,1.0};
  double centerlineOnly[4] = {0.0,0.5,0.33,1};


  /*
  vtkSmartPointer<vtkDICOMImageReader> reader = vtkSmartPointer<vtkDICOMImageReader>::New();
  reader->SetDirectoryName("E://patientData//WU_AMAO");
  reader->Update();
  */
  //读取原始的图像切片数据
   vtkSmartPointer<vtkXMLImageDataReader> reader =
    vtkSmartPointer<vtkXMLImageDataReader>::New();
  reader->SetFileName("E:\\image_volume_voi.vti");
  reader->Update();
  vtkImageData* original = reader->GetOutput();
  vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
  //centerline
  Vector3 direction(2,3,1);
  int extent[6];
  double spacing[3];
  double origin[3];
  original->GetExtent(extent);
  original->GetSpacing(spacing);
  original->GetOrigin(origin);
  double center[3];
  center[0] = origin[0] + spacing[0] * 0.5 * (extent[0] + extent[1]);
  center[1] = origin[1] + spacing[1] * 0.5 * (extent[2] + extent[3]);
  center[2] = origin[2] + spacing[2] * 0.5 * (extent[4] + extent[5]);
  Vector3 fixedPoint(center[0],center[1],center[2]);

  computeOblique(original,direction,fixedPoint,data);

  int extentr[6];
  double spacingr[6];
  double originr[6];

  data->GetExtent(extentr);
  data->GetSpacing(spacingr);
  data->GetOrigin(originr);

  // Setup renderer展示图像的render
  vtkSmartPointer<vtkRenderer> renderer = 
    vtkSmartPointer<vtkRenderer>::New();
  renderer->SetViewport(tangentPicture);
   renderer->SetBackground(0.5,0.5,0.5);
     vtkSmartPointer<vtkImageActor> actor = 
    vtkSmartPointer<vtkImageActor>::New();
	   actor->GetMapper()->SetInputData(data);
  renderer->AddActor(actor);
  renderer->ResetCamera();  
  window->AddRenderer(renderer);
  renderer->SetViewPoint(tangentPicture);
  actor->Update();
























    // Read all the data from the file读取血管模型数据
  vtkSmartPointer<vtkXMLPolyDataReader> readerModel =
    vtkSmartPointer<vtkXMLPolyDataReader>::New();
  readerModel->SetFileName("E:\\model0927smooth.vtp");
  readerModel->Update();

  //读取中心线模型数据
  vtkSmartPointer<vtkXMLPolyDataReader> centerlineModelReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  centerlineModelReader->SetFileName("E:\\model0927centerline.vtp");
  centerlineModelReader->Update();
  vtkPolyData* centerline = centerlineModelReader->GetOutput();
  computeNormalByPoints(centerline);

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(readerModel->GetOutputPort());
  vtkSmartPointer<vtkPolyDataMapper> centerlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
 //centerlineMapper->SetInputConnection(centerlineModelReader->GetOutputPort());
  centerlineMapper->SetInputData(centerline);
 
  vtkSmartPointer<vtkActor> actorModel =
    vtkSmartPointer<vtkActor>::New();
  actorModel->GetProperty()->SetOpacity(0.1);
  actorModel->GetProperty()->SetColor(0.5,0.5,0.5);
  actorModel->SetMapper(mapper);

  vtkSmartPointer<vtkActor> centerlineModel = vtkSmartPointer<vtkActor>::New();
  centerlineModel->SetMapper(centerlineMapper);
  
  vtkSmartPointer<vtkRenderer> rendererModel =
    vtkSmartPointer<vtkRenderer>::New();
  double back[3] = {1.0,1.0,1.0};
  rendererModel->SetBackground(back);
  rendererModel->SetViewport(modelRegin);
  rendererModel->AddActor(actorModel);
  rendererModel->AddActor(centerlineModel);
  window->AddRenderer(rendererModel);




  //中心线数据我不拷贝，下面是单独显示中心线的操作
  vtkSmartPointer<vtkPolyDataMapper> centerLineMapperOnly = vtkSmartPointer<vtkPolyDataMapper>::New();
  centerLineMapperOnly->SetInputData(centerline);
  vtkSmartPointer<vtkActor> centerlineModelOnlyActor = vtkSmartPointer<vtkActor>::New();
  centerlineModelOnlyActor->SetMapper(centerLineMapperOnly);
  vtkSmartPointer<vtkRenderer> centerlineOnlyRender = vtkSmartPointer<vtkRenderer>::New();
  centerlineOnlyRender->AddActor(centerlineModelOnlyActor);
  centerlineOnlyRender->SetViewport(centerlineOnly);
  window->AddRenderer(centerlineOnlyRender);

  interactor->Start();
  return EXIT_SUCCESS;
}

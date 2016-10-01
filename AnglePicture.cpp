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

int main(int argc, char* argv[])
{

   // Setup render window
  vtkSmartPointer<vtkRenderWindow> window = 
    vtkSmartPointer<vtkRenderWindow>::New();
    // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(window);
  double tangentPicture[4] = {0.0,0.0,0.33,0.5};
  double modelRegin[4]= {0.33,0.5,1.0,1.0};


  /*
  vtkSmartPointer<vtkDICOMImageReader> reader = vtkSmartPointer<vtkDICOMImageReader>::New();
  reader->SetDirectoryName("E://patientData//WU_AMAO");
  reader->Update();
  */
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

  // Setup renderer
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
























    // Read all the data from the file
  vtkSmartPointer<vtkXMLPolyDataReader> readerModel =
    vtkSmartPointer<vtkXMLPolyDataReader>::New();
  readerModel->SetFileName("E:\\model0927.vtp");
  readerModel->Update();
 
  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(readerModel->GetOutputPort());
 
  vtkSmartPointer<vtkActor> actorModel =
    vtkSmartPointer<vtkActor>::New();
  actorModel->SetMapper(mapper);
 
  vtkSmartPointer<vtkRenderer> rendererModel =
    vtkSmartPointer<vtkRenderer>::New();

  rendererModel->SetViewport(modelRegin);
  rendererModel->AddActor(actorModel);
  window->AddRenderer(rendererModel);

  interactor->Start();
  return EXIT_SUCCESS;
}

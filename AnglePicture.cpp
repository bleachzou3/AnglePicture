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
		U.x(),V.x(),W.x(),fixedPoint.x(),
		U.y(),V.y(),W.y(),fixedPoint.y(),
		U.z(),V.z(),W.z(),fixedPoint.z(),
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


  // Read the image
  /*
  vtkSmartPointer<vtkJPEGReader> reader = 
      vtkSmartPointer<vtkJPEGReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();
  */
  // Create an actor
  //prepare for reading picture
  vtkSmartPointer<vtkDICOMImageReader> reader = vtkSmartPointer<vtkDICOMImageReader>::New();
  reader->SetDirectoryName("E://patientData//WU_AMAO");
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








  vtkSmartPointer<vtkImageActor> actor = 
    vtkSmartPointer<vtkImageActor>::New();
  actor->GetMapper()->SetInputData(data);
  actor->Update();
  // Setup renderer
  vtkSmartPointer<vtkRenderer> renderer = 
    vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(actor);
  renderer->ResetCamera();

  // Setup render window
  vtkSmartPointer<vtkRenderWindow> window = 
    vtkSmartPointer<vtkRenderWindow>::New();
  window->AddRenderer(renderer);

  // Setup render window interactor
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = 
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(window);

  // Setup interactor style (this is what implements the zooming, panning and brightness adjustment functionality)
  vtkSmartPointer<vtkInteractorStyleImage> style = 
    vtkSmartPointer<vtkInteractorStyleImage>::New();
  interactor->SetInteractorStyle(style);
  
  // Render and start interaction
  interactor->Start();

  return EXIT_SUCCESS;
}

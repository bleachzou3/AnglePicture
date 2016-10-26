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
#include "AnglePictureUtility.h"
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
  //vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
  vtkImageData* data = vtkImageData::New();
  //centerline
























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
  AnglePictureUtility::computeNormalByPoints(centerline);

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











  //中心点   
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
  Vector3 direction(0,0,1);
  AnglePictureUtility::computeOblique(original,direction,fixedPoint,data);
  unordered_map<vtkIdType,vtkImageData*> allImage;
  AnglePictureUtility::computeAllAngleImages(centerline,allImage,original);

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
 //    vtkSmartPointer<vtkImageActor> actor = 
 //   vtkSmartPointer<vtkImageActor>::New();
       vtkImageActor*actor = vtkImageActor::New();
	   actor->GetMapper()->SetInputData(data);
	   actor->Update();

  renderer->AddActor(actor);
  renderer->ResetCamera();  
  window->AddRenderer(renderer);
  renderer->SetViewPoint(tangentPicture);
  actor->Update();




















  //中心线数据我不拷贝，下面是单独显示中心线的操作
  vtkSmartPointer<vtkPolyDataMapper> centerLineMapperOnly = vtkSmartPointer<vtkPolyDataMapper>::New();
  centerLineMapperOnly->SetInputData(centerline);
  vtkSmartPointer<vtkActor> centerlineModelOnlyActor = vtkSmartPointer<vtkActor>::New();
  centerlineModelOnlyActor->SetMapper(centerLineMapperOnly);
  vtkSmartPointer<vtkRenderer> centerlineOnlyRender = vtkSmartPointer<vtkRenderer>::New();
  centerlineOnlyRender->AddActor(centerlineModelOnlyActor);
  centerlineOnlyRender->SetViewport(centerlineOnly);
  window->AddRenderer(centerlineOnlyRender);




  interactorStyle->setCenterLineOnlyRenderer(centerlineOnlyRender);
  interactorStyle->setVascularRenderer(rendererModel);
  interactorStyle->setImageRenderer(renderer);
  interactorStyle->setCenterLineData(centerline);
  interactorStyle->setOriginalImage(original);
  interactorStyle->setImageActor(actor);
  interactorStyle->setCurrentOblique(data);

  //不是把所有点的切片计算出来，可以把这个注释掉
  //interactorStyle->setAllAnglesImages(allImage);

  interactor->SetInteractorStyle(interactorStyle);
  interactor->Start();
  if(data != 0 && data->GetReferenceCount() == 1)
  {
	  data->Delete();
  }
  return EXIT_SUCCESS;
}

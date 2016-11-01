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
#include <log4cpp/Category.hh>
#include <log4cpp/PropertyConfigurator.hh>
int main(int argc, char* argv[])
{


	try
	{
		log4cpp::PropertyConfigurator::configure("log4cpp.conf");
	}
	catch (log4cpp::ConfigureFailure& f)
	{
		std::cout << "Configure Problem " << f.what() << std::endl;
		return -1;
	}
	log4cpp::Category& rootLog  = log4cpp::Category::getRoot();
	rootLog.info("����ʼִ��");
	log4cpp::Category& subLog = log4cpp::Category::getInstance(std::string("sub1"));
	subLog.info("����ʼִ��");
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
  //��ȡԭʼ��ͼ����Ƭ����
   vtkSmartPointer<vtkXMLImageDataReader> originalImageDataReader =
    vtkSmartPointer<vtkXMLImageDataReader>::New();
  originalImageDataReader->SetFileName("E:\\image_volume_voi.vti");
  originalImageDataReader->Update();
  vtkImageData* original = originalImageDataReader->GetOutput();
  //vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
  vtkImageData* initialObliqueData = vtkImageData::New();
  //centerline
























    // Read all the data from the file��ȡѪ��ģ������
  vtkSmartPointer<vtkXMLPolyDataReader> vascularModelReader =
    vtkSmartPointer<vtkXMLPolyDataReader>::New();
  vascularModelReader->SetFileName("E:\\model0927smooth.vtp");
  vascularModelReader->Update();

  //��ȡ������ģ������
  vtkSmartPointer<vtkXMLPolyDataReader> centerlineModelReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  centerlineModelReader->SetFileName("E:\\model0927centerline.vtp");
  centerlineModelReader->Update();
  vtkPolyData* centerline = centerlineModelReader->GetOutput();
  AnglePictureUtility::computeNormalByPoints(centerline);

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> centerLineMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  centerLineMapper->SetInputConnection(vascularModelReader->GetOutputPort());
  vtkSmartPointer<vtkPolyDataMapper> centerlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
 //centerlineMapper->SetInputConnection(centerlineModelReader->GetOutputPort());
  centerlineMapper->SetInputData(centerline);
 
  vtkSmartPointer<vtkActor> actorModel =
    vtkSmartPointer<vtkActor>::New();
  actorModel->GetProperty()->SetOpacity(0.1);
  actorModel->GetProperty()->SetColor(0.5,0.5,0.5);
  actorModel->SetMapper(centerLineMapper);

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











  //���ĵ�   
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
  AnglePictureUtility::computeOblique(original,direction,fixedPoint,initialObliqueData);
  unordered_map<vtkIdType,vtkImageData*> allImage;
  AnglePictureUtility::computeAllAngleImages(centerline,allImage,original);

  int extentr[6];
  double spacingr[6];
  double originr[6];

  initialObliqueData->GetExtent(extentr);
  initialObliqueData->GetSpacing(spacingr);
  initialObliqueData->GetOrigin(originr);

  // Setup rendererչʾͼ���render
  vtkSmartPointer<vtkRenderer> obliqueImageRenderer = 
    vtkSmartPointer<vtkRenderer>::New();
  obliqueImageRenderer->SetViewport(tangentPicture);
   obliqueImageRenderer->SetBackground(0.5,0.5,0.5);
 //    vtkSmartPointer<vtkImageActor> actor = 
 //   vtkSmartPointer<vtkImageActor>::New();
       vtkImageActor*obliqueImageActor = vtkImageActor::New();
	   obliqueImageActor->GetMapper()->SetInputData(initialObliqueData);
	   obliqueImageActor->Update();

  obliqueImageRenderer->AddActor(obliqueImageActor);
  obliqueImageRenderer->ResetCamera();  
  window->AddRenderer(obliqueImageRenderer);
  obliqueImageRenderer->SetViewPoint(tangentPicture);
  obliqueImageActor->Update();

 


















  //�����������Ҳ������������ǵ�����ʾ�����ߵĲ���
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
  interactorStyle->setImageRenderer(obliqueImageRenderer);
  interactorStyle->setCenterLineData(centerline);
  interactorStyle->setOriginalImage(original);
  interactorStyle->setImageActor(obliqueImageActor);
  interactorStyle->setCurrentOblique(initialObliqueData);

  //���ǰ����е����Ƭ������������԰����ע�͵�
  //interactorStyle->setAllAnglesImages(allImage);

  interactor->SetInteractorStyle(interactorStyle);
  interactor->Start();
  if(initialObliqueData != 0 && initialObliqueData->GetReferenceCount() == 1)
  {
	  initialObliqueData->Delete();
  }
  return EXIT_SUCCESS;
}

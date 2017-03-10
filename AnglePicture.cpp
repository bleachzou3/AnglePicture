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
#include <vtkLabeledDataMapper.h>
#include <vtkPointSource.h>
#include <vtkActor2D.h>
#include <vtkCellArray.h>
#include "Common.h"
#include <vtkPolyLine.h>
#include <vtkLine.h>
#include <vtkTriangleFilter.h>
#include <vtkExtractEdges.h>
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
	subLog.info("程序开始执行......");
	rootLog.info("程序开始执行......");
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
  double addedInformation[4]={0.33,0.0,1.0,0.5};


  /*
  vtkSmartPointer<vtkDICOMImageReader> reader = vtkSmartPointer<vtkDICOMImageReader>::New();
  reader->SetDirectoryName("E://patientData//WU_AMAO");
  reader->Update();
  */
  //读取原始的图像切片数据
   vtkSmartPointer<vtkXMLImageDataReader> originalImageDataReader =
    vtkSmartPointer<vtkXMLImageDataReader>::New();
  originalImageDataReader->SetFileName("E:\\image_volume_voi.vti");
  originalImageDataReader->Update();
  vtkImageData* original = originalImageDataReader->GetOutput();
  //vtkSmartPointer<vtkImageData> data = vtkSmartPointer<vtkImageData>::New();
  vtkImageData* initialObliqueData = vtkImageData::New();
  //centerline
























    // Read all the data from the file读取血管模型数据
  vtkSmartPointer<vtkXMLPolyDataReader> vascularModelReader =
    vtkSmartPointer<vtkXMLPolyDataReader>::New();
  vascularModelReader->SetFileName("E:\\model0927smooth.vtp");
  vascularModelReader->Update();

  //读取中心线模型数据
  vtkSmartPointer<vtkXMLPolyDataReader> centerlineModelReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

  //数据中心线湿地9+
  centerlineModelReader->SetFileName("E:\\model0927centerline.vtp");
  //centerlineModelReader->SetFileName("E:\\VMTKCenterlinesOut.vtp");
  centerlineModelReader->Update();





  vtkPolyData* centerline = centerlineModelReader->GetOutput();
  cout << "numberofCells" <<centerline->GetNumberOfCells()<< endl;
  cout << "numberofPoints" << centerline->GetNumberOfPoints() << endl;
  centerline->GetLines()->InitTraversal();
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();

  //获得某一段内的所有点
  /**
  centerline->GetCellPoints(0,idList);
  for(vtkIdType pointId = 0; pointId < idList->GetNumberOfIds(); pointId++)
      {
      std::cout << idList->GetId(pointId) << " ";
      }
*/
  //遍历一个血管中心线所有图元，每一个图元是一个vtkPolyLine或者是vtkLine，然后遍历每个图元里的所有点
  /**
  while(centerline->GetLines()->GetNextCell(idList))
    {
    std::cout << "Line has " << idList->GetNumberOfIds() << " points." << std::endl;
 
    for(vtkIdType pointId = 0; pointId < idList->GetNumberOfIds(); pointId++)
      {
      std::cout << idList->GetId(pointId) << " ";
      }
    std::cout << std::endl;
    }
	*/
   

  //这种方法是失败的
  /**
    vtkSmartPointer<vtkTriangleFilter> triangleFilter =
      vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(centerline);
  triangleFilter->Update();

  vtkSmartPointer<vtkExtractEdges> extractEdges =
    vtkSmartPointer<vtkExtractEdges>::New();
  extractEdges->SetInputConnection(triangleFilter->GetOutputPort());
  extractEdges->Update();

  vtkSmartPointer<vtkPolyData> mesh = extractEdges->GetOutput();
  cout <<"mesh.............."<< mesh->GetNumberOfPoints()<<endl;
  cout <<"mesh............." << mesh->GetNumberOfCells()<<endl;
  cout << "mesh............." << mesh->GetCellType(2) << endl;
  */
   




  AnglePictureUtility::computeNormalByPoints(centerline);
  

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> vascularModelMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  vascularModelMapper->SetInputConnection(vascularModelReader->GetOutputPort());
  vtkSmartPointer<vtkPolyDataMapper> centerlineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
 //centerlineMapper->SetInputConnection(centerlineModelReader->GetOutputPort());
  centerlineMapper->SetInputData(centerline);
 
  vtkSmartPointer<vtkActor> actorVascularModel =
    vtkSmartPointer<vtkActor>::New();
  actorVascularModel->GetProperty()->SetOpacity(0.1);
  actorVascularModel->GetProperty()->SetColor(0.5,0.5,0.5);
  actorVascularModel->SetMapper(vascularModelMapper);

  vtkSmartPointer<vtkActor> centerlineModel = vtkSmartPointer<vtkActor>::New();
  centerlineModel->SetMapper(centerlineMapper);
  
  vtkSmartPointer<vtkRenderer> rendererModel =
    vtkSmartPointer<vtkRenderer>::New();
  double back[3] = {1.0,1.0,1.0};
  rendererModel->SetBackground(back);
  rendererModel->SetViewport(modelRegin);
  rendererModel->AddActor(actorVascularModel);
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
  AnglePictureUtility::computeOblique(original,direction,fixedPoint,initialObliqueData);
  unordered_map<vtkIdType,vtkImageData*> allImage;
  //计算所有角度的图片
 // AnglePictureUtility::computeAllAngleImages(centerline,allImage,original);

  int extentr[6];
  double spacingr[6];
  double originr[6];

  initialObliqueData->GetExtent(extentr);
  initialObliqueData->GetSpacing(spacingr);
  initialObliqueData->GetOrigin(originr);

  // Setup renderer展示图像的render
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

 


















  //中心线数据我不拷贝，下面是单独显示中心线的操作
  vtkSmartPointer<vtkPolyDataMapper> centerLineMapperOnly = vtkSmartPointer<vtkPolyDataMapper>::New();
  centerLineMapperOnly->SetInputData(centerline);
  vtkSmartPointer<vtkActor> centerlineModelOnlyActor = vtkSmartPointer<vtkActor>::New();
  centerlineModelOnlyActor->SetMapper(centerLineMapperOnly);
  vtkSmartPointer<vtkRenderer> centerlineOnlyRender = vtkSmartPointer<vtkRenderer>::New();
  centerlineOnlyRender->AddActor(centerlineModelOnlyActor);





  /**
  *
  vtkSmartPointer<vtkLabeledDataMapper> labelMapper = 
    vtkSmartPointer<vtkLabeledDataMapper>::New();
  labelMapper->SetInputData(centerline);
  vtkSmartPointer<vtkActor2D> labelActor = 
    vtkSmartPointer<vtkActor2D>::New();
  labelActor->SetMapper(labelMapper);
  centerlineOnlyRender->AddActor2D(labelActor);
  */

  centerlineOnlyRender->SetViewport(centerlineOnly);
  window->AddRenderer(centerlineOnlyRender);
  







  //在右下方格里面显示额外的信息
 









  
  interactorStyle->setCenterLineOnlyRenderer(centerlineOnlyRender);
  interactorStyle->setVascularRenderer(rendererModel);
  interactorStyle->setImageRenderer(obliqueImageRenderer);
  interactorStyle->setCenterLineData(centerline);
  interactorStyle->setOriginalImage(original);
  interactorStyle->setImageActor(obliqueImageActor);
  interactorStyle->setCurrentOblique(initialObliqueData);

  //不是把所有点的切片计算出来，可以把这个注释掉
  //interactorStyle->setAllAnglesImages(allImage);









  interactor->SetInteractorStyle(interactorStyle);
  interactor->Start();
  if(initialObliqueData != 0 && initialObliqueData->GetReferenceCount() == 1)
  {
	  initialObliqueData->Delete();
  }
  return EXIT_SUCCESS;
}

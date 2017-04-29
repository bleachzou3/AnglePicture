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
	AnglePictureUtility::segment("E:\\patientData\\WU_AMAO","E:\\patientData\\WU_AMAO_SEGMENTED",307,166,121,300,1000);
	//AnglePictureUtility::coronaryVoxelRender("E:\\patientData\\WU_AMAO_SEGMENTED");
	//AnglePictureUtility::showHistogram("E:\\patientData\\WU_AMAO");
  return EXIT_SUCCESS;
}

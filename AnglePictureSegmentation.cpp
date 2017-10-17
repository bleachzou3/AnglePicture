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
#include <DisplayVoxelUtility.h>
#include "MeanShiftAlgo.h"
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
	//AnglePictureUtility::segment("E:\\patientData\\WU_AMAO","E:\\patientData\\WU_AMAO_SEGMENTED",307,166,121,300,1000);
	//AnglePictureUtility::coronaryVoxelRender("E:\\patientData\\WU_AMAO_SEGMENTED");
	//AnglePictureUtility::showHistogram("E:\\patientData\\WU_AMAO");
	//AnglePictureUtility::SegmentBloodVessels("E:\\patientData\\WU_AMAO","E:\\patientData\\WU_AMAO_SegmentBloodVessels",0,0,0);
	//AnglePictureUtility::WatershedSegmentation("E:\\patientData\\WU_AMAO_image_volume_voi.vti","E:\\patientData\\WU_AMAO_image_volume_voi_WatershedSegmentation");
	//#AnglePictureUtility::SegmentBloodVesselsFromVti("E:\\patientData\\WU_AMAO_image_volume_voi.vti","E:\\patientData\\hessian_WU_AMAO_image_volume_2.vti",0,0,0);
	//DisplayVoxelUtility::displaySegmentBloodVesselsFromVti("E:\\patientData\\WU_AMAO_image_volume.vti");
	//AnglePictureUtility::CurvesLevelSetImage("E:\\patientData\\WU_AMAO_image_volume_voi.vti",120,111,182);
	
	//AnglePictureUtility::SegmentBloodVesselsWithMultiScaleHessianBasedMeasure("E:\\patientData\\WU_AMAO_image_volume_voi.vti","E:\\patientData\\WU_AMAO_SegmentBloodVesselsWithMultiScaleHessianBasedMeasure.vti",1,10.0,10);


	//公司在跑实验室的数据的时候，一定会有某个参考标准，而我这里的确是没有，我现在只想跑出来，看看效果
	AnglePictureUtility::segmentMeanShiftClustering("E:\\patientData\\WU_AMAO_image_volume_voi.vti","E:\\patientData\\WU_AMAO_segmentMeanShiftClustering.vti",10,20,10,30,1000);

    return EXIT_SUCCESS;
}

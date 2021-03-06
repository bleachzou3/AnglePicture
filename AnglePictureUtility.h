#ifndef ANGLE_PICTURE_UTILITY_H_
#define ANGLE_PICTURE_UTILITY_H_
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include "Vector3.h"
#include <vtkSmartPointer.h>
#include <unordered_map>
class AnglePictureUtility
{
public:
	AnglePictureUtility();
	/**
	*@param original 原始的图像切片, direction 中心线的方向, fixedPoint图像切片中心经过的点
	*
	*/
	static void computeOblique(vtkImageData* original,Vector3 direction,Vector3 fixedPoint,vtkImageData*res);

	/**
	*计算每个切面的法向量obliqueNormal
	* 计算左右相邻的两个点，哪个点近就用哪个
	*/
	static void computeNormalByPoints(vtkPolyData * data);

	/**
	*images里面vtkImageData的销毁交给调用者
	*/
	static void computeAllAngleImages(vtkPolyData*data,unordered_map<vtkIdType,vtkImageData*>& images,vtkImageData* originalImage);

	/**
	*用itk bfs切割图像
	*/
	static bool segment(string directoryName,string outputFileName,int x,int y, int z,float lowerThreshold,float uppperThreshold);

	/**
	*选择存放已经切割好的血管的文件夹，显示血管
	*/
	static void coronaryVoxelRender(string directoryName);

	/**
	*显示histogram
	*/
	static void showHistogram(string directoryName);

	/**
	*SegmentBloodVesselsWithMultiScaleHessianBasedMeasure这个方法
	*/
	static void SegmentBloodVesselsWithMultiScaleHessianBasedMeasure(string inputDirectoryName,string outputDirectoryName);

	/**
	*SegmentBloodVessels这个方法
	*/
	static void SegmentBloodVessels(string inputDirectoryName,string outputDirectoryName,double sigma,double alpha1,double alpha2);

	/**
	*试一下water shed效果
	*/
	static void WatershedSegmentation(string inputDirectoryName,string outputDirectoryName);

	/**
	*从vti读取图像,进行图像切割
	*/
	static void SegmentBloodVesselsFromVti(string fileName,string outputFileName,double sigma,double alpha1,double alpha2);


	/**
	*从文件里面读取图像
	*/
	static void CurvesLevelSetImage(string inputFileName,int x,int y,int z);
	static const string  OBLIQUE_NORMAL;
	static const string  MAXIMUM_INSCRIBED_SPHERE_RADIUS;  

	//这个是多尺度
	static void SegmentBloodVesselsWithMultiScaleHessianBasedMeasure(string inputFileName,string outputFileName,double sigmaMinimum,double sigmaMaximum,unsigned int numberOfSigmaSteps);


	/**
	*mean-shift用于图像切割
	*/
	static void segmentMeanShiftClustering(string inputFileName,string outFileName,int _radius,float _grayDistance, int _numOfCenters,int _itertionNums,int _clusterPointSize);
};
#endif
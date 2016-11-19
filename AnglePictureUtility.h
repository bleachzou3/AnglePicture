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


	static const string  OBLIQUE_NORMAL;
	static const string  MAXIMUM_INSCRIBED_SPHERE_RADIUS;                          
};
#endif
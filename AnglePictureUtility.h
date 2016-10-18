#ifndef ANGLE_PICTURE_UTILITY_H_
#define ANGLE_PICTURE_UTILITY_H_
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include "Vector3.h"
#include <vtkSmartPointer.h>
class AnglePictureUtility
{
public:
	AnglePictureUtility();
	/**
	*@param original ԭʼ��ͼ����Ƭ, direction �����ߵķ���, fixedPointͼ����Ƭ���ľ����ĵ�
	*
	*/
	static void computeOblique(vtkImageData* original,Vector3 direction,Vector3 fixedPoint,vtkImageData*res);

	/**
	*����ÿ������ķ�����obliqueNormal
	* �����������ڵ������㣬�ĸ���������ĸ�
	*/
	static void computeNormalByPoints(vtkPolyData * data);
};
#endif
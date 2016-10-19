#include "AnglePictureUtility.h"
#include <vtkImageData.h>
#include "Vector3.h"
#include <vtkSmartPointer.h>
#include <vtkMatrix4x4.h>
#include <vtkImageReslice.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkDataObject.h>
#include "ONB.h"

AnglePictureUtility::AnglePictureUtility()
{
	
};
const string AnglePictureUtility::OBLIQUE_NORMAL =  "obliqueNormal";
const string AnglePictureUtility::MAXIMUM_INSCRIBED_SPHERE_RADIUS = "MaximumInscribedSphereRadius";
void AnglePictureUtility::computeOblique(vtkImageData* original,Vector3 direction,Vector3 fixedPoint,vtkImageData*res)
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

   double axialElements[16] = {
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


void AnglePictureUtility::computeNormalByPoints(vtkPolyData * data)
{
	vtkDataArray* normals = data->GetPointData()->GetArray("obliqueNormal");

	if(normals != 0)
	{
		return;
	}
	int N = data->GetNumberOfPoints();
	int line = data->GetNumberOfLines();
	if(N < 3)
	{
		return;
	}
	vtkSmartPointer<vtkFloatArray> newNormals = vtkSmartPointer<vtkFloatArray>::New();
	newNormals->SetName("obliqueNormal");

	newNormals->SetNumberOfComponents(3);
	newNormals->SetNumberOfTuples(N);

	double firstNode[3]; 
	data->GetPoint(0,firstNode);
	double secondNode[3];
	data->GetPoint(1,secondNode);
	newNormals->SetComponent(0,0,firstNode[0] - secondNode[0]);
	newNormals->SetComponent(0,1,firstNode[1] - secondNode[1]);
	newNormals->SetComponent(0,2,firstNode[2] - secondNode[2]);

	for(int i = 1; i < N - 1; i++)
	{
		double currentNode[3];
		data->GetPoint(i,currentNode);
		double previousNode[3];
		double nextNode[3] ;
		data->GetPoint(i - 1,previousNode);
		data->GetPoint(i + 1,nextNode);

		double prevCurr = sqrt(vtkMath::Distance2BetweenPoints(currentNode,previousNode));
		double nextCurr = sqrt(vtkMath::Distance2BetweenPoints(currentNode,nextNode));

		double curNormal[3];
		if(prevCurr < nextCurr)
		{
			curNormal[0] = previousNode[0] - currentNode[0];
			curNormal[1] = previousNode[1] - currentNode[1];
			curNormal[2] = previousNode[2] - currentNode[2];
		}else
		{
			curNormal[0] = currentNode[0] - nextNode[0];
			curNormal[1] = currentNode[1] - nextNode[1];
			curNormal[2] = currentNode[2] - nextNode[2];
		}
		newNormals->SetComponent(i,0,curNormal[0]);
		newNormals->SetComponent(i,1,curNormal[1]);
	    newNormals->SetComponent(i,2,curNormal[2]);
	}

	double* rFirstNode = data->GetPoint(N - 1);
	double* rSecondNode = data->GetPoint(N - 2);

    newNormals->SetComponent(N - 1,0,firstNode[0] - secondNode[0]);
	newNormals->SetComponent(N - 1,1,firstNode[1] - secondNode[1]);
	newNormals->SetComponent(N - 1,2,firstNode[2] - secondNode[2]);

	data->GetPointData()->AddArray(newNormals);

}

void AnglePictureUtility::computeAllAngleImages(vtkPolyData*data,unordered_map<vtkIdType,vtkImageData*>&images,vtkImageData* originalImage)
{
	int N = data->GetNumberOfPoints();
	vtkDataArray* normals = data->GetPointData()->GetArray(OBLIQUE_NORMAL.c_str());
	for(int i = 0; i < N; i++)
	{
		double center[3];
		data->GetPoint(i,center);
		Vector3 fixedPoint(center[0],center[1],center[2]);

		double* normal = normals->GetTuple(i);
		Vector3 direction(normal[0],normal[1],normal[2]);
		vtkImageData* temp = vtkImageData::New();
		computeOblique(originalImage,direction,fixedPoint,temp);
		images[i] = temp;
	}

}
#include "DisplayVoxelUtility.h"
#include <vtkXMLImageDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkImageCast.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkImageData.h>
#include <vtkVolume.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
DisplayVoxelUtility::DisplayVoxelUtility()
{
}
DisplayVoxelUtility::~DisplayVoxelUtility()
{
}

void DisplayVoxelUtility::displaySegmentBloodVesselsFromVti(string fileName)
{

	 vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	 reader->SetFileName(fileName.c_str());
	 reader->Update();



   vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
   cast->SetInputData(reader->GetOutput());
   cast->SetOutputScalarTypeToUnsignedShort();
   cast->Update();






   	vtkSmartPointer<vtkVolumeRayCastCompositeFunction> rayCastFun = 
		vtkSmartPointer<vtkVolumeRayCastCompositeFunction>::New();

	vtkSmartPointer<vtkVolumeRayCastMapper> volumeMapper = 
		vtkSmartPointer<vtkVolumeRayCastMapper>::New();
	volumeMapper->SetInputData(cast->GetOutput());
	volumeMapper->SetVolumeRayCastFunction(rayCastFun);

	//���ù��߲�������
	//volumeMapper->SetSampleDistance(volumeMapper->GetSampleDistance()*4);
	//����ͼ���������
	//volumeMapper->SetAutoAdjustSampleDistances(0);
	//volumeMapper->SetImageSampleDistance(4);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty = 
		vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetInterpolationTypeToLinear();
	//volumeProperty->ShadeOn();  //�򿪻��߹ر���Ӱ����
	volumeProperty->SetAmbient(0.4);
	volumeProperty->SetDiffuse(0.6);
	volumeProperty->SetSpecular(0.2);

	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = 
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	//compositeOpacity->AddPoint(1998,   0.00);
	compositeOpacity->AddPoint(100,0.0);
	compositeOpacity->AddPoint(255,   0.70);
	volumeProperty->SetScalarOpacity(compositeOpacity); //���ò�͸���ȴ��亯��

	//�������ز�������,�ԱȲ�ͬ������
	//compositeOpacity->AddPoint(120,  0.00);
	//compositeOpacity->AddPoint(180,  0.60);
	//volumeProperty->SetScalarOpacity(compositeOpacity);

	
	vtkSmartPointer<vtkPiecewiseFunction> volumeGradientOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	volumeGradientOpacity->AddPoint(10,  0.0);
	volumeGradientOpacity->AddPoint(90,  0.5);
	volumeGradientOpacity->AddPoint(100, 1.0);
	//volumeProperty->SetGradientOpacity(volumeGradientOpacity);//�����ݶȲ�͸����Ч���Ա�

	vtkSmartPointer<vtkColorTransferFunction> color = 
		vtkSmartPointer<vtkColorTransferFunction>::New();
	color->AddRGBPoint(0.000,  0.00, 0.00, 0.00);
	color->AddRGBPoint(255,  1.00, 1, 1);
	//color->AddRGBPoint(190.0,  1.00, 1.00, 1.00);
	//color->AddRGBPoint(220.0,  0.20, 0.20, 0.20);
	volumeProperty->SetColor(color);

	vtkSmartPointer<vtkVolume> volume = 
		vtkSmartPointer<vtkVolume>::New();
	volume->SetMapper(volumeMapper);
	volume->SetProperty(volumeProperty);

	vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
	ren->SetBackground(1.0, 1.0, 1.0);
	ren->AddVolume( volume ); 

	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren);
	renWin->SetSize(640, 480);
	renWin->Render();
	renWin->SetWindowName("VolumeRenderingApp");

	vtkSmartPointer<vtkRenderWindowInteractor> iren = 
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);
	ren->ResetCamera();

	renWin->Render();
	iren->Start();

}
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

#include <itkImage.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include <itkImageSeriesReader.h>
#include <itkImageFileWriter.h>
#include <itkConnectedThresholdImageFilter.h>
#include <itkImageSeriesWriter.h>
#include <vtkDICOMImageReader.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVolumeProperty.h>
#include <vtkAxesActor.h>
#include <vtkImageShiftScale.h>
#include <vtkImageCast.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkImageAccumulate.h>
#include <vtkBarChartActor.h>
#include <vtkActor.h>
#include <vtkBarChartActor.h>
#include <vtkFieldData.h>
#include <vtkImageAccumulate.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkJPEGReader.h>
#include <vtkLegendBoxActor.h>
#include <vtkProperty2D.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTextProperty.h>

#include <itkHessian3DToVesselnessMeasureImageFilter.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkWatershedImageFilter.h>
#include <itkScalarToRGBPixelFunctor.h>
#include <itkUnaryFunctorImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <vtkXMLImageDataReader.h>
#include <itkVTKImageToImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <vtkXMLImageDataWriter.h>


#include <itkCurvatureAnisotropicDiffusionImageFilter.h>
#include <itkGradientMagnitudeRecursiveGaussianImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkFastMarchingImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkCurvesLevelSetImageFilter.h>
#include <log4cpp/Category.hh>
#include <log4cpp/PropertyConfigurator.hh>
extern log4cpp::Category& rootLog;
extern log4cpp::Category& subLog ;

using namespace itk;
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

bool AnglePictureUtility::segment(string directoryName,string outputDirectory,int x,int y,int z,float lowerThreshold,float uppperThreshold)
{
	// Software Guide : BeginLatex
//
// We define the pixel type and dimension of the image to be read. In this
// particular case, the dimensionality of the image is 3, and we assume a
// \code{signed short} pixel type that is commonly used for X-Rays CT scanners.
//
// The image orientation information contained in the direction cosines
// of the DICOM header are read in and passed correctly down the image processing
// pipeline.
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  typedef itk::Image< PixelType, Dimension >         ImageType;
// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
// We use the image type for instantiating the type of the series reader and
// for constructing one object of its type.
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
// A GDCMImageIO object is created and connected to the reader. This object is
// the one that is aware of the internal intricacies of the DICOM format.
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  reader->SetImageIO( dicomIO );
// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
// Now we face one of the main challenges of the process of reading a DICOM
// series: to identify from a given directory the set of filenames
// that belong together to the same volumetric image. Fortunately for us, GDCM
// offers functionalities for solving this problem and we just need to invoke
// those functionalities through an ITK class that encapsulates a communication
// with GDCM classes. This ITK object is the GDCMSeriesFileNames. Conveniently,
// we only need to pass to this class the name of the directory where
// the DICOM slices are stored. This is done with the \code{SetDirectory()}
// method. The GDCMSeriesFileNames object will explore the directory and will
// generate a sequence of filenames for DICOM files for one study/series.
// In this example, we also call the \code{SetUseSeriesDetails(true)} function
// that tells the GDCMSeriesFileNames object to use additional DICOM
// information to distinguish unique volumes within the directory.  This is
// useful, for example, if a DICOM device assigns the same SeriesID to
// a scout scan and its 3D volume; by using additional DICOM information
// the scout scan will not be included as part of the 3D volume.  Note that
// \code{SetUseSeriesDetails(true)} must be called prior to calling
// \code{SetDirectory()}. By default \code{SetUseSeriesDetails(true)} will use
// the following DICOM tags to sub-refine a set of files into multiple series:
//
// \begin{description}
// \item[0020 0011] Series Number
// \item[0018 0024] Sequence Name
// \item[0018 0050] Slice Thickness
// \item[0028 0010] Rows
// \item[0028 0011] Columns
// \end{description}
//
// If this is not enough for your specific case you can always add some more
// restrictions using the \code{AddSeriesRestriction()} method. In this example we will use
// the DICOM Tag: 0008 0021 DA 1 Series Date, to sub-refine each series. The format
// for passing the argument is a string containing first the group then the element
// of the DICOM tag, separated by a pipe ($|$) sign.
//
//
// \index{itk::GDCMSeriesFileNames!SetDirectory()}
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetDirectory( directoryName );
// Software Guide : EndCodeSnippet
  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
	std::cout << std::endl << directoryName << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;
// Software Guide : BeginLatex
//
// The GDCMSeriesFileNames object first identifies the list of DICOM series
// present in the given directory. We receive that list in a reference
// to a container of strings and then we can do things like print out all
// the series identifiers that the generator had found. Since the process of
// finding the series identifiers can potentially throw exceptions, it is
// wise to put this code inside a \code{try/catch} block.
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
    typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
      }
// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
// Given that it is common to find multiple DICOM series in the same directory,
// we must tell the GDCM classes what specific series we want to read. In
// this example we do this by checking first if the user has provided a series
// identifier in the command line arguments. If no series identifier has been
// passed, then we simply use the first series found during the exploration of
// the directory.
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
    std::string seriesIdentifier;
  //  if( argc > 3 ) // If no optional series identifier
      {
  //    seriesIdentifier = argv[3];
      }
    //else
    //  {
      seriesIdentifier = seriesUID.begin()->c_str();
     // }
// Software Guide : EndCodeSnippet
    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;
// Software Guide : BeginLatex
//
// We pass the series identifier to the name generator and ask for all the
// filenames associated to that series. This list is returned in a container of
// strings by the \code{GetFileNames()} method.
//
// \index{itk::GDCMSeriesFileNames!GetFileNames()}
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;
    fileNames = nameGenerator->GetFileNames( seriesIdentifier );
// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
//
// The list of filenames can now be passed to the \doxygen{ImageSeriesReader}
// using the \code{SetFileNames()} method.
//
//  \index{itk::ImageSeriesReader!SetFileNames()}
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
    reader->SetFileNames( fileNames );
// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
// Finally we can trigger the reading process by invoking the \code{Update()}
// method in the reader. This call as usual is placed inside a \code{try/catch}
// block.
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
      }

	ImageType* image = reader->GetOutput();
	//解决图像的大小

	ImageType::RegionType region = image->GetLargestPossibleRegion();  
    ImageType::SizeType size  = region.GetSize();
	cout << size[0] << " " << size[1] << " " << size[2] << endl;
	cout << typeid(ImageType::PixelType).raw_name() << endl;
	signed short high = -20000;
	signed short low = 20000;
	ImageType::IndexType pixelIndex;
	
	for(int z = 0; z < size[2] ; z++)
	{
		for(int x = 0; x < size[0]; x++)
		{
			for(int y = 0; y < size[1]; y++)
			{
				pixelIndex[0] = x;
				pixelIndex[1] = y;
				pixelIndex[2] = z;
				if(image->GetPixel(pixelIndex) == 1252)
				{
					cout <<"position:  " << x << "  " << y <<" " <<  z << endl;
				}
				high = std::max(high,image->GetPixel(pixelIndex));
				low = std::min(low,image->GetPixel(pixelIndex));
			}
		}
	}
	cout << "hight" << high << "   low:" << low << endl;
	

	
	pixelIndex[0] = x;
	pixelIndex[1] = y;
	pixelIndex[2] = z;
	
	

	
	//打印种子点的灰度值
	cout << "当前种子点的灰度值 pixelIndex" << image->GetPixel(pixelIndex) << endl;

	//今天下午
	typedef itk::ConnectedThresholdImageFilter< ImageType,
		ImageType > ConnectedFilterType;

	 ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
	 connectedThreshold->SetInput(image);
	 connectedThreshold->SetSeed(pixelIndex);
	 connectedThreshold->SetLower(lowerThreshold);
	 connectedThreshold->SetUpper(uppperThreshold);
	 connectedThreshold->SetReplaceValue(2000);
	 connectedThreshold->Update();







	  itksys::SystemTools::MakeDirectory( outputDirectory );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;
  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
	 typedef itk::ImageSeriesWriter<
		 ImageType,Image2DType>  SeriesWriterType;
  // Software Guide : EndCodeSnippet
  //  Software Guide : BeginLatex
  //
  //  We construct a series writer and connect to its input the output from the
  //  reader. Then we pass the GDCM image IO object in order to be able to write
  //  the images in DICOM format.
  //
  //  the writer filter.  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput(connectedThreshold->GetOutput() );
  seriesWriter->SetImageIO( dicomIO );
  // Software Guide : EndCodeSnippet
  //  Software Guide : BeginLatex
  //
  //  It is time now to setup the GDCMSeriesFileNames to generate new filenames
  //  using another output directory.  Then simply pass those newly generated
  //  files to the series writer.
  //
  //  \index{GDCMSeriesFileNames!SetOutputDirectory()}
  //  \index{GDCMSeriesFileNames!GetOutputFileNames()}
  //  \index{ImageSeriesWriter!SetFileNames()}
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  nameGenerator->SetOutputDirectory( outputDirectory );
  seriesWriter->SetFileNames( nameGenerator->GetOutputFileNames() );
  // Software Guide : EndCodeSnippet
  //  Software Guide : BeginLatex
  //
  //  The following line of code is extremely important for this process to work
  //  correctly.  The line is taking the MetaDataDictionary from the input reader
  //  and passing it to the output writer. This step is important because the
  //  MetaDataDictionary contains all the entries of the input DICOM header.
  //
  //  \index{itk::ImageSeriesReader!GetMetaDataDictionaryArray()}
  //  \index{itk::ImageSeriesWriter!SetMetaDataDictionaryArray()}
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  seriesWriter->SetMetaDataDictionaryArray(
                        reader->GetMetaDataDictionaryArray() );
  // Software Guide : EndCodeSnippet
  // Software Guide : BeginLatex
  //
  // Finally we trigger the writing process by invoking the \code{Update()} method
  // in the series writer. We place this call inside a \code{try/catch} block,
  // in case any exception is thrown during the writing process.
  //
  // Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  try
    {
    seriesWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }
















    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
    }
  // Software Guide : BeginLatex
  //
  // Note that in addition to writing the volumetric image to a file we could
  // have used it as the input for any 3D processing pipeline. Keep in mind that
  // DICOM is simply a file format and a network protocol. Once the image data
  // has been loaded into memory, it behaves as any other volumetric dataset that
  // you could have loaded from any other file format.
  //
  // Software Guide : EndLatex

}

void AnglePictureUtility::coronaryVoxelRender(string directoryName)
{
	   // Read all the DICOM files in the specified directory.
   vtkSmartPointer<vtkDICOMImageReader> reader =
      vtkSmartPointer<vtkDICOMImageReader>::New();
   reader->SetDirectoryName(directoryName.c_str());
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

	//设置光线采样距离
	//volumeMapper->SetSampleDistance(volumeMapper->GetSampleDistance()*4);
	//设置图像采样步长
	//volumeMapper->SetAutoAdjustSampleDistances(0);
	//volumeMapper->SetImageSampleDistance(4);

	vtkSmartPointer<vtkVolumeProperty> volumeProperty = 
		vtkSmartPointer<vtkVolumeProperty>::New();
	volumeProperty->SetInterpolationTypeToLinear();
	//volumeProperty->ShadeOn();  //打开或者关闭阴影测试
	volumeProperty->SetAmbient(0.4);
	volumeProperty->SetDiffuse(0.6);
	volumeProperty->SetSpecular(0.2);

	vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = 
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	compositeOpacity->AddPoint(1998,   0.00);
	compositeOpacity->AddPoint(2000,   0.40);
	volumeProperty->SetScalarOpacity(compositeOpacity); //设置不透明度传输函数

	//测试隐藏部分数据,对比不同的设置
	//compositeOpacity->AddPoint(120,  0.00);
	//compositeOpacity->AddPoint(180,  0.60);
	//volumeProperty->SetScalarOpacity(compositeOpacity);

	
	vtkSmartPointer<vtkPiecewiseFunction> volumeGradientOpacity =
		vtkSmartPointer<vtkPiecewiseFunction>::New();
	volumeGradientOpacity->AddPoint(10,  0.0);
	volumeGradientOpacity->AddPoint(90,  0.5);
	volumeGradientOpacity->AddPoint(100, 1.0);
	//volumeProperty->SetGradientOpacity(volumeGradientOpacity);//设置梯度不透明度效果对比

	vtkSmartPointer<vtkColorTransferFunction> color = 
		vtkSmartPointer<vtkColorTransferFunction>::New();
	color->AddRGBPoint(0.000,  0.00, 0.00, 0.00);
	color->AddRGBPoint(1999,  1.00, 0.52, 0.30);
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

void AnglePictureUtility::showHistogram(string directoryName)
{
   vtkSmartPointer<vtkDICOMImageReader> reader =
      vtkSmartPointer<vtkDICOMImageReader>::New();
   reader->SetDirectoryName(directoryName.c_str());
   reader->Update();

   vtkImageData * image = reader->GetOutput();

   //把灰度值的范围找出来,遍历访问就行了
    int dims[3];
	image->GetDimensions(dims);
	signed short top = -3000;
	signed short down = 4000;
    for(int k = 0; k < dims[2]; k++)
	{
		for(int j = 0; j < dims[1]; j++)
		{
			for(int i = 0; i < dims[0]; i++)
			{
				signed short* pixel = (signed  short*)image->GetScalarPointer(i,j,k);
				//cout << pixel[0] << endl;
				top = std::max(pixel[0],top);
				down = std::min(pixel[0],down);
			}
		}
	}
	cout << "top: " << top << "  down: " << down << endl;
	int bins   = 10;
	int comps  = 1;

	vtkSmartPointer<vtkImageAccumulate> histogram =
		vtkSmartPointer<vtkImageAccumulate>::New();
	histogram->SetInputData(image);
	histogram->SetComponentExtent(0, bins-1, 0, 0, 0, 0);
	histogram->SetComponentOrigin(down, 0, 0);
	histogram->SetComponentSpacing((top - down)/bins, 0, 0);
	histogram->Update();

	int* output = static_cast<int*>(histogram->GetOutput()->GetScalarPointer());

	vtkSmartPointer<vtkIntArray> frequencies = 
		vtkSmartPointer<vtkIntArray>::New();
	frequencies->SetNumberOfComponents(1);

	for(int j = 0; j < bins; ++j)
	{
		for(int i=0; i<comps; i++)
		{
			cout << *output << endl;
			frequencies->InsertNextTuple1(*output++);
		}
	}

	vtkSmartPointer<vtkDataObject> dataObject = 
		vtkSmartPointer<vtkDataObject>::New();
	dataObject->GetFieldData()->AddArray( frequencies );

	vtkSmartPointer<vtkBarChartActor> barChart = 
		vtkSmartPointer<vtkBarChartActor>::New();
	barChart->SetInput(dataObject);
	barChart->SetTitle("Histogram");
	barChart->GetPositionCoordinate()->SetValue(0.05,0.05,0.0);
	barChart->GetPosition2Coordinate()->SetValue(0.95,0.95,0.0);
	barChart->GetProperty()->SetColor(0,0,0);
	barChart->GetTitleTextProperty()->SetColor(0,0,0);
	barChart->GetLabelTextProperty()->SetColor(0,0,0);
	barChart->GetLegendActor()->SetNumberOfEntries(dataObject->GetFieldData()->GetArray(0)->GetNumberOfTuples());
	barChart->LegendVisibilityOn();
	barChart->LabelVisibilityOn();

	double colors[3][3] = {
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 } };

	int count = 0;
	for( int i = 0; i < bins; ++i )
	{
		for( int j = 0; j < comps; ++j )
		{
			barChart->SetBarColor( count++, colors[j] );
		}
	}

	vtkSmartPointer<vtkRenderer> renderer = 
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(barChart);
	renderer->SetBackground(1.0, 1.0, 1.0);

	vtkSmartPointer<vtkRenderWindow> renderWindow = 
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize(640, 480);
	renderWindow->Render();
	renderWindow->SetWindowName("ImageAccumulateExample");

	vtkSmartPointer<vtkRenderWindowInteractor> interactor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	interactor->Initialize();
	interactor->Start();

	
}

void AnglePictureUtility::SegmentBloodVesselsWithMultiScaleHessianBasedMeasure(string inputDirectoryName,string outputDirectoryName)
{
	// Software Guide : BeginCodeSnippet
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  typedef itk::Image< PixelType, Dimension >         ImageType;

  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  reader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetDirectory( inputDirectoryName );
// Software Guide : EndCodeSnippet
  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
	std::cout << std::endl << inputDirectoryName << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;

    typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
      }

    std::string seriesIdentifier;

      seriesIdentifier = seriesUID.begin()->c_str();

    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;

    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;
    fileNames = nameGenerator->GetFileNames( seriesIdentifier );

    reader->SetFileNames( fileNames );

    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
 
      }

	ImageType* image = reader->GetOutput();
	























































	itksys::SystemTools::MakeDirectory( outputDirectoryName );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;
  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
	 typedef itk::ImageSeriesWriter<
		 ImageType,Image2DType>  SeriesWriterType;
  // Software Guide : EndCodeSnippet
  //  Software Guide : BeginLatex
  //
  //  We construct a series writer and connect to its input the output from the
  //  reader. Then we pass the GDCM image IO object in order to be able to write
  //  the images in DICOM format.
  //
  //  the writer filter.  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  //seriesWriter->SetInput(connectedThreshold->GetOutput() );
  seriesWriter->SetImageIO( dicomIO );
  // Software Guide : EndCodeSnippet
  //  Software Guide : BeginLatex
  //
  //  It is time now to setup the GDCMSeriesFileNames to generate new filenames
  //  using another output directory.  Then simply pass those newly generated
  //  files to the series writer.
  //
  //  \index{GDCMSeriesFileNames!SetOutputDirectory()}
  //  \index{GDCMSeriesFileNames!GetOutputFileNames()}
  //  \index{ImageSeriesWriter!SetFileNames()}
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  nameGenerator->SetOutputDirectory( outputDirectoryName );
  seriesWriter->SetFileNames( nameGenerator->GetOutputFileNames() );
  // Software Guide : EndCodeSnippet
  //  Software Guide : BeginLatex
  //
  //  The following line of code is extremely important for this process to work
  //  correctly.  The line is taking the MetaDataDictionary from the input reader
  //  and passing it to the output writer. This step is important because the
  //  MetaDataDictionary contains all the entries of the input DICOM header.
  //
  //  \index{itk::ImageSeriesReader!GetMetaDataDictionaryArray()}
  //  \index{itk::ImageSeriesWriter!SetMetaDataDictionaryArray()}
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  seriesWriter->SetMetaDataDictionaryArray(
                        reader->GetMetaDataDictionaryArray() );
  // Software Guide : EndCodeSnippet
  // Software Guide : BeginLatex
  //
  // Finally we trigger the writing process by invoking the \code{Update()} method
  // in the series writer. We place this call inside a \code{try/catch} block,
  // in case any exception is thrown during the writing process.
  //
  // Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  try
    {
    seriesWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    return;
    }
















    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return ;
    }
  // Software Guide : BeginLatex
  //
  // Note that in addition to writing the volumetric image to a file we could
  // have used it as the input for any 3D processing pipeline. Keep in mind that
  // DICOM is simply a file format and a network protocol. Once the image data
  // has been loaded into memory, it behaves as any other volumetric dataset that
  // you could have loaded from any other file format.
  //
  // Software Guide : EndLatex

}

void AnglePictureUtility::SegmentBloodVessels(string inputDirectoryName,string outputDirectoryName,double sigma,double alpha1,double alpha2)
{
	// Software Guide : BeginCodeSnippet
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  typedef itk::Image< PixelType, Dimension >         ImageType;

  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  reader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetDirectory( inputDirectoryName );
// Software Guide : EndCodeSnippet
  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
	std::cout << std::endl << inputDirectoryName << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;

    typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
      }

    std::string seriesIdentifier;

      seriesIdentifier = seriesUID.begin()->c_str();

    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;

    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;
    fileNames = nameGenerator->GetFileNames( seriesIdentifier );

    reader->SetFileNames( fileNames );

    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
 
      }

	ImageType* image = reader->GetOutput();
	















	 cout << "2222222222222222222222222222222222222222222" << endl;
	typedef itk::HessianRecursiveGaussianImageFilter< ImageType >
    HessianFilterType;
  HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
  hessianFilter->SetInput( image );
  hessianFilter->Update();
  cout << "1111111111111111111111111111111111111" << endl;

  if( sigma )
    {
    hessianFilter->SetSigma( sigma  );
    }
  typedef itk::Hessian3DToVesselnessMeasureImageFilter< PixelType >
    VesselnessMeasureFilterType;
  VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();
  vesselnessFilter->SetInput( hessianFilter->GetOutput() );
  if( alpha1 )
    {
    vesselnessFilter->SetAlpha1( alpha1  );
    }
  if( alpha2 )
    {
    vesselnessFilter->SetAlpha2(  alpha2 );
    }

























  cout << "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh"  << endl;
























	itksys::SystemTools::MakeDirectory( outputDirectoryName );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;
  typedef itk::Image< OutputPixelType, OutputDimension >    Image2DType;
	 typedef itk::ImageSeriesWriter<
		 ImageType,Image2DType>  SeriesWriterType;

  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput(vesselnessFilter->GetOutput());
  seriesWriter->SetImageIO( dicomIO );

  nameGenerator->SetOutputDirectory( outputDirectoryName );
  seriesWriter->SetFileNames( nameGenerator->GetOutputFileNames() );

  seriesWriter->SetMetaDataDictionaryArray(
                        reader->GetMetaDataDictionaryArray() );

  try
    {
    seriesWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    return;
    }
















    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return ;
    }

}

void   AnglePictureUtility::WatershedSegmentation(string inputDirectoryName,string outputDirectoryName)
{
	
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  typedef itk::Image< PixelType, Dimension >         ImageType;

  typedef itk::ImageSeriesReader< ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

  typedef itk::GDCMImageIO       ImageIOType;
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  reader->SetImageIO( dicomIO );

  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  nameGenerator->SetUseSeriesDetails( true );
  nameGenerator->AddSeriesRestriction("0008|0021" );
  nameGenerator->SetDirectory( inputDirectoryName );
// Software Guide : EndCodeSnippet
  try
    {
    std::cout << std::endl << "The directory: " << std::endl;
	std::cout << std::endl << inputDirectoryName << std::endl << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;

    typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
      }

    std::string seriesIdentifier;
  //  if( argc > 3 ) // If no optional series identifier
      {
  //    seriesIdentifier = argv[3];
      }
    //else
    //  {
      seriesIdentifier = seriesUID.begin()->c_str();
     // }
// Software Guide : EndCodeSnippet
    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl;
    std::cout << std::endl << std::endl;

    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames;
    fileNames = nameGenerator->GetFileNames( seriesIdentifier );

    reader->SetFileNames( fileNames );

    try
      {
      reader->Update();
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return ;
      }

	ImageType* image = reader->GetOutput();
	//解决图像的大小
 typedef itk::GradientMagnitudeImageFilter<
	 ImageType, ImageType >  FilterType;
  FilterType::Pointer gradientFilter = FilterType::New();
  gradientFilter->SetInput( image );
  gradientFilter->Update();
  cout << "11111111111111111111111111111" << endl;

	 typedef  itk::WatershedImageFilter<
		 ImageType
                                            > WatershedFilterType;
  WatershedFilterType::Pointer watershedFilter = WatershedFilterType::New();
  watershedFilter->SetInput( gradientFilter->GetOutput() );
  watershedFilter->SetThreshold( 300 );
  watershedFilter->SetLevel(     500 );
   watershedFilter->Update();
    cout << "11111111111111111111111111111111111111" << endl;
  //  Instantiate the filter that will encode the label image
  //  into a color image (random color attribution).
  //
  typedef itk::Functor::ScalarToRGBPixelFunctor<
	  unsigned long
                                                    > ColorMapFunctorType;
  typedef WatershedFilterType::OutputImageType  LabeledImageType;
  typedef itk::RGBPixel<unsigned short>      RGBPixelType;
  typedef itk::Image< RGBPixelType,       Dimension >  RGBImageType;
  typedef itk::UnaryFunctorImageFilter<
                                LabeledImageType,
                                RGBImageType,
                                ColorMapFunctorType
                                                > ColorMapFilterType;
  ColorMapFilterType::Pointer colorMapFilter = ColorMapFilterType::New();
 
  colorMapFilter->SetInput(  watershedFilter->GetOutput() );
  colorMapFilter->Update();
  cout << "222222222222222222222222222222222222222222222222222222222222222222" << endl;










































	itksys::SystemTools::MakeDirectory( outputDirectoryName );
  typedef signed short    OutputPixelType;
  const unsigned int      OutputDimension = 2;
  typedef itk::Image< RGBPixelType, OutputDimension >    Image2DType;
	 typedef itk::ImageSeriesWriter<
		 RGBImageType,Image2DType>  SeriesWriterType;

  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  seriesWriter->SetInput(colorMapFilter->GetOutput());
  seriesWriter->SetImageIO( dicomIO );

  nameGenerator->SetOutputDirectory( outputDirectoryName );
  seriesWriter->SetFileNames( nameGenerator->GetOutputFileNames() );

  seriesWriter->SetMetaDataDictionaryArray(
                        reader->GetMetaDataDictionaryArray() );

  try
    {
    seriesWriter->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the series " << std::endl;
    std::cerr << excp << std::endl;
    return ;
    }
















    }
  catch (itk::ExceptionObject &ex)
    {
    std::cout << ex << std::endl;
    return ;
    }
  // Software Guide : BeginLatex
  //
  // Note that in addition to writing the volumetric image to a file we could
  // have used it as the input for any 3D processing pipeline. Keep in mind that
  // DICOM is simply a file format and a network protocol. Once the image data
  // has been loaded into memory, it behaves as any other volumetric dataset that
  // you could have loaded from any other file format.
  //
  // Software Guide : EndLatex



}

 void AnglePictureUtility::SegmentBloodVesselsFromVti(string fileName,string outputFileName,double sigma,double alpha1,double alpha2)
 {
	 vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	 reader->SetFileName(fileName.c_str());
	 reader->Update();

	 const unsigned int Dimension = 3;
	 typedef signed short                    PixelType;
     typedef itk::Image< PixelType, Dimension > ImageType;
 

	 typedef itk::VTKImageToImageFilter< ImageType > FilterType;
     FilterType::Pointer filter = FilterType::New();
	 reader->Update();
	 filter->SetInput( reader->GetOutput() );
	 filter->Update();
	 cout << "55555555555555555555555555555555555555" << endl;
	 ImageType* image = filter->GetOutput();

	 typedef itk::HessianRecursiveGaussianImageFilter< ImageType >
    HessianFilterType;
  HessianFilterType::Pointer hessianFilter = HessianFilterType::New();
    if( sigma )
    {
    //hessianFilter->SetSigma( 1  );
    }
  hessianFilter->SetInput( image );
  hessianFilter->Update();
  cout << "1111111111111111111111111111111111111" << endl;


  typedef itk::Hessian3DToVesselnessMeasureImageFilter< PixelType >
    VesselnessMeasureFilterType;
  VesselnessMeasureFilterType::Pointer vesselnessFilter = VesselnessMeasureFilterType::New();
  vesselnessFilter->SetInput( hessianFilter->GetOutput() );
  
  if( alpha1 )
    {
    //vesselnessFilter->SetAlpha1( 1  );
    }
  if( alpha2 )
    {
    //vesselnessFilter->SetAlpha2(  1 );
    }
   
  vesselnessFilter->Update();
   cout << "222222222222222222222222222222222222222222" << endl;
   /**
    typedef itk::GDCMImageIO           ImageIOType;
  ImageIOType::Pointer gdcmImageIO = ImageIOType::New();
  typedef itk::ImageFileWriter< ImageType >  Writer1Type;
  Writer1Type::Pointer writer1 = Writer1Type::New();
  writer1->SetImageIO(gdcmImageIO);
  writer1->UseInputMetaDataDictionaryOff();
  writer1->SetInput(vesselnessFilter->GetOutput());
  writer1->SetFileName(outputFileName);
  writer1->Update();
  */
    typedef itk::ImageToVTKImageFilter< ImageType > FilterTypeToVtk;
  FilterTypeToVtk::Pointer filterTovtk = FilterTypeToVtk::New();
  filterTovtk->SetInput( vesselnessFilter->GetOutput() );
  filterTovtk->Update();

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(outputFileName.c_str());
  writer->SetInputData(filterTovtk->GetOutput());
  writer->Update();

 }

 void AnglePictureUtility::CurvesLevelSetImage(string inputFileName,int x,int y,int z)
 {

	subLog.info("AnglePictureUtility::CurvesLevelSetImage()开始执行");
	rootLog.info("AnglePictureUtility::CurvesLevelSetImage()开始执行");
	   //  Software Guide : BeginLatex
  //
  //  We now define the image type using a particular pixel type and
  //  dimension. In this case the \code{float} type is used for the pixels
  //  due to the requirements of the smoothing filter.
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet

  //第一个疑问，InternalPixelType，它代表哪个方向来的数据
  typedef   float           InternalPixelType;
  const     unsigned int    Dimension = 3;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
  // Software Guide : EndCodeSnippet
  //  The following lines instantiate the thresholding filter that will
  //  process the final level set at the output of the
  //  CurvesLevelSetImageFilter.
  //
  typedef unsigned char                            OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;


  //是给谁用的
  typedef itk::BinaryThresholdImageFilter<
                        InternalImageType,
                        OutputImageType    >       ThresholdingFilterType;


  //thresholder是在最后的时候用，关于它的参数设置，要怎么弄，水平集跑完之后弄的，等于对应的是最后出来的结果
  ThresholdingFilterType::Pointer thresholder = ThresholdingFilterType::New();
  thresholder->SetLowerThreshold( -1000.0 );
  thresholder->SetUpperThreshold(     0.0 );
  thresholder->SetOutsideValue(  0  );
  thresholder->SetInsideValue(  255 );
  // We instantiate reader and writer types in the following lines.
  //

  /**
  typedef  itk::ImageFileReader< InternalImageType > ReaderType;
  typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  */



  //  The RescaleIntensityImageFilter type is declared below. This filter will
  //  renormalize image before sending them to writers.
  //
  typedef itk::RescaleIntensityImageFilter<
                               InternalImageType,
                               OutputImageType >   CastFilterType;
  //  The \doxygen{CurvatureAnisotropicDiffusionImageFilter} type is
  //  instantiated using the internal image type.
  //
  typedef   itk::CurvatureAnisotropicDiffusionImageFilter<
                               InternalImageType,
                               InternalImageType >  SmoothingFilterType;
  SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();
  //  The types of the
  //  GradientMagnitudeRecursiveGaussianImageFilter and
  //  SigmoidImageFilter are instantiated using the internal image 
  //  type.
  //
  typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                               InternalImageType,
                               InternalImageType >  GradientFilterType;
  typedef   itk::SigmoidImageFilter<
                               InternalImageType,
                               InternalImageType >  SigmoidFilterType;
  GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
  SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
  //  The minimum and maximum values of the SigmoidImageFilter output
  //  are defined with the methods \code{SetOutputMinimum()} and
  //  \code{SetOutputMaximum()}. In our case, we want these two values to be
  //  $0.0$ and $1.0$ respectively in order to get a nice speed image to feed
  //  the \code{FastMarchingImageFilter}. Additional details on the user of the
  //  \doxygen{SigmoidImageFilter} are presented in
  //  section~\ref{sec:IntensityNonLinearMapping}.
  sigmoid->SetOutputMinimum(  0.0  );
  sigmoid->SetOutputMaximum(  1.0  );
  //  We declare now the type of the FastMarchingImageFilter that
  //  will be used to generate the initial level set in the form of a distance
  //  map.
  //
  typedef  itk::FastMarchingImageFilter<
                              InternalImageType,
                              InternalImageType >    FastMarchingFilterType;
  //  Next we construct one filter of this class using the \code{New()}
  //  method.
  //
  FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();
  //  Software Guide : BeginLatex
  //
  //  In the following lines we instantiate the type of the
  //  CurvesLevelSetImageFilter and create an object of this
  //  type using the \code{New()} method.
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  typedef  itk::CurvesLevelSetImageFilter< InternalImageType,
                InternalImageType >    CurvesFilterType;
  CurvesFilterType::Pointer geodesicActiveContour =
                                     CurvesFilterType::New();
  // Software Guide : EndCodeSnippet
  //  Software Guide : BeginLatex
  //
  //  For the CurvesLevelSetImageFilter, scaling parameters
  //  are used to trade off between the propagation (inflation), the
  //  curvature (smoothing) and the advection terms. These parameters are set
  //  using methods \code{SetPropagationScaling()},
  //  \code{SetCurvatureScaling()} and \code{SetAdvectionScaling()}. In this
  //  example, we will set the curvature and advection scales to one and let
  //  the propagation scale be a command-line argument.
  //
  //  \index{itk::Geodesic\-Active\-Contour\-LevelSet\-Image\-Filter!SetPropagationScaling()}
  //  \index{itk::Segmentation\-Level\-Set\-Image\-Filter!SetPropagationScaling()}
  //  \index{itk::Geodesic\-Active\-Contour\-LevelSet\-Image\-Filter!SetCurvatureScaling()}
  //  \index{itk::Segmentation\-Level\-Set\-Image\-Filter!SetCurvatureScaling()}
  //  \index{itk::Geodesic\-Active\-Contour\-LevelSet\-Image\-Filter!SetAdvectionScaling()}
  //  \index{itk::Segmentation\-Level\-Set\-Image\-Filter!SetAdvectionScaling()}
  //
  //  Software Guide : EndLatex
  //原来要设置的参数
  const double propagationScaling = 1.0;//之前没有设置，我到现在都看不懂
  //  Software Guide : BeginCodeSnippet
  geodesicActiveContour->SetPropagationScaling( propagationScaling );
  geodesicActiveContour->SetCurvatureScaling( 1.0 );
  geodesicActiveContour->SetAdvectionScaling( 1.0 );
  //  Software Guide : EndCodeSnippet
  //  Once activiated the level set evolution will stop if the convergence
  //  criteria or if the maximum number of iterations is reached.  The
  //  convergence criteria is defined in terms of the root mean squared (RMS)
  //  change in the level set function. The evolution is said to have
  //  converged if the RMS change is below a user specified threshold.  In a
  //  real application is desirable to couple the evolution of the zero set
  //  to a visualization module allowing the user to follow the evolution of
  //  the zero set. With this feedback, the user may decide when to stop the
  //  algorithm before the zero set leaks through the regions of low gradient
  //  in the contour of the anatomical structure to be segmented.
  geodesicActiveContour->SetMaximumRMSError( 0.02 );
  geodesicActiveContour->SetNumberOfIterations( 800 );
  //  Software Guide : BeginLatex
  //
  //  The filters are now connected in a pipeline indicated in
  //  Figure~\ref{fig:CurvessCollaborationDiagram} using the
  //  following lines:
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet



  //这里就是关键性的代码
  
 // itk::SmartPointer<InternalImageType> lineimage;
  
 // {
  vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  reader->SetFileName(inputFileName.c_str());
  reader->Update();
  	subLog.info("AnglePictureUtility::CurvesLevelSetImage()从vti中读取图像完成");
	rootLog.info("AnglePictureUtility::CurvesLevelSetImage()从vti中读取图像完成");
	

	//这里我大概还记的是要做数据类型的转换
	vtkImageData* original = reader->GetOutput();
	vtkSmartPointer<vtkImageCast> castToFloat = vtkSmartPointer<vtkImageCast>::New();
	castToFloat->SetInputData(original);
	castToFloat->SetOutputScalarTypeToFloat();
	castToFloat->Update();
	vtkImageData* imageData = castToFloat->GetOutput(); 
	  	subLog.info("AnglePictureUtility::CurvesLevelSetImage()把图像数据类型从signed short转化到float");
	rootLog.info("AnglePictureUtility::CurvesLevelSetImage()把图像数据类型从signed short转化到float");

	 typedef itk::VTKImageToImageFilter< InternalImageType > FilterType;
     FilterType::Pointer filter = FilterType::New();
	 filter->SetInput( imageData );
	 filter->Update();
	// InternalImageType* image = filter->GetOutput();
//	 cout << filter->GetOutput()->GetReferenceCount() << endl;
	 InternalImageType*          lineimage = filter->GetOutput();
//	 cout << lineimage->GetReferenceCount() << endl;
	subLog.info("AnglePictureUtility::CurvesLevelSetImage()转换成itk图像");
	rootLog.info("AnglePictureUtility::CurvesLevelSetImage()转换成itk图像");
 // }



  //writer->SetInput( thresholder->GetOutput() );








  // Software Guide : EndCodeSnippet
  //  The CurvatureAnisotropicDiffusionImageFilter requires a couple of
  //  parameter to be defined. The following are typical values for $2D$
  //  images. However they may have to be adjusted depending on the amount of
  //  noise present in the input image. This filter has been discussed in
  //  section~\ref{sec:GradientAnisotropicDiffusionImageFilter}.

  //每一个参数给出的意义
 /*Conductance Parameter
 * The conductance parameter controls the sensitivity of the conductance term
 * in the basic anisotropic diffusion equation.  It affects the conductance term
 * in different ways depending on the particular variation on the basic
 * equation. As a general rule, the lower the value, the more strongly the
 * diffusion equation preserves image features (such as high gradients or
 * curvature).  A high value for conductance will cause the filter to diffuse
 * image features more readily.  Typical values range from 0.5 to 2.0 for data
 * like the Visible Human color data, but the correct value for your
 * application is wholly dependent on the results you want from a specific data
 * set and the number or iterations you perform.
 */
  smoothing->SetTimeStep( 0.02 );
  smoothing->SetNumberOfIterations(  5 );
  smoothing->SetConductanceParameter( 2.0 );
  //  The GradientMagnitudeRecursiveGaussianImageFilter performs the
  //  equivalent of a convolution with a Gaussian kernel, followed by a
  //  derivative operator. The sigma of this Gaussian can be used to control
  //  the range of influence of the image edges. This filter has been discussed
  //  in Section~\ref{sec:GradientMagnitudeRecursiveGaussianImageFilter}.


  //上面的注释里有说明
  const double sigma = 0.4; //原来是参数
  gradientMagnitude->SetSigma(  sigma  );
  //  The SigmoidImageFilter requires two parameters that define the linear
  //  transformation to be applied to the sigmoid argument. This parameters
  //  have been discussed in Sections~\ref{sec:IntensityNonLinearMapping} and
  //  \ref{sec:FastMarchingImageFilter}.

  //alpha,beta的参数选择怎么做
  const double alpha =  5;
  const double beta  =  10;
  sigmoid->SetAlpha( alpha );
  sigmoid->SetBeta(  beta  );
  //  The FastMarchingImageFilter requires the user to provide a seed
  //  point from which the level set will be generated. The user can actually
  //  pass not only one seed point but a set of them. Note the the
  //  FastMarchingImageFilter is used here only as a helper in the
  //  determination of an initial level set. We could have used the
  //  \doxygen{DanielssonDistanceMapImageFilter} in the same way.
  //
  //  The seeds are passed stored in a container. The type of this
  //  container is defined as \code{NodeContainer} among the
  //  FastMarchingImageFilter traits.
  //
  typedef FastMarchingFilterType::NodeContainer  NodeContainer;
  typedef FastMarchingFilterType::NodeType       NodeType;
  NodeContainer::Pointer seeds = NodeContainer::New();
  InternalImageType::IndexType  seedPosition;
  seedPosition[0] = x;
  seedPosition[1] = y;
  seedPosition[2] = z;
  //  Nodes are created as stack variables and initialized with a value and an
  //  \doxygen{Index} position. Note that here we assign the value of minus the
  //  user-provided distance to the unique node of the seeds passed to the
  //  FastMarchingImageFilter. In this way, the value will increment
  //  as the front is propagated, until it reaches the zero value corresponding
  //  to the contour. After this, the front will continue propagating until it
  //  fills up the entire image. The initial distance is taken here from the
  //  command line arguments. The rule of thumb for the user is to select this
  //  value as the distance from the seed points at which she want the initial
  //  contour to be.

  //这个值如何选择,spacing,是多少，我可以参照不同数量级，
  const double initialDistance = 0.05;
  NodeType node;
  const double seedValue = - initialDistance;
  node.SetValue( seedValue );
  node.SetIndex( seedPosition );
  //  The list of nodes is initialized and then every node is inserted using
  //  the \code{InsertElement()}.
  seeds->Initialize();
  seeds->InsertElement( 0, node );
  //  The set of seed nodes is passed now to the
  //  FastMarchingImageFilter with the method
  //  \code{SetTrialPoints()}.
  //
  fastMarching->SetTrialPoints(  seeds  );
  //  Since the FastMarchingImageFilter is used here just as a
  //  Distance Map generator. It does not require a speed image as input.
  //  Instead the constant value $1.0$ is passed using the
  //  \code{SetSpeedConstant()} method.
  //
  fastMarching->SetSpeedConstant( 1.0 );
  //  Here we configure all the writers required to see the intermediate
  //  outputs of the pipeline. This is added here only for
  //  pedagogical/debugging purposes. These intermediate output are normaly not
  //  required. Only the output of the final thresholding filter should be
  //  relevant.  Observing intermediate output is helpful in the process of
  //  fine tuning the parameters of filters in the pipeline.
  //

  smoothing->SetInput( lineimage );
  smoothing->Update();
  subLog.info("AnglePictureUtility::CurvesLevelSetImage()完成图像平滑的操作");
  rootLog.info("AnglePictureUtility::CurvesLevelSetImage()完成图像平滑的操作");


  gradientMagnitude->SetInput( smoothing->GetOutput() );
  gradientMagnitude->Update();
  subLog.info("AnglePictureUtility::CurvesLevelSetImage()完成图像梯度的操作");
  rootLog.info("AnglePictureUtility::CurvesLevelSetImage()完成图像梯度的操作");


  sigmoid->SetInput( gradientMagnitude->GetOutput() );
  sigmoid->Update();
  subLog.info("AnglePictureUtility::CurvesLevelSetImage()完成图像sigmoid的操作");
  rootLog.info("AnglePictureUtility::CurvesLevelSetImage()完成图像sigmoid的操作");

  fastMarching->Update();
  subLog.info("AnglePictureUtility::CurvesLevelSetImage()  完成图像fastMarching操作");
  rootLog.info("AnglePictureUtility::CurvesLevelSetImage() 完成图像fastMarching的操作");


   sigmoid->GetOutput()->SetSpacing(fastMarching->GetOutput()->GetSpacing());
   //sigmoid->GetOutput()->SetRequestedRegion(sigmoid->GetOutput()->GetLargestPossibleRegion());
   //fastMarching->GetOutput()->SetRequestedRegion(fastMarching->GetOutput()->GetLargestPossibleRegion());
   //sigmoid->GetOutput()->SetRequestedRegion(fastMarching->GetOutput()->GetRequestedRegion());
  //fastMarching->GetOutput()->SetSpacing(sigmoid->GetOutput()->GetSpacing());
  //fastMarching->GetOutput()->SetRequestedRegion(sigmoid->GetOutput()->GetRequestedRegion());
  // fastMarching->UpdateLargestPossibleRegion();
 //  sigmoid->UpdateLargestPossibleRegion();

  geodesicActiveContour->SetInput(  fastMarching->GetOutput() );
  geodesicActiveContour->SetFeatureImage( sigmoid->GetOutput() );
  
 
  try
  {
	geodesicActiveContour->UpdateLargestPossibleRegion();
    geodesicActiveContour->Update();
  }catch(itk::ExceptionObject&exp)
  {
	  string based = "AnglePictureUtility::CurvesLevelSetImage()";
	  string res = based +exp.GetDescription();
	  subLog.info(res);
	  rootLog.info(res);
	  return;
  }
  //geodesicActiveContour->GetOutput();curveLevelSet,有时候成功，有时候失败，为什么失败，不知道
  subLog.info("AnglePictureUtility::CurvesLevelSetImage()完成levelset的操作");
  rootLog.info("AnglePictureUtility::CurvesLevelSetImage()完成levelset的操作");
  thresholder->SetInput( geodesicActiveContour->GetOutput() );  
  thresholder->Update();
    subLog.info("AnglePictureUtility::CurvesLevelSetImage()完成thresholder的操作");
  rootLog.info("AnglePictureUtility::CurvesLevelSetImage()完成thresholder的操作");

  /**
  CastFilterType::Pointer caster1 = CastFilterType::New();
  CastFilterType::Pointer caster2 = CastFilterType::New();
  CastFilterType::Pointer caster3 = CastFilterType::New();
  CastFilterType::Pointer caster4 = CastFilterType::New();
  WriterType::Pointer writer1 = WriterType::New();
  WriterType::Pointer writer2 = WriterType::New();
  WriterType::Pointer writer3 = WriterType::New();
  WriterType::Pointer writer4 = WriterType::New();
  caster1->SetInput( smoothing->GetOutput() );
  writer1->SetInput( caster1->GetOutput() );
  writer1->SetFileName("CurvesImageFilterOutput1.png");
  caster1->SetOutputMinimum(   0 );
  caster1->SetOutputMaximum( 255 );
  writer1->Update();
  caster2->SetInput( gradientMagnitude->GetOutput() );
  writer2->SetInput( caster2->GetOutput() );
  writer2->SetFileName("CurvesImageFilterOutput2.png");
  caster2->SetOutputMinimum(   0 );
  caster2->SetOutputMaximum( 255 );
  writer2->Update();
  caster3->SetInput( sigmoid->GetOutput() );
  writer3->SetInput( caster3->GetOutput() );
  writer3->SetFileName("CurvesImageFilterOutput3.png");
  caster3->SetOutputMinimum(   0 );
  caster3->SetOutputMaximum( 255 );
  writer3->Update();
  caster4->SetInput( fastMarching->GetOutput() );
  writer4->SetInput( caster4->GetOutput() );
  writer4->SetFileName("CurvesImageFilterOutput4.png");
  caster4->SetOutputMinimum(   0 );
  caster4->SetOutputMaximum( 255 );
  //  The FastMarchingImageFilter requires the user to specify the
  //  size of the image to be produced as output. This is done using the
  //  \code{SetOutputSize()}. Note that the size is obtained here from the
  //  output image of the smoothing filter. The size of this image is valid
  //  only after the \code{Update()} methods of this filter has been called
  //  directly or indirectly.
  //
  fastMarching->SetOutputSize(
           reader->GetOutput()->GetBufferedRegion().GetSize() );
  //  Software Guide : BeginLatex
  //
  //  The invocation of the \code{Update()} method on the writer triggers the
  //  execution of the pipeline.  As usual, the call is placed in a
  //  \code{try/catch} block should any errors occur or exceptions be thrown.
  //
  //  Software Guide : EndLatex
  // Software Guide : BeginCodeSnippet
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  // Software Guide : EndCodeSnippet
  // Print out some useful information
  std::cout << std::endl;
  std::cout << "Max. no. iterations: " << geodesicActiveContour->GetNumberOfIterations() << std::endl;
  std::cout << "Max. RMS error: " << geodesicActiveContour->GetMaximumRMSError() << std::endl;
  std::cout << std::endl;
  std::cout << "No. elpased iterations: " << geodesicActiveContour->GetElapsedIterations() << std::endl;
  std::cout << "RMS change: " << geodesicActiveContour->GetRMSChange() << std::endl;
  writer4->Update();
  // The following writer type is used to save the output of the time-crossing
  // map in a file with apropiate pixel representation. The advantage of saving
  // this image in native format is that it can be used with a viewer to help
  // determine an appropriate threshold to be used on the output of the
  // fastmarching filter.
  //
  typedef itk::ImageFileWriter< InternalImageType > InternalWriterType;
  InternalWriterType::Pointer mapWriter = InternalWriterType::New();
  mapWriter->SetInput( fastMarching->GetOutput() );
  mapWriter->SetFileName("CurvesImageFilterOutput4.mha");
  mapWriter->Update();
  InternalWriterType::Pointer speedWriter = InternalWriterType::New();
  speedWriter->SetInput( sigmoid->GetOutput() );
  speedWriter->SetFileName("CurvesImageFilterOutput3.mha");
  speedWriter->Update();
  InternalWriterType::Pointer gradientWriter = InternalWriterType::New();
  gradientWriter->SetInput( gradientMagnitude->GetOutput() );
  gradientWriter->SetFileName("CurvesImageFilterOutput2.mha");
  gradientWriter->Update();
  //  Software Guide : BeginLatex
  //
  //  Let's now run this example using as input the image
  //  \code{BrainProtonDensitySlice.png} provided in the directory
  //  \code{Examples/Data}. We can easily segment the major anatomical
  //  structures by providing seeds in the appropriate locations.
  //  Table~\ref{tab:CurvesImageFilterOutput2} presents the
  //  parameters used for some structures.
  //
  //  \begin{table}
  //  \begin{center}
  //  \begin{tabular}{|l|c|c|c|c|c|c|c|c|}
  //  \hline
  //  Structure    & Seed Index &  Distance   &   $\sigma$  &
  //  $\alpha$     &  $\beta$   & Propag. & Output Image \\  \hline
  //  Left Ventricle  & $(81,114)$ & 5.0 & 1.0 & -0.5 & 3.0  &  2.0 & First   \\  \hline
  //  Right Ventricle & $(99,114)$ & 5.0 & 1.0 & -0.5 & 3.0  &  2.0 & Second  \\  \hline
  //  White matter    & $(56, 92)$ & 5.0 & 1.0 & -0.3 & 2.0  & 10.0 & Third   \\  \hline
  //  Gray matter     & $(40, 90)$ & 5.0 & 0.5 & -0.3 & 2.0  & 10.0 & Fourth  \\  \hline
  //  \end{tabular}
  //  \end{center}
  //  \itkcaption[Curves segmentation example parameters]{Parameters used
  //  for segmenting some brain structures shown in
  //  Figure~\ref{fig:CurvesImageFilterOutput2} using the filter
  //  CurvesLevelSetImageFilter.
  //  \label{tab:CurvesImageFilterOutput2}}
  //  \end{table}
  //
  //  Figure~\ref{fig:CurvesImageFilterOutput} presents the
  //  intermediate outputs of the pipeline illustrated in
  //  Figure~\ref{fig:CurvessCollaborationDiagram}. They are
  //  from left to right: the output of the anisotropic diffusion filter, the
  //  gradient magnitude of the smoothed image and the sigmoid of the gradient
  //  magnitude which is finally used as the edge potential for the
  //  CurvesLevelSetImageFilter.
  //
  // \begin{figure} \center
  // \includegraphics[height=0.40\textheight]{BrainProtonDensitySlice}
  // \includegraphics[height=0.40\textheight]{CurvesImageFilterOutput1}
  // \includegraphics[height=0.40\textheight]{CurvesImageFilterOutput2}
  // \includegraphics[height=0.40\textheight]{CurvesImageFilterOutput3}
  // \itkcaption[CurvesLevelSetImageFilter intermediate
  // output]{Images generated by the segmentation process based on the
  // CurvesLevelSetImageFilter. From left to right and top to
  // bottom: input image to be segmented, image smoothed with an
  // edge-preserving smoothing filter, gradient magnitude of the smoothed
  // image, sigmoid of the gradient magnitude. This last image, the sigmoid, is
  // used to compute the speed term for the front propagation.}
  // \label{fig:CurvesImageFilterOutput} \end{figure}
  //
  //  Segmentations of the main brain structures are presented in
  //  Figure~\ref{fig:CurvesImageFilterOutput2}. The results
  //  are quite similar to those obtained with the
  //  ShapeDetectionLevelSetImageFilter in
  //  Section~\ref{sec:ShapeDetectionLevelSetFilter}.
  //
  //  Note that a relatively larger propagation scaling value was required to
  //  segment the white matter. This is due to two factors: the lower
  //  contrast at the border of the white matter and the complex shape of the
  //  structure. Unfortunately the optimal value of these scaling parameters
  //  can only be determined by experimentation. In a real application we
  //  could imagine an interactive mechanism by which a user supervises the
  //  contour evolution and adjusts these parameters accordingly.
  //
  // \begin{figure} \center
  // \includegraphics[width=0.24\textwidth]{CurvesImageFilterOutput5}
  // \includegraphics[width=0.24\textwidth]{CurvesImageFilterOutput6}
  // \includegraphics[width=0.24\textwidth]{CurvesImageFilterOutput7}
  // \includegraphics[width=0.24\textwidth]{CurvesImageFilterOutput8}
  // \itkcaption[CurvesImageFilter segmentations]{Images generated by the
  // segmentation process based on the CurvesImageFilter. From left to
  // right: segmentation of the left ventricle, segmentation of the right
  // ventricle, segmentation of the white matter, attempt of segmentation of
  // the gray matter.}
  // \label{fig:CurvesImageFilterOutput2}
  // \end{figure}
  //
  //  Software Guide : EndLatex
  */

 }
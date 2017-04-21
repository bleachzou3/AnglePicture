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

bool AnglePictureUtility::segment(string directoryName,string outputFileName,int x,int y,int z,float lowerThreshold,float uppperThreshold)
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
	/**
	for(int z = 0; z < size[2] ; z++)
	{
		for(int x = 0; x < size[0]; x++)
		{
			for(int y = 0; y < size[1]; y++)
			{
				pixelIndex[0] = x;
				pixelIndex[1] = y;
				pixelIndex[2] = z;
				high = std::max(high,image->GetPixel(pixelIndex));
				low = std::min(low,image->GetPixel(pixelIndex));
			}
		}
	}
	cout << "hight" << high << "   low:" << low << endl;
	*/

	
	pixelIndex[0] = x;
	pixelIndex[1] = y;
	pixelIndex[2] = z;
	
	//今天下午
	typedef itk::ConnectedThresholdImageFilter< ImageType,
		ImageType > ConnectedFilterType;

	 ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
	 connectedThreshold->SetSeed(pixelIndex);
	 connectedThreshold->SetLower(lowerThreshold);
	 connectedThreshold->SetUpper(uppperThreshold);























// Software Guide : EndCodeSnippet
// Software Guide : BeginLatex
//
// At this point, we have a volumetric image in memory that we can access by
// invoking the \code{GetOutput()} method of the reader.
//
// Software Guide : EndLatex
// Software Guide : BeginLatex
//
// We proceed now to save the volumetric image in another file, as specified by
// the user in the command line arguments of this program. Thanks to the
// ImageIO factory mechanism, only the filename extension is needed to identify
// the file format in this case.
//
// Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
	writer->SetFileName( outputFileName );
	writer->SetInput( connectedThreshold->GetOutput() );
// Software Guide : EndCodeSnippet
    std::cout  << "Writing the image as " << std::endl << std::endl;
	std::cout  << outputFileName << std::endl << std::endl;
// Software Guide : BeginLatex
//
// The process of writing the image is initiated by invoking the
// \code{Update()} method of the writer.
//
// Software Guide : EndLatex
    try
      {
// Software Guide : BeginCodeSnippet
     // writer->Update();
// Software Guide : EndCodeSnippet
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
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
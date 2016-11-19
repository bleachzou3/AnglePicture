#include <vtkSmartPointer.h>
#include <vtkVersion.h>

#include <vtkParametricFunctionSource.h>
#include <vtkTupleInterpolator.h>
#include <vtkTubeFilter.h>
#include <vtkParametricSpline.h>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <map>
#include <unordered_map>
#include <vtkLine.h>
using namespace std;
void computeNewRadius(vtkPolyData * data)
{
	vtkFloatArray* radius =static_cast<vtkFloatArray*>( data->GetPointData()->GetArray("MaximumInscribedSphereRadius"));


	int N = data->GetNumberOfPoints();

	vtkSmartPointer<vtkFloatArray> newRadius = vtkSmartPointer<vtkFloatArray>::New();
	newRadius->SetName("NewMaximumInscribedSphereRadius");

	newRadius->SetNumberOfComponents(1);
	newRadius->SetNumberOfTuples(N);


	for(int i = 0; i < N; i++)
	{
		newRadius->SetTuple1(1,radius->GetTuple1(i)*10);
	}



	data->GetPointData()->AddArray(newRadius);

}
void computeSegment(vtkPolyData*centerLine,unordered_map<int,vtkPolyData*>& centerLineCell);
int main(int, char *[])
{
   

 /**
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->InsertPoint(0,1,0,0);
  
  points->InsertPoint(1,2,0,0);
  points->InsertPoint(2,3,1,0);
  points->InsertPoint(3,4,1,0);
  points->InsertPoint(4,5,0,0);
  points->InsertPoint(5,6,0,0);

  // Fit a spline to the points
  vtkSmartPointer<vtkParametricSpline> spline =
    vtkSmartPointer<vtkParametricSpline>::New();
  spline->SetPoints(points);
  vtkSmartPointer<vtkParametricFunctionSource> functionSource =
    vtkSmartPointer<vtkParametricFunctionSource>::New();
  functionSource->SetParametricFunction(spline);
  functionSource->SetUResolution(10 * points->GetNumberOfPoints());
  functionSource->Update();

  // Interpolate the scalars
  double rad;
  vtkSmartPointer<vtkTupleInterpolator> interpolatedRadius =
    vtkSmartPointer<vtkTupleInterpolator> ::New();
  interpolatedRadius->SetInterpolationTypeToLinear();
  interpolatedRadius->SetNumberOfComponents(1);
  rad = .2; interpolatedRadius->AddTuple(0,&rad);
  rad = .2; interpolatedRadius->AddTuple(1,&rad);
  rad = .2; interpolatedRadius->AddTuple(2,&rad);
  rad = .1; interpolatedRadius->AddTuple(3,&rad);
  rad = .1; interpolatedRadius->AddTuple(4,&rad);
  rad = .1; interpolatedRadius->AddTuple(5,&rad);

  // Generate the radius scalars
  vtkSmartPointer<vtkDoubleArray> tubeRadius =
    vtkSmartPointer<vtkDoubleArray>::New();
  unsigned int n = functionSource->GetOutput()->GetNumberOfPoints();
  tubeRadius->SetNumberOfTuples(n);
  tubeRadius->SetName("TubeRadius");
  double tMin = interpolatedRadius->GetMinimumT();
  double tMax = interpolatedRadius->GetMaximumT();
  double r;
  for (unsigned int i = 0; i < n; ++i)
    {
    double t = (tMax - tMin) / (n - 1) * i + tMin;
    interpolatedRadius->InterpolateTuple(t, &r);
    tubeRadius->SetTuple1(i, r);
    }
	

  // Add the scalars to the polydata

  vtkSmartPointer<vtkPolyData> tubePolyData =
    vtkSmartPointer<vtkPolyData>::New();
  tubePolyData = functionSource->GetOutput();
  tubePolyData->GetPointData()->AddArray(tubeRadius);
  tubePolyData->GetPointData()->SetActiveScalars("TubeRadius");
  */



    vtkSmartPointer<vtkXMLPolyDataReader> centerlineModelReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	//centerlineModelReader->SetFileName("E:\\model0927centerline.vtp");
	centerlineModelReader->SetFileName("E:\\VMTKCenterlinesOut.vtp");
  centerlineModelReader->Update();
  vtkSmartPointer<vtkPolyData>centerline = centerlineModelReader->GetOutput();
  //computeNewRadius(centerline);
  centerline->GetPointData()->SetActiveScalars("MaximumInscribedSphereRadius");

  /*
  // Create the tubes
  vtkSmartPointer<vtkTubeFilter> tuber =
    vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  tuber->SetInput(centerline);
#else
  tuber->SetInputData(centerline);
#endif
  //tuber->SetRadius(20);
  //tuber->SetRadiusFactor(20);
  tuber->SetNumberOfSides(100);
  tuber->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
  tuber->Update();
  //--------------
  // Setup actors and mappers
  vtkSmartPointer<vtkPolyDataMapper> lineMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  lineMapper->SetInput(centerline);
#else
  lineMapper->SetInputData(centerline);
#endif
  lineMapper->SetScalarRange(centerline->GetScalarRange());

  vtkSmartPointer<vtkPolyDataMapper> tubeMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  vtkPolyData*tubeData = tuber->GetOutput();
  cout << "tubeData" << tubeData->GetNumberOfCells()<<endl;
  tubeMapper->SetInputConnection(tuber->GetOutputPort());
  double range[2];
  centerline->GetScalarRange(range);
  cout << "range............" << range[0] <<"  "<< range[1] << endl;
  tubeMapper->SetScalarRange(centerline->GetScalarRange());

  vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
  lineActor->SetMapper(lineMapper);
  vtkSmartPointer<vtkActor> tubeActor = vtkSmartPointer<vtkActor>::New();
  tubeActor->SetMapper(tubeMapper);
  */

    vtkSmartPointer<vtkActor> lineActor = vtkSmartPointer<vtkActor>::New();
	vtkSmartPointer<vtkPolyDataMapper> lineMapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
	lineMapper->SetInputData(centerline);
    lineActor->SetMapper(lineMapper);
	lineMapper->SetScalarRange(centerline->GetScalarRange());



  // Setup render window, renderer, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  unordered_map<int,vtkPolyData*> segments;
	computeSegment(centerline,segments);
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
  	for(int i = 0; i < segments.size(); i++)
	{
		  vtkSmartPointer<vtkTubeFilter> tuber = vtkSmartPointer<vtkTubeFilter>::New();
		  tuber->SetNumberOfSides(20);
          tuber->SetVaryRadiusToVaryRadiusByAbsoluteScalar();
		  segments[i]->GetPointData()->SetActiveScalars("MaximumInscribedSphereRadius");
		  tuber->SetInputData(segments[i]);
		  tuber->Update();
		  cout << "i" << i << "numberOfCells" << segments[i]->GetNumberOfCells()<<endl;
		  tuber->Update();
		  vtkSmartPointer<vtkPolyDataMapper> tubeMapper =
          vtkSmartPointer<vtkPolyDataMapper>::New();
		  tubeMapper->SetInputData(tuber->GetOutput());
		  tubeMapper->SetScalarRange(segments[i]->GetScalarRange());
		  vtkSmartPointer<vtkActor> tubeActor = vtkSmartPointer<vtkActor>::New();
          tubeActor->SetMapper(tubeMapper);
		  renderer->AddActor(tubeActor);

	}
  renderer->AddActor(lineActor);
  //renderer->AddActor(tubeActor);
  renderer->SetBackground(.4, .5, .6);
  renderWindow->Render();
  renderWindowInteractor->Start();
  for(int i = 0; i < segments.size(); i++)
  {
	  if(segments[i]->GetReferenceCount() == 1)
	  {
		  segments[i]->Delete();
	  }
    
  }
  return EXIT_SUCCESS;
}
/**
*给出血管中心线
*/
void computeSegment(vtkPolyData*centerLine,unordered_map<int,vtkPolyData*>& centerLineCell)
{
 
  centerLine->GetLines()->InitTraversal();
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  int numCell = 0;
  while(centerLine->GetLines()->GetNextCell(idList))
    {
         std::cout << "Line has " << idList->GetNumberOfIds() << " points." << std::endl;
	     vtkSmartPointer<vtkPoints> curPoints = vtkSmartPointer<vtkPoints>::New();
	     vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
		 for(vtkIdType pointId = 0; pointId < idList->GetNumberOfIds(); pointId++)
		 {
		   std::cout << idList->GetId(pointId) << " ";
		   double singlePoint[3];
		   centerLine->GetPoint(idList->GetId(pointId),singlePoint);
		  // curPoints->SetPoint(idList->GetId(pointId),singlePoint);
		   curPoints->InsertPoint(idList->GetId(pointId),singlePoint);
		   if(pointId < idList->GetNumberOfIds() - 1)
		   {
				vtkSmartPointer<vtkLine> line =vtkSmartPointer<vtkLine>::New();
				line->GetPointIds()->SetId(0,idList->GetId(pointId));
				line->GetPointIds()->SetId(1,idList->GetId(pointId + 1));
				lines->InsertNextCell(line);
		   }
		 }
		 vtkPolyData* curSegment = vtkPolyData::New();
		 curSegment->SetPoints(curPoints);
		 curSegment->SetLines(lines);
		 curSegment->GetPointData()->AddArray(centerLine->GetPointData()->GetArray("MaximumInscribedSphereRadius"));
		 centerLineCell[numCell++] = curSegment;

    }
  
}
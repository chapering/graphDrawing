
#include "vtkDataSetAttributes.h"
#include "vtkGraphLayoutView.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkRandomGraphSource.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStringArray.h"
#include "vtkTree.h"
#include "vtkViewTheme.h"
#include "vtkHierarchicalGraphView.h"
#include "vtkSmartPointer.h"
#include "vtkCosmicTreeLayoutStrategy.h"
#include "vtkSplineGraphEdges.h"
#include "vtkRenderedHierarchyRepresentation.h"

#include <boost/lexical_cast.hpp>
#include "tree.h"
#include "newick_file.h"
#include <iostream>
#include <string>


#include "vtkActor.h"
#include "vtkActor2D.h"
#include "vtkDynamic2DLabelMapper.h"
#include "vtkGlyph3D.h"
#include "vtkGlyphSource2D.h"
#include "vtkGraphLayout.h"
#include "vtkGraphToPolyData.h"
#include "vtkInteractorStyleImage.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStringToNumeric.h"
#include "vtkTextProperty.h"
#include "vtkTreeLayoutStrategy.h"
#include "vtkXMLTreeReader.h"

#include "vtkIdTypeArray.h"
#include "vtkStringToNumeric.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkLookupTable.h"


#include "vtkCamera.h"

using vtkstd::string;

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()



#define CON

using namespace BiRC::treelib;
using namespace std;

int main(int argc, char* argv[])
{
  if ( argc < 2 ) return -1;

  std::auto_ptr<Tree> my_tree;
  my_tree = parse_newick_file(argv[1]);

  int nodeNum = my_tree.get()->size();

  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  vtkSmartPointer<vtkPoints> clrvalues = vtkSmartPointer<vtkPoints>::New();
  //tree->print(std::cout);
  //vtkIntArray* colors = NULL;
  vtkUnsignedCharArray* colors = NULL;

  if ( argc >= 3 ) {
	  //colors = vtkIntArray::New();
	  colors = vtkUnsignedCharArray::New();
	  colors->SetName("NodeColor");

	  ifstream ifs (argv[2]);
	  if (!ifs.is_open()) {
		  cerr << "can not load color values from " << argv[2] << ", ignored.\n";
		  return -1;
	  }

	  double vals[3];
	  int ilabTotal = 0;
	  while (ifs) {
		  ifs >> vals[0] >> vals[1] >> vals[2];
		  /*
		  cout << vals[0] << "," << vals[1] << "," << vals[2] << "\n";
		  colors->InsertNextTuple3( int(vals[0]*255), int(vals[1]*255), int(vals[2]*255) );
		  */
		  clrvalues->InsertNextPoint( vals );
		  ilabTotal ++;
	  }
	  ifs.close();

	  cerr << ilabTotal << " colors loaded from file: " << argv[2] << "\n";

	  vtkIdType total = clrvalues->GetNumberOfPoints();
	  lut->SetNumberOfTableValues( total  );

	  for (vtkIdType idx = 0; idx < total; idx++ ) {
		  clrvalues->GetPoint( idx, vals );
		  lut->SetTableValue(idx, vals);
	  }
  }

  cerr << "tree # of nodes: " << nodeNum << endl;

  VTK_CREATE(vtkIdTypeArray, dist);
  dist->SetName("distance");

  vtkMutableDirectedGraph* graph = vtkMutableDirectedGraph::New();
  int nodenum = my_tree.get()->size();
  vector< vector<int> > allEdge = my_tree.get()->get_edge();
  vtkIdType* node = new vtkIdType[nodenum];
  for ( int i = 0; i < nodenum; i ++ )
  {
    node[i] = graph->AddVertex();
  }

  for ( int i = 0; i < nodenum; i ++ )
  {
    for ( int j = 0; j < allEdge[i].size(); j ++ )
    {
        graph->AddGraphEdge( node[i], node[ allEdge[i][j] ] );
    }
  }

  if ( colors ) {
	  colors->SetNumberOfValues( nodenum );
  }

  vtkStringArray* labels = vtkStringArray::New();
  labels->SetName("Label");
  int leavecnt = 1;
  for ( int i = 0; i < nodenum; i ++ )
  {
    if ( my_tree.get()->is_leaf( i ) )
    {
      labels->InsertValue(node[i], my_tree.get()->label(i) );
	  dist->InsertNextValue(my_tree.get()->length_to_parent(i));

	  /*
	  if ( colors ) {
		  //colors->InsertNextValue( i+1 );
		  double dv[3];
		  clrvalues->GetPoint( leavecnt, dv );
		  int vals[3] = { int(dv[0]*255),  int(dv[1]*255), int(dv[2])*255 };
		  //colors->InsertValue(node[i], leavecnt );
		  //colors->InsertNextTupleValue( vals );
		  colors->InsertTupleValue(node[i], vals );
	  }

	  leavecnt ++;
	  */

    }
    else
    {
      string str = boost::lexical_cast<string>(i);
      //labels->InsertValue(node[i], str );
      labels->InsertValue(node[i], "");
	  dist->InsertNextValue(my_tree.get()->length_to_parent(i));

	  /*
	  if ( colors ) {
		  int vals[3] = {255,255,255};
		  //colors->InsertValue(node[i], 0);
		  //colors->InsertNextTupleValue( vals );
		  colors->InsertTupleValue(node[i], vals );
	  }
	  */
    }
  }

  //cout << "count of leaves: " << leavecnt << "\n";

  graph->GetVertexData()->AddArray(labels);

  graph->GetEdgeData()->AddArray(dist);

  vtkTree* tree = vtkTree::New();
  tree->CheckedShallowCopy(graph);

  cout << "Root id: " << tree->GetRoot() << "\n";
  for ( int i = 0; i < nodenum; i ++ )
  {
    if ( tree->IsLeaf(i) )
    {
	  if ( colors ) {
		  //colors->InsertNextValue( i+1 );
		  double dv[3];
		  clrvalues->GetPoint( leavecnt, dv );
		  unsigned char vals[3] = { int(dv[0]*255),  int(dv[1]*255), int(dv[2])*255 };
		  //colors->InsertValue(node[i], leavecnt );
		  //colors->InsertNextTupleValue( vals );
		  colors->InsertTupleValue(i, vals );
	  }

	  leavecnt ++;

    }
    else
    {
      string str = boost::lexical_cast<string>(i);
      //labels->InsertValue(node[i], str );
      labels->InsertValue(node[i], "");
	  dist->InsertNextValue(my_tree.get()->length_to_parent(i));

	  if ( colors ) {
		  unsigned char vals[3] = {255,255,255};
		  //colors->InsertValue(node[i], 0);
		  //colors->InsertNextTupleValue( vals );
		  colors->InsertTupleValue(i, vals );
	  }
    }
  }

  VTK_CREATE(vtkStringToNumeric, numeric);
  numeric->SetInput(tree);
  numeric->Update();

  if ( colors ) {
	  //graph->GetVertexData()->AddArray( colors );
	  tree->GetVertexData()->AddArray( colors );
	  tree->GetVertexData()->SetActiveScalars("NodeColor");
  }

//////////////////////////////
//////////////////////////////
  /*
  colors = vtkIntArray::New();
  colors->SetName("node colors");
  colors->SetNumberOfComponents(3);
  colors->SetNumberOfTuples( tree->GetNumberOfVertices() );

  for ( int i = 0; i < tree->GetNumberOfVertices(); i ++ )
  {
	  colors->SetTuple3(i, 0, 0, 0);
  }

  cout << "#tuples: " << colors->GetNumberOfTuples() << "\n";
  cout << "#vertices: " << tree->GetNumberOfVertices() << "\n";

  if ( colors ) {
	  //tree->GetVertexData()->AddArray( colors );
	  //tree->GetVertexData()->SetActiveScalars("node colors");
	  tree->GetVertexData()->SetScalars( colors );
  }
  */

  // Layout the tree using vtkGraphLayout.
  vtkGraphLayout* layout = vtkGraphLayout::New();
  layout->SetInputConnection(numeric->GetOutputPort());

  vtkGraphLayoutView* view = vtkGraphLayoutView::New();
  //view->AddRepresentationFromInputConnection( tree->GetProducerPort() );
  view->AddRepresentationFromInput( tree );

  // Specify that we want to use the tree layout strategy.
  vtkTreeLayoutStrategy* strategy = vtkTreeLayoutStrategy::New();
  strategy->RadialOn();              // Radial layout (as opposed to standard top-down layout)
  strategy->SetAngle(360.0);         // The tree fills a full circular arc.
  layout->SetLayoutStrategy(strategy);

  view->SetLayoutStrategyToTree();
  view->SetVertexLabelArrayName( "Label" );
  view->SetVertexLabelVisibility(true);

  view->DisplayHoverTextOn();
  
  if ( colors ) {
	  //view->GetRenderer()->GetActors()->GetNextItem()->GetMapper()->SetLookupTable( lut );
	  vtkViewTheme* theme = vtkViewTheme::CreateMellowTheme();
	  //theme->SetPointLookupTable( lut );
	  //theme->SetPointColor(0,1,0);
	  theme->SetLineWidth(2);
	  view->ApplyViewTheme(theme);
	  theme->Delete();

	  view->ColorVerticesOn();
	  //view->SetVertexColorArrayName("VertexDegree");
	  view->SetVertexColorArrayName("NodeColor");
  }

  /*
  // vtkGraphToPolyData converts a graph or tree to polydata.
  vtkGraphToPolyData* graphToPoly = vtkGraphToPolyData::New();
  graphToPoly->SetInputConnection(layout->GetOutputPort());

  // Create the standard VTK polydata mapper and actor
  // for the connections (edges) in the tree.
  vtkPolyDataMapper* edgeMapper = vtkPolyDataMapper::New();
  edgeMapper->SetInputConnection(graphToPoly->GetOutputPort());
  vtkActor* edgeActor = vtkActor::New();
  edgeActor->SetMapper(edgeMapper);
  edgeActor->GetProperty()->SetColor(0.0, 0.5, 1.0);
  
  // Glyph the points of the tree polydata to create
  // VTK_VERTEX cells at each vertex in the tree.
  vtkGlyph3D* vertGlyph = vtkGlyph3D::New();
  vertGlyph->SetInputConnection(0, graphToPoly->GetOutputPort());
  vtkGlyphSource2D* glyphSource = vtkGlyphSource2D::New();
  glyphSource->SetGlyphTypeToVertex();
  vertGlyph->SetInputConnection(1, glyphSource->GetOutputPort());
  
  // Create a mapper for the vertices, and tell the mapper
  // to use the specified color array.
  vtkPolyDataMapper* vertMapper = vtkPolyDataMapper::New();
  vertMapper->SetInputConnection(vertGlyph->GetOutputPort());
    
  // Create an actor for the vertices.  Move the actor forward
  // in the z direction so it is drawn on top of the edge actor.
  vtkActor* vertActor = vtkActor::New();
  vertActor->SetMapper(vertMapper);
  vertActor->GetProperty()->SetPointSize(5);
  vertActor->SetPosition(0, 0, 0.001);
  
  // Use a dynamic label mapper to draw the labels.  This mapper
  // does not allow labels to overlap, as long as the camera is
  // not rotated from pointing down the z axis.
  vtkDynamic2DLabelMapper* labelMapper = vtkDynamic2DLabelMapper::New();
  labelMapper->SetInputConnection(graphToPoly->GetOutputPort());
  labelMapper->GetLabelTextProperty()->SetJustificationToLeft();
  labelMapper->GetLabelTextProperty()->SetColor(0, 0, 0);
  //if (labelArray)
    //{
    //labelMapper->SetLabelModeToLabelFieldData();
    //labelMapper->SetFieldDataName(labels);
    //}
  //vtkActor2D* labelActor = vtkActor2D::New();
  //labelActor->SetMapper(labelMapper);
  
  // Add the edges, vertices, and labels to the renderer.
  vtkRenderer* ren = vtkRenderer::New();
  ren->SetBackground(0.8, 0.8, 0.8);
  ren->AddActor(edgeActor);
  ren->AddActor(vertActor);
  //ren->AddActor(labelActor);
  
  // Setup the render window and interactor.
  vtkRenderWindow* win = vtkRenderWindow::New();
  win->AddRenderer(ren);
  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(win);
  
  // Constrain movement to zoom and pan using the image interactor style.
  vtkInteractorStyleImage* style = vtkInteractorStyleImage::New();
  iren->SetInteractorStyle(style);
  
  // Start the main application loop.
  iren->Initialize();
  iren->Start();
  
  // Clean up.
  style->Delete();
  iren->Delete();
  win->Delete();
  ren->Delete();
  //labelActor->Delete();
  labelMapper->Delete();
  vertActor->Delete();
  vertMapper->Delete();
  glyphSource->Delete();
  vertGlyph->Delete();
  edgeMapper->Delete();
  edgeActor->Delete();
  graphToPoly->Delete();
  strategy->Delete();
  layout->Delete();
  numeric->Delete();
  */

  view->GetRenderWindow()->SetSize(1920, 1080);
  view->ResetCamera();
  //view->GetRenderer()->RotateY(90);
  view->Render();
  view->GetInteractor()->Start();

 
  return 0;
}

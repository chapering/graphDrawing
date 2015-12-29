
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



#include "vtkCosmicTreeLayoutStrategy.h"
#include "vtkDataRepresentation.h"
#include "vtkHierarchicalGraphView.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRegressionTestImage.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderedHierarchyRepresentation.h"
#include "vtkSelection.h"
#include "vtkSplineGraphEdges.h"
#include "vtkTestUtilities.h"
#include "vtkViewTheme.h"
#include "vtkXMLTreeReader.h"

#include "vtkIdTypeArray.h"
#include "vtkStringToNumeric.h"

#include "vtkTreeRingView.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPoints.h"
#include "vtkLookupTable.h"

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

  //tree->print(std::cout);
  /*
  vtkUnsignedCharArray* colors = NULL;
  if ( argc >= 3 ) {
	  colors = vtkUnsignedCharArray::New();
	  colors->SetNumberOfComponents(3);
	  ifstream ifs (argv[2]);
	  if (!ifs.is_open()) {
		  cerr << "can not load color values from " << argv[2] << ", ignored.\n";
		  return -1;
	  }

	  double vals[3];
	  int ilabTotal = 0;
	  while (ifs) {
		  ifs >> vals[0] >> vals[1] >> vals[2];
		  colors->InsertNextTuple3( int(vals[0]*255), int(vals[1]*255), int(vals[2]*255) );
		  ilabTotal ++;
	  }
	  ifs.close();

	  cerr << ilabTotal << " colors loaded from file: " << argv[2] << "\n";
  }
  */

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

  vtkStringArray* labels = vtkStringArray::New();
  labels->SetName("Label");
  for ( int i = 0; i < nodenum; i ++ )
  {
    if ( my_tree.get()->is_leaf( i ) )
    {
      labels->InsertValue(node[i], my_tree.get()->label(i) );
	  dist->InsertNextValue(my_tree.get()->length_to_parent(i));
    }
    else
    {
      string str = boost::lexical_cast<string>(i);
      //labels->InsertValue(node[i], str );
      labels->InsertValue(node[i], "");
	  dist->InsertNextValue(my_tree.get()->length_to_parent(i));
    }
  }

  graph->GetVertexData()->AddArray(labels);
  graph->GetEdgeData()->AddArray(dist);
 
  // add pedigree ids for vertices and edges
  // Create vertex pedigree ids
  vtkSmartPointer<vtkIdTypeArray> vids = vtkSmartPointer<vtkIdTypeArray>::New();
  vids->SetName("vertex pedigree id");
  vtkIdType numVerts = graph->GetNumberOfVertices();
  vids->SetNumberOfTuples(numVerts);
  for (vtkIdType i = 0; i < numVerts; ++i)
  {
	  vids->SetValue(i, i);
  }
  graph->GetVertexData()->SetPedigreeIds(vids);

  // Create vertex pedigree ids
  vtkSmartPointer<vtkIdTypeArray> eids = vtkSmartPointer<vtkIdTypeArray>::New();
  eids->SetName("edge pedigree id");
  vtkIdType numEdges = graph->GetNumberOfEdges();
  eids->SetNumberOfTuples(numEdges);
  for (vtkIdType i = 0; i < numEdges; ++i)
  {
	  eids->SetValue(i, i);
  }
  graph->GetEdgeData()->SetPedigreeIds(eids);

  vtkTree* tree = vtkTree::New();
  tree->CheckedShallowCopy(graph);

  int leavecnt = 1;
  for ( int i = 0; i < numVerts; i ++ )
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

  tree->GetVertexData()->SetPedigreeIds(vids);
  tree->GetEdgeData()->SetPedigreeIds(eids);

  vtkSmartPointer<vtkTreeRingView> view = vtkSmartPointer<vtkTreeRingView>::New();
  view->SetTreeFromInput(tree);
  view->SetGraphFromInput(graph);
  view->SetAreaColorArrayName("VertexDegree");
  /*
  view->SetAreaHoverArrayName("vertex pedigree id");
  view->SetAreaLabelArrayName("vertex pedigree id");
  */
  view->SetAreaHoverArrayName("Label");
  view->SetAreaLabelArrayName("Label");
  view->SetAreaLabelVisibility(true);
  view->SetShrinkPercentage(0.02);
  view->SetBundlingStrength(.5);
  view->Update();
  view->SetEdgeColorArrayName("edge pedigree id");
  //view->SetColorEdges(true);
  view->RootAtCenterOn();
  view->SetRootAngles (0,360);
  view->SetLayerThickness(10.5);
  view->ColorAreasOn();
  //view->ColorEdgesOn();
  //view->SetEdgeScalarBarVisibility(true);
  //view->EdgeLabelVisibilityOn();
  if ( colors ) {
	  tree->GetVertexData()->AddArray( colors );
	  tree->GetVertexData()->SetActiveScalars("NodeColor");
	  view->SetAreaColorArrayName("NodeColor");
  }
  view->EdgeLabelVisibilityOff();

  /*
  VTK_CREATE(vtkStringToNumeric, numeric);
  numeric->SetInput(tree);
  
  // Graph layout view
   VTK_CREATE(vtkHierarchicalGraphView, view);
  view->DisplayHoverTextOff();
  view->GetRenderWindow()->SetMultiSamples(0);
  view->SetHierarchyFromInputConnection(numeric->GetOutputPort());
  //view->SetGraphFromInputConnection(reader2->GetOutputPort());
  //view->SetVertexColorArrayName("VertexDegree");
  //view->SetColorVertices(true);
  view->SetVertexLabelArrayName("Label");
  view->SetVertexLabelVisibility(true);
  //view->SetScalingArrayName("TreeRadius");

  view->Update(); // Needed for now
  //view->SetGraphEdgeColorArrayName("graph edge");
  //view->SetColorGraphEdgesByArray(true);
  vtkRenderedHierarchyRepresentation::SafeDownCast(view->GetRepresentation())->SetGraphSplineType(vtkSplineGraphEdges::CUSTOM, 0);

  VTK_CREATE(vtkCosmicTreeLayoutStrategy, ct);
  //ct->SetNodeSizeArrayName("VertexDegree");
  //ct->SetSizeLeafNodesOnly(true);
  view->SetLayoutStrategy(ct);
  */

  
  // Apply a theme to the views
  vtkViewTheme* const theme = vtkViewTheme::CreateMellowTheme();
  theme->SetLineWidth(2);
  theme->SetSelectedPointColor ( 0,0,0);
  theme->SetSelectedCellColor( 0,0,1 );
  view->ApplyViewTheme(theme);
  theme->Delete();
 
  view->ResetCamera();

  // JC added this
  view->Render();
  view->GetInteractor()->Initialize();
  view->GetInteractor()->Start();
 
  return 0;
}

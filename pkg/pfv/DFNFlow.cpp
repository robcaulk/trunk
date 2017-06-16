 
/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#ifdef FLOW_ENGINE

#include<pkg/dem/JointedCohesiveFrictionalPM.hpp>
// #include<pkg/dem/ScGeom.hpp>

//keep this #ifdef for commited versions unless you really have stable version that should be compiled by default
//it will save compilation time for everyone else
//when you want it compiled, you can pass -DDFNFLOW to cmake, or just uncomment the following line

#define DFNFLOW

#ifdef DFNFLOW
#include "FlowEngine_DFNFlowEngineT.hpp"

class DFNCellInfo : public FlowCellInfo_DFNFlowEngineT
{
	public:
//	Real anotherVariable;
	bool crack;
	Real crackArea;// the volume of cracks
    	bool fractureTip; //  The cell is the fracture tip
	Real cellHalfWidth; // distance from cell to injection point
	int breakType = 2;

// 	bool preExistingJoint;
// 	void anotherFunction() {};
// 	DFNCellInfo() : FlowCellInfo(),crack(false)  {}
//	DFNCellInfo() : crack(false), crackArea(0) {}
};

class DFNVertexInfo : public FlowVertexInfo_DFNFlowEngineT {
	public:
	//same here if needed
};

typedef CGT::_Tesselation<CGT::TriangulationTypes<DFNVertexInfo,DFNCellInfo> > DFNTesselation;
#ifdef LINSOLV
class DFNBoundingSphere : public CGT::FlowBoundingSphereLinSolv<DFNTesselation>
#else
class DFNBoundingSphere : public CGT::FlowBoundingSphere<DFNTesselation>
#endif
{
public:
  void saveVtk(const char* folder)
  {
      Tesselation& Tes = T[noCache?(!currentTes):currentTes];
      RTriangulation& Tri = Tes.Triangulation();

        static unsigned int number = 0;
        char filename[80];
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(filename,"%s/out_%d.vtk",folder,number++);
	int firstReal=-1;

	//count fictious vertices and cells
	vtkInfiniteVertices=vtkInfiniteCells=0;
 	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (!isDrawable) vtkInfiniteCells+=1;
	}
	for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
                if (!v->info().isReal()) vtkInfiniteVertices+=1;
                else if (firstReal==-1) firstReal=vtkInfiniteVertices;
	}

        basicVTKwritter vtkfile((unsigned int) Tri.number_of_vertices()-vtkInfiniteVertices, (unsigned int) Tri.number_of_finite_cells()-vtkInfiniteCells);

        vtkfile.open(filename,"test");

        vtkfile.begin_vertices();
        double x,y,z;
        for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
		if (v->info().isReal()){
		x = (double)(v->point().point()[0]);
                y = (double)(v->point().point()[1]);
                z = (double)(v->point().point()[2]);
                vtkfile.write_point(x,y,z);}
        }
        vtkfile.end_vertices();

        vtkfile.begin_cells();
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
        	if (isDrawable){vtkfile.write_cell(cell->vertex(0)->info().id()-firstReal, cell->vertex(1)->info().id()-firstReal, cell->vertex(2)->info().id()-firstReal, cell->vertex(3)->info().id()-firstReal);}
        }
        vtkfile.end_cells();

	if (permeabilityMap){
	vtkfile.begin_data("Permeability",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
	 if (isDrawable){vtkfile.write_data((cell->info().kNorm()[0] + cell->info().kNorm()[1] +cell->info().kNorm()[2] +cell->info().kNorm()[3]) / 4 );}
        }  // rcaulk since .s was not working, trying to record the average permeability of each cell
	vtkfile.end_data();}
	else{
	vtkfile.begin_data("Pressure",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().p());}
	}
	vtkfile.end_data();}

	if (1){
	averageRelativeCellVelocity();
	vtkfile.begin_data("Velocity",CELL_DATA,VECTORS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().averageVelocity()[0],cell->info().averageVelocity()[1],cell->info().averageVelocity()[2]);}
	}
	vtkfile.end_data();}
	
	if(1){
	vtkfile.begin_data("fracturedCells",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().crack);}
	}
	vtkfile.end_data();}

	if(1){
	vtkfile.begin_data("fractureTip",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().fractureTip);}
	}
	vtkfile.end_data();}

	if (1){
	vtkfile.begin_data("volumeChange",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().dv());}
	}
	vtkfile.end_data();}

	if(1){
	vtkfile.begin_data("breakType",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().breakType);}
	}
	vtkfile.end_data();}

////	can we parallelize these?
//	if(1){
//	vtkfile.begin_data("fractureTip",CELL_DATA,SCALARS,FLOAT);
//	#ifdef YADE_OPENMP
//    	const long size = Tes.cellHandles.size();
//		#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
//    	for (long i=0; i<size; i++){
//			CellHandle& cell = Tes.cellHandles[i];
//        #else
//            FOREACH(CellHandle& newCell, Tes.cellHandles){
//       #endif
//		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
//		if (isDrawable){vtkfile.write_data(cell->info().fractureTip);}
//	}
//	vtkfile.end_data();}


  }
};

typedef TemplateFlowEngine_DFNFlowEngineT<DFNCellInfo,DFNVertexInfo, DFNTesselation,DFNBoundingSphere> DFNFlowEngineT;
REGISTER_SERIALIZABLE(DFNFlowEngineT);
YADE_PLUGIN((DFNFlowEngineT));
class DFNFlowEngine : public DFNFlowEngineT
{
	public :
	void trickPermeability(Solver* flow);
    void interpolateCrack(Tesselation& Tes, Tesselation& NewTes);
	void trickPermeability (RTriangulation::Facet_circulator& facet,Real fracturePerm, RTriangulation::Finite_edges_iterator& edge, Solver* flow);
	void trickPermeability (RTriangulation::Finite_edges_iterator& edge,Real fracturePerm, Solver* flow);
	void setPositionsBuffer(bool current);
	void checkOldTesselation(Tesselation& oldTes, Tesselation& newTes, const CellHandle& newCell1, const CellHandle& newCell2);
	Real stepCrackHalfWidth; // step crack halfwidth
	Real fractureHalfWidth = 0; // full fracture halfwidth
    	Real averageAperture;
	Real averageFracturePermeability;
    	Real maxAperture;
	Real crackArea;
	Real branchIntensity;
	int breakType;
	Real leakOffRate;
	Real getCrackHalfWidth() {return fractureHalfWidth;}
    	Real getAverageAperture() {return averageAperture;}
    	Real getMaxAperture() {return maxAperture;}
	Real getCrackArea() {return crackArea;}
	Real getBranchIntensity() {return branchIntensity;}
	Real getLeakOffRate() {return leakOffRate;}
	Point injectionCellCenter;
// 	void computeTotalFractureArea(Real totalFracureArea,bool printFractureTotalArea);/// Trying to get fracture's surface
//	Real totalFracureArea; /// Trying to get fracture's surface

//	CELL_SCALAR_GETTER(double,.crackArea,crackArea)

	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(DFNFlowEngine,DFNFlowEngineT,"This is an enhancement of the FlowEngine for intact and fractured rocks that takes into acount pre-existing discontinuities and bond breakage between particles. The local conductivity around the broken link is calculated according to parallel plates model",
	((Real, jointResidual, 0,,"calibration parameter for residual aperture of joints"))
	((Real, inducedResidual, 0,,"calibration parameter for residual aperture of induced cracks"))
	((bool, updatePositions, false,,"update particles positions when rebuilding the mesh (experimental)"))
 	((bool, printFractureTotalArea, 0,,"The final fracture area computed through the network")) /// Trying to get fracture's surface
	((bool, calcCrackArea, true,,"The amount of crack per pore () is updated if calcCrackArea=True")) /// Trying to get fracture's surface((bool))
    	((bool, calcCrackHalfWidth, true,,"Calculate the fracture halfwidth")) /// Trying to get fracture's surface
	,,,
//	.def("getCrackArea",&DFNFlowEngine::crackArea,(boost::python::arg("id")),"get the cracked area within cell 'id'.")
	.def("getCrackHalfWidth", &DFNFlowEngine::getCrackHalfWidth, "report the current fracture half width")
    	.def("getAverageAperture", &DFNFlowEngine::getAverageAperture, "report the current average aperture")
    	.def("getMaxAperture", &DFNFlowEngine::getMaxAperture, "report the max aperture")
	.def("getCrackArea", &DFNFlowEngine::getCrackArea, "report the crack area")
	.def("getBranchIntensity", &DFNFlowEngine::getBranchIntensity, "report branch intensity (fracture area/half length)")
	.def("getLeakOffRate", &DFNFlowEngine::getLeakOffRate, "report leak-off rate")
// 	.def("computeTotalFractureArea",&DFNFlowEngineT::computeTotalFractureArea," Compute and print the total fracture area of the network") /// Trying to get fracture's surface
// 	.def("trickPermeability",&DFNFlowEngineT::trickPermeability,"measure the mean trickPermeability in the period")
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(DFNFlowEngine);
YADE_PLUGIN((DFNFlowEngine));

//In this version, we never update positions when !updatePositions, i.e. keep triangulating the same positions
void DFNFlowEngine::setPositionsBuffer(bool current)
{
	vector<posData>& buffer = current? positionBufferCurrent : positionBufferParallel;
	if (!updatePositions && buffer.size()>0) return;
	buffer.clear();
	buffer.resize(scene->bodies->size());
	shared_ptr<Sphere> sph ( new Sphere );
        const int Sph_Index = sph->getClassIndexStatic();
	FOREACH ( const shared_ptr<Body>& b, *scene->bodies ) {
                if (!b || ignoredBody==b->getId()) continue;
                posData& dat = buffer[b->getId()];
		dat.id=b->getId();
		dat.pos=b->state->pos;
		dat.isSphere= (b->shape->getClassIndex() ==  Sph_Index);
		if (dat.isSphere) dat.radius = YADE_CAST<Sphere*>(b->shape.get())->radius;
		dat.exists=true;
	}
}

// function allows us to interpolate information about fractured/non fractured cells so we can identify newly fractured cells, monitor half width, and identify fracture tip
void DFNFlowEngine::interpolateCrack(Tesselation& Tes, Tesselation& NewTes){
        RTriangulation& Tri = Tes.Triangulation();
		RTriangulation& newTri = NewTes.Triangulation();
		FiniteCellsIterator cellEnd = newTri.finite_cells_end();
	#ifdef YADE_OPENMP
    	const long size = NewTes.cellHandles.size();
	#pragma omp parallel for num_threads(ompThreads>0 ? ompThreads : 1)
    	for (long i=0; i<size; i++){
		CellHandle& newCell = NewTes.cellHandles[i];
        #else
	FOREACH(CellHandle& newCell, NewTes.cellHandles){
        #endif
		CVector center (0,0,0);
		if (newCell->info().fictious()==0) for ( int k=0;k<4;k++ ) center= center + 0.25* (Tes.vertex(newCell->vertex(k)->info().id())->point()-CGAL::ORIGIN);		
		
		CellHandle oldCell = Tri.locate(Point(center[0],center[1],center[2]));
		newCell->info().crack = oldCell->info().crack;
		newCell->info().fractureTip = oldCell->info().fractureTip;
		newCell->info().cellHalfWidth = oldCell->info().cellHalfWidth;
//		newCell->info().id = oldCell->info().id;


		if (oldCell->info().crack && !oldCell->info().fictious()){
			Real facetFlowRate = 0;
			facetFlowRate -= oldCell->info().dv();
			for (int k=0; k<4;k++) {
				
				if (!oldCell->neighbor(k)->info().crack){
					facetFlowRate= oldCell->info().kNorm()[k]*(oldCell->info().shiftedP()-oldCell->neighbor(k)->info().shiftedP());
					leakOffRate += facetFlowRate;
				}
			}
		}
	}


    }

// function used to change fractureTip info to false for cells neighboring new fractureTip = true cells. 
void DFNFlowEngine::checkOldTesselation(Tesselation& oldTes, Tesselation& newTes, const CellHandle& newCell1, const CellHandle& newCell2){
	RTriangulation& Tri = oldTes.Triangulation();	
	CVector center1 (0,0,0);
	CVector center2 (0,0,0);

	if (newCell1->info().fictious()==0) for ( int k=0;k<4;k++ ) center1= center1 + 0.25* (oldTes.vertex(newCell1->vertex(k)->info().id())->point()-CGAL::ORIGIN);	
	CellHandle oldCell1 = Tri.locate(Point(center1[0],center1[1],center1[2]));

	if (newCell2->info().fictious()==0) for ( int k=0;k<4;k++ ) center2= center2 + 0.25* (oldTes.vertex(newCell2->vertex(k)->info().id())->point()-CGAL::ORIGIN);	
	CellHandle oldCell2 = Tri.locate(Point(center2[0],center2[1],center2[2]));

       	for (int i = 0; i<=3; i++) { 
		CellHandle oldNeighborCell = oldCell1->neighbor(i); 
		if (oldNeighborCell->info().fractureTip==1) {
		CVector neighborCenter (0,0,0);
		for ( int k=0;k<4;k++ ) neighborCenter= neighborCenter + 0.25* (newTes.vertex(oldNeighborCell->vertex(k)->info().id())->point()-CGAL::ORIGIN);
		CellHandle newNeighborCell = newTes.Triangulation().locate(Point(neighborCenter[0], neighborCenter[1], neighborCenter[2]));
	newNeighborCell->info().fractureTip=0;
		}
	}

for (int i = 0; i<=3; i++) { 
		CellHandle oldNeighborCell = oldCell2->neighbor(i); 
		if (oldNeighborCell->info().fractureTip==1) {
		CVector neighborCenter (0,0,0);
		for ( int k=0;k<4;k++ ) neighborCenter= neighborCenter + 0.25* (newTes.vertex(oldNeighborCell->vertex(k)->info().id())->point()-CGAL::ORIGIN);
		CellHandle newNeighborCell = newTes.Triangulation().locate(Point(neighborCenter[0], neighborCenter[1], neighborCenter[2]));
	newNeighborCell->info().fractureTip=0;

		}
	}
}


void DFNFlowEngine::trickPermeability(RTriangulation::Facet_circulator& facet, Real fracturePerm, RTriangulation::Finite_edges_iterator& ed_it, Solver* flow)
{

	const RTriangulation::Facet& currentFacet = *facet; // seems verbose but facet->first was declaring a junk cell and crashing program (https://bugs.launchpad.net/yade/+bug/1666339)
	const RTriangulation& Tri = flow->T[flow->currentTes].Triangulation();
	const CellHandle& cell1 = currentFacet.first;
	const CellHandle& cell2 = currentFacet.first->neighbor(facet->second);
	if ( Tri.is_infinite(cell1) || Tri.is_infinite(cell2)) cerr<<"Infinite cell found in trickPermeability, should be handled somehow, maybe"<<endl;
	cell1->info().kNorm()[currentFacet.second]=cell2->info().kNorm()[Tri.mirror_index(cell1, currentFacet.second)] = fracturePerm;
	if (first) {
		cell1->info().crack= 1;
		cell2->info().crack= 1;
	}
	cell1->info().breakType=cell2->info().breakType=breakType;
	//For vtk recorder:
    // if the fractured cell is newly fractured, must be at tip and therefore gives info about the fracture dimensions and stresses
	if (!first && !cell1->info().crack && !cell1->info().isFictious && calcCrackHalfWidth && !flow->imposedF.empty()) {        
	
	Point& cellCenter = cell1->info();
        CVector halfWidthVector = cellCenter - injectionCellCenter;
        Real halfWidth = sqrt(halfWidthVector.squared_length());
        cell1->info().cellHalfWidth = halfWidth;
		if (halfWidth > fractureHalfWidth){
			stepCrackHalfWidth = halfWidth;
//			cell1->info().fractureTip = 1;
		}
	checkOldTesselation(solver->T[solver->currentTes], flow->T[flow->currentTes], cell1, cell2);
//	cout<<"DFN ----- checked old tesselation" <<endl;	
	cell1->info().fractureTip = cell2->info().fractureTip =1;	
	cell1->info().crack = cell2->info().crack= 1;
	
        }
	cell2->info().blocked = cell1->info().blocked = cell2->info().Pcondition = cell1->info().Pcondition = false;//those ones will be included in the flow problem
	Point& CellCentre1 = cell1->info(); /// Trying to get fracture's surface 
	Point& CellCentre2 = cell2->info(); /// Trying to get fracture's surface 
//	CVector networkFractureLength = CellCentre1 - CellCentre2; /// Trying to get fracture's surface
//	double networkFractureDistance = sqrt(networkFractureLength.squared_length()); /// Trying to get fracture's surface
//	Real networkFractureArea = pow(networkFractureDistance,2);  /// Trying to get fracture's surface
//	totalFracureArea += networkFractureArea; /// Trying to get fracture's surface
// 	cout <<" ------------------ The total surface area up to here is --------------------" << totalFracureArea << endl;
// 	printFractureTotalArea = totalFracureArea; /// Trying to get fracture's surface 
	if (calcCrackArea) {
			CVector edge = ed_it->first->vertex(ed_it->second)->point().point() - ed_it->first->vertex(ed_it->third)->point().point();
			CVector unitV = edge*(1./sqrt(edge.squared_length()));
			Point p3 = ed_it->first->vertex(ed_it->third)->point().point() + unitV*(cell1->info() - ed_it->first->vertex(ed_it->third)->point().point())*unitV;
			Real halfCrackArea = 0.25*sqrt(std::abs(cross_product(CellCentre1-p3,CellCentre2-p3).squared_length()));//
			cell1->info().crackArea += halfCrackArea;
			cell2->info().crackArea += halfCrackArea;
			crackArea += 2*halfCrackArea;
		}
}

void DFNFlowEngine::trickPermeability(RTriangulation::Finite_edges_iterator& edge, Real fracturePerm,Solver* flow)
{
	const RTriangulation& Tri = flow->T[flow->currentTes].Triangulation();
//	cout << "DFN ---- defined the Triangulation to be used in facet circulator" << endl;
	
	RTriangulation::Facet_circulator facet1 = Tri.incident_facets(*edge);
	RTriangulation::Facet_circulator facet0=facet1++;
	trickPermeability(facet0, fracturePerm, edge, flow);
	while ( facet1!=facet0 ) {trickPermeability(facet1, fracturePerm, edge, flow); facet1++;}
	/// Needs the fracture surface for this edge?
// 	double edgeArea = solver->T[solver->currentTes].computeVFacetArea(edge); cout<<"edge area="<<edgeArea<<endl;
}

void DFNFlowEngine::trickPermeability(Solver* flow)
{	
	leakOffRate = 0;	
	const RTriangulation& Tri = flow->T[flow->currentTes].Triangulation();
//	cout << "DFN --- Assigned tri" <<endl;
    if (!first) interpolateCrack(solver->T[solver->currentTes], flow->T[flow->currentTes]);
//	cout << "DFN ---- interpolated crack" << endl;
    const JCFpmPhys* jcfpmphys;
	const shared_ptr<InteractionContainer> interactions = scene->interactions;
	int numberOfCrackedOrJoinedInteractions = 0;
	int numberOfTrickedCells = 0;
	Real sumOfPermeability = 0;
	Real SumOfApertures = 0.;
	averageAperture = 0;
	crackArea = 0;
	branchIntensity = 0;
	maxAperture = 0;
	
	
//	Real totalFracureArea=0; /// Trying to get fracture's surface
	if (!flow->imposedF.empty()) injectionCellCenter = flow->IFCells[0]->info();
//	cout << "DFN ---- defined injection Cell center" << endl;
	stepCrackHalfWidth = 0; // The fracture half width
// 	const shared_ptr<IGeom>& ig;
// 	const ScGeom* geom; // = static_cast<ScGeom*>(ig.get());
	FiniteEdgesIterator edge = Tri.finite_edges_begin();
	for( ; edge!= Tri.finite_edges_end(); ++edge) {

		const VertexInfo& vi1=(edge->first)->vertex(edge->second)->info();
		const VertexInfo& vi2=(edge->first)->vertex(edge->third)->info();
		const shared_ptr<Interaction>& interaction=interactions->find( vi1.id(),vi2.id() );
		
		if (interaction && interaction->isReal()) {
			if (edge->first->info().isFictious) continue; // avoid trick permeability for fictitious cells
			jcfpmphys = YADE_CAST<JCFpmPhys*>(interaction->phys.get());
			
			if ( jcfpmphys->isOnJoint || jcfpmphys->isBroken ) {
				
				breakType = jcfpmphys->breakType;
				numberOfCrackedOrJoinedInteractions +=1;
				
				// here are some workarounds
// 				Real residualAperture = jcfpmphys->isOnJoint? jointResidual : 0.1*jointResidual;
				Real residualAperture = jcfpmphys->isOnJoint? jointResidual : inducedResidual; // do we define a residual aperture for induced cracks?
// 				Real residualAperture = jointResidual;
// 				cout<<"residual aperture = " << residualAperture <<endl;
				
				if (!jcfpmphys->isOnJoint && jcfpmphys->crackJointAperture<=0) continue; // avoid trick permeability for closed induced crack (permeability=matrix permeability)
				
				Real aperture = jcfpmphys->crackJointAperture;
				Real fracturePerm = pow((aperture+residualAperture),3)/(12*viscosity);
				numberOfTrickedCells += 1;
				sumOfPermeability += fracturePerm;
// 				cout<<"aperture = " << aperture <<endl;
                	if (aperture > maxAperture) maxAperture=aperture;
				SumOfApertures += aperture;
				trickPermeability(edge, fracturePerm, flow);
// 				trickPermeability(edge, jcfpmphys->crackJointAperture, residualAperture); // we should be able to use this line
				
			};	
		}
	}
	averageAperture = SumOfApertures/numberOfCrackedOrJoinedInteractions; /// DEBUG
	averageFracturePermeability = sumOfPermeability/numberOfTrickedCells;
	if (debug) cout << "DFN --- Average Fracture Permeability =" << averageFracturePermeability << endl;

	if (stepCrackHalfWidth > fractureHalfWidth) fractureHalfWidth = stepCrackHalfWidth;
	branchIntensity = crackArea/fractureHalfWidth;
//	cout << "fracture half width = " << fractureHalfWidth << endl;
// 	cout << " Average aperture in joint ( -D ) = " << AverageAperture << endl; /// DEBUG
}

// Real DFNFlowEngine::computeTotalFractureArea(totalFracureArea,printFractureTotalArea) /// Trying to get fracture's surface
// {
// 	 if (printFractureTotalArea >0) {
// 		cout<< " The total fracture area computed from the Network is: " << totalFracureArea <<endl;
// 	 }
// }

#endif //DFNFLOW
#endif //FLOWENGINE

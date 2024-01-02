// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <memory>
#include <iomanip>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/common/initialize.hh>
#include <dumux/io/grid/gridmanager_foam.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <test/geometry/writetriangulation.hh>

#include <dumux/microcirculation/rbcnetwork.hh>
#include <dumux/microcirculation/geometry/rbcentityset.hh>

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI+X
    Dumux::initialize(argc, argv);

    // initialize parameters
    Parameters::init([](auto& params){
        params["Grid.File"] = "../../examples/grids/honeycomb.dgf";
        params["Grid.Refinement"] = "3";
    });

    // construct a grid from file
    using Grid = Dune::FoamGrid<1, 3>;
    GridManager<Grid> gridManager; gridManager.init();
    const auto leafGridView = gridManager.grid().leafGridView();
    const auto level0GridView = gridManager.grid().levelGridView(0);
    const auto gridData = gridManager.getGridData();

    Dune::VTKWriter<Grid::LeafGridView> vtkWriter(leafGridView);
    vtkWriter.write("grid");

    // construct a bounding box tree of the grid
    using GridEntitySet = GridViewGeometricEntitySet<Grid::LeafGridView, 0>;
    auto gridBBoxTree = std::make_shared<BoundingBoxTree<GridEntitySet>>(std::make_shared<GridEntitySet>(leafGridView));

    // construct an rbc network on this grid
    std::vector<double> radii(level0GridView.size(0));
    for (const auto& element : elements(level0GridView))
        radii[leafGridView.indexSet().index(element)] = gridData->parameters(element)[0];
    auto rbcNetwork = std::make_shared<RBCNetwork<Grid::LevelGridView>>(level0GridView, radii);

    // construct a bounding box tree of the rbcs
    using RBCEntitySet = Dumux::RBCEntitySet<RBCNetwork<Grid::LevelGridView>>;
    auto rbcBBoxTree = std::make_shared<BoundingBoxTree<RBCEntitySet>>(std::make_shared<RBCEntitySet>(rbcNetwork));

    // intersect the bounding box trees
    Dune::Timer timer;
    const auto treeIntersections = intersectingEntities(*gridBBoxTree, *rbcBBoxTree);
    std::cout << "Computed " << treeIntersections.size() << " tree intersections in " << timer.elapsed() << std::endl;
    timer.reset();

    // write out all the intersections
    // convert format
    using Point = Dune::FieldVector<double, 3>;
    std::vector<std::vector<Point>> intersections;
    intersections.reserve(treeIntersections.size());
    for (const auto& is : treeIntersections)
        intersections.emplace_back(std::vector<Point>(is.corners()));
    std::cout << "Converted to output format in " << timer.elapsed() << " seconds." << std::endl;
    timer.reset();

    std::cout << "Writing " << intersections.size() << " intersections to file ...";
    Dumux::writeVTKPolyDataLines(intersections, "intersections");
    std::cout << " done ( " << timer.elapsed() << " seconds)." << std::endl;

    return 0;
}
catch (const Dumux::ParameterException &e) {
    std::cerr << e << ". Abort!\n";
    return 1;
}
catch (const Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 3;
}

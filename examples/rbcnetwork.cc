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

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/io/grid/gridmanager_foam.hh>
#include <dumux/io/pointcloudvtkwriter.hh>

#include <dumux/microcirculation/rbcnetwork.hh>
#include <dumux/microcirculation/rbcnetworkvtkwriter.hh>
#include <dumux/microcirculation/flowsolver.hh>

#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI+X
    Dumux::initialize(argc, argv);

    // read parameters from input file
    Parameters::init(argc, argv);

    // construct a grid from file
    using Grid = Dune::FoamGrid<1, 3>;
    GridManager<Grid> gridManager; gridManager.init();
    const auto leafGridView = gridManager.grid().leafGridView();
    const auto gridData = gridManager.getGridData();

    // create vector of radii
    std::vector<double> radii(leafGridView.size(0));
    for (const auto& element : elements(leafGridView))
        radii[leafGridView.indexSet().index(element)] = gridData->parameters(element)[0];

    // construct an RBC network
    auto rbcNetwork = std::make_shared<RBCNetwork<Grid::LeafGridView>>(leafGridView, radii);

    // write the initial positions to file
    RBCNetworkPVDWriter rbcWriter(rbcNetwork, "rbcs");
    rbcWriter.write(0.0);

    // the grid geometry
    auto gridGeometry = std::make_shared<FlowSolver::GridGeometry>(leafGridView);

    // create the flow solver
    FlowSolver flowSolver(gridGeometry, *gridData);

    // update effective resistance and solve for flow field
    flowSolver.updateResistance(*rbcNetwork);
    flowSolver.solve(0.0);

    // the main time loop
    const double timeStep = getParam<double>("TimeLoop.DtInitial", 1e-3);
    const double tEnd = getParam<double>("TimeLoop.TEnd", 4.0);
    TimeLoop<double> timeLoop(0.0, timeStep, tEnd);
    timeLoop.start(); do
    {
        Dune::Timer rbcTimer;
        rbcNetwork->advance(timeStep, *gridGeometry, flowSolver.volumeFluxes());
        timeLoop.advanceTimeStep();
        std::cout << "RBC move took " << rbcTimer.elapsed() << " seconds." << std::endl;

        // write vtk
        rbcWriter.write(timeLoop.time());

        // flow solver
        flowSolver.updateResistance(*rbcNetwork);
        flowSolver.solve(timeLoop.time());

        // report statistics of this time step
        timeLoop.reportTimeStep();

    } while (!timeLoop.finished());

    timeLoop.finalize(leafGridView.comm());

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

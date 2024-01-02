// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Microcirculation
 * \author Timo Koch
 * \brief the flow solver
 */
#ifndef DUMUX_RBC_FLOWSOLVER_HH
#define DUMUX_RBC_FLOWSOLVER_HH

#include <iostream>
#include <memory>

#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/io/grid/gridmanager_foam.hh>

#include "bloodflowproblem.hh"
#include "bloodflowspatialparams.hh"

namespace Dumux {

//! Set the properties of the flow solver
namespace Properties {

namespace TTag { struct BloodFlow { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; }; }

template<class TypeTag>
struct Grid<TypeTag, TTag::BloodFlow> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::BloodFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct SolutionDependentAdvection<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentMolecularDiffusion<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };
template<class TypeTag>
struct SolutionDependentHeatConduction<TypeTag, TTag::BloodFlow> { static constexpr bool value = false; };

template<class TypeTag>
struct Problem<TypeTag, TTag::BloodFlow>
{ using type = BloodFlowProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BloodFlow>
{ using type = BloodFlowSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>, double>; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::BloodFlow>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BloodFlow>
{ using type = FluidSystems::OnePLiquid<double, Components::Constant<1, double> >; };

} // end namespace Properties

class FlowSolver
{
    using TypeTag = Properties::TTag::BloodFlow;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VtkOutputFields = GetPropType<TypeTag, Properties::IOFields>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
public:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
private:
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using LinearPDESolver = Dumux::LinearPDESolver<Assembler, LinearSolver>;
public:
    //! construct a new flow solver
    template<class GridData>
    FlowSolver(std::shared_ptr<const GridGeometry> gridGeometry,
               const GridData& gridData)
    : gridGeometry_(std::move(gridGeometry))
    {
        initialize_(gridData);
    }

    //! construct a new flow solver with a grid view
    template<class GridData>
    FlowSolver(const typename GridGeometry::GridView& gridView,
               const GridData& gridData)
    {
        gridGeometry_ = std::make_shared<GridGeometry>(gridView);
        initialize_(gridData);
    }

    //! update the resistance according to the current RBC distribution
    template<class RBCNetwork>
    void updateResistance(const RBCNetwork& rbcs)
    {
        // the resistance is treated as a spatial parameter for the flow problem
        problem_->spatialParams().updateResistance(rbcs);
        // update the transmissibilities
        gridVariables_->update(sol_, true);
    }

    //! solve for the pressure distribution and upate the volume fluxes
    void solve(double simulationTime)
    {
        // solve
        solver_->solve(sol_);

        // update volume fluxes
        computeVolumeFluxes_();

        // write output
        if (enableVtkOutput)
            outputModule_->write(simulationTime);
    }

    //! get the current volume fluxes
    const std::vector<double>& volumeFluxes() const
    { return volumeFluxes_; }

    //! get the grid geometry
    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

private:
    /*!
     * \brief Initialize all simulation components
     */
    template<class GridData>
    void initialize_(const GridData& gridData)
    {
        Dune::Timer timer;

        // the problem class
        problem_ = std::make_shared<Problem>(gridGeometry_, gridData);
        problem_->applyInitialSolution(sol_);

        // the volume and flux variables
        gridVariables_ = std::make_shared<GridVariables>(problem_, gridGeometry_);
        gridVariables_->init(sol_);

        // initialize the volume fluxes
        volumeFluxes_.resize(gridGeometry_->gridView().size(0), 0.0);

        // initialize vtk output
        enableVtkOutput = getParam<bool>("Problem.EnableVtkOutput", true);
        outputModule_ = std::make_shared<VtkOutputModule<GridVariables, SolutionVector>>(*gridVariables_, sol_, problem_->name());
        VtkOutputFields::initOutputModule(*outputModule_);
        problem_->addVtkOutputFields(*outputModule_);
        outputModule_->addField(volumeFluxes_, "Q (m^3/s)");

        // assembler & linear solver
        auto assembler = std::make_shared<Assembler>(problem_, gridGeometry_, gridVariables_);
        auto linearSolver = std::make_shared<LinearSolver>(gridGeometry_->gridView(), gridGeometry_->dofMapper());
        solver_ = std::make_shared<LinearPDESolver>(assembler, linearSolver);

        std::cout << "Flow solver setup took " << timer.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Compute the volume fluxes, averaged for each segment
     */
    void computeVolumeFluxes_()
    {
        volumeFluxes_.assign(gridGeometry_->gridView().size(0), 0.0);
        auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };

        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            const auto eIdx = gridGeometry_->gridView().indexSet().index(element);

            auto fvGeometry = localView(*gridGeometry_);
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVariables_->curGridVolVars());
            elemVolVars.bind(element, fvGeometry, sol_);

            auto elemFluxVars = localView(gridVariables_->gridFluxVarsCache());
            elemFluxVars.bind(element, fvGeometry, elemVolVars);

            const auto geo = element.geometry();
            auto orientation = geo.corner(1)-geo.corner(0);
            orientation /= orientation.two_norm();

            for (const auto& scvf : scvfs(fvGeometry))
            {
                // sign convection positive fluxes from 0->1
                const double sign = (scvf.unitOuterNormal()*orientation > 0.0) ? 1.0 : -1.0;
                if (!scvf.boundary())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                    volumeFluxes_[eIdx] += 0.5*sign*fluxVars.advectiveFlux(0, upwindTerm);
                }
                else
                {
                    const auto bcTypes = problem_->boundaryTypes(element, scvf);
                    if (bcTypes.hasOnlyDirichlet()) // Dirichlet
                    {
                        FluxVariables fluxVars;
                        fluxVars.init(*problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                        volumeFluxes_[eIdx] += 0.5*sign*fluxVars.advectiveFlux(0, upwindTerm);
                    }

                    else // Neumann
                    {
                        volumeFluxes_[eIdx] += 0.5*sign*problem_->neumann(element, fvGeometry, elemVolVars, elemFluxVars, scvf)[0]
                                                  * scvf.area() * elemVolVars[scvf.insideScvIdx()].extrusionFactor()
                                                  / elemVolVars[scvf.insideScvIdx()].density(); // volume flux from mass flux
                    }
                }
            }
        }
    }

    std::shared_ptr<const GridGeometry> gridGeometry_;
    std::shared_ptr<Problem> problem_;
    std::shared_ptr<GridVariables> gridVariables_;
    SolutionVector sol_;

    bool enableVtkOutput;
    std::shared_ptr<VtkOutputModule<GridVariables, SolutionVector>> outputModule_;

    std::shared_ptr<LinearPDESolver> solver_;
    std::vector<double> volumeFluxes_;
};

} // end namespace Dumux

#endif

// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Microcirculation
 * \brief The spatial parameters class blood flow problem
 */
#ifndef DUMUX_BlOOD_FLOW_SPATIALPARAMS_HH
#define DUMUX_BlOOD_FLOW_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/microcirculation/effectiveresistance.hh>
#include <dumux/microcirculation/redbloodcell.hh>

namespace Dumux {

/*!
 * \ingroup Microcirculation
 * \brief Definition of the spatial parameters for the blood flow problem
 */
template<class GridGeometry, class Scalar>
class BloodFlowSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, BloodFlowSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = BloodFlowSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    BloodFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here makes extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto radius = this->radius(eIdx);
        return M_PI*radius*radius;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperatureAtPos(const GlobalPosition&) const
    { return 273.15 + 37.0; } // Body temperature

    /*!
     * \brief the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const Scalar r = radius(scv.dofIndex());
        const Scalar muEff = relEffViscosity(scv.dofIndex());
        return r*r/8.0/muEff;
    }

    /*!
     * \brief return the radius of the circular pipe for the current sub-control volume in [m].
     * \param the index of the element
     */
    Scalar radius(std::size_t eIdxGlobal) const
    { return radius_[eIdxGlobal];}

    /*!
     * \brief the velocity estimate.
     * \param the index of the element
     */
    Scalar velocityEstimate(std::size_t eIdxGlobal) const
    { return velocityEstimate_[eIdxGlobal]; }

    /*!
     * \brief the relative effective viscosity
     * \param the index of the element
     */
    Scalar relEffViscosity(std::size_t eIdxGlobal) const
    { return relEffViscosity_[eIdxGlobal]; }

    /*!
     * \brief returns the porosity \f$[-]\f$
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! get the radii vector for e.g. output
    const std::vector<Scalar>& getRadii() const
    { return radius_; }

    //! get the rel. eff. viscosity vector for e.g. output
    const std::vector<Scalar>& getRelEffViscosity() const
    { return relEffViscosity_; }

    //! get the tube hematocrit vector for e.g. output
    const std::vector<Scalar>& getTubeHematocrit() const
    { return tubeHematocrit_; }

    //! read params from dgf
    template<class GridData>
    void readGridParams(const GridData& gridData)
    {
        const auto& gg = this->gridGeometry();
        auto numElements = gg.gridView().size(0);
        radius_.resize(numElements);
        velocityEstimate_.resize(numElements);
        relEffViscosity_.resize(numElements, 1.0);
        tubeHematocrit_.resize(numElements, 0.2);

        // gridview is a leafGridView. Parameters are only set on level 0.
        // elements have to inherit spatial parameters from their father.
        for (const auto& element : elements(gg.gridView()))
        {
            auto level0element = element;
            while(level0element.hasFather())
                level0element = level0element.father();

            auto eIdx = gg.elementMapper().index(element);
            radius_[eIdx] = gridData.parameters(level0element)[0];
            velocityEstimate_[eIdx] = gridData.parameters(level0element)[1];
        }
    }

    //! update the relative effective viscosity according to RBC distribution
    template<class RBCNetwork>
    void updateResistance(const RBCNetwork& rbcs)
    {
        const auto& gg = this->gridGeometry();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.gridView().indexSet().index(element);
            const auto& segmentRBCs = rbcs.segmentRBCs(eIdx);

            const auto r = radius_[eIdx];
            const auto volume = element.geometry().volume()*M_PI*r*r;
            tubeHematocrit_[eIdx] = segmentRBCs.size()*RedBloodCell::volume()/volume;
            const auto dischargeHematocrit = dischargeHematocritInVitro(r, tubeHematocrit_[eIdx]);
            relEffViscosity_[eIdx] = relEffViscosityInVitro(r, dischargeHematocrit);
        }
    }

private:
    std::vector<Scalar> radius_;
    std::vector<Scalar> velocityEstimate_;
    std::vector<Scalar> relEffViscosity_;
    std::vector<Scalar> tubeHematocrit_;
};

} // end namespace Dumux

#endif

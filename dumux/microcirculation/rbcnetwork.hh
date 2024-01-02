// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup Microcirculation
 * \brief A red blood cell network (collection of RBCs living on a grid)
 * \author Timo Koch
 */

#ifndef DUMUX_MICROCIRCULATION_RBC_NETWORK_HH
#define DUMUX_MICROCIRCULATION_RBC_NETWORK_HH

#include <cmath>
#include <vector>
#include <random>
#include <deque>

#include <dune/common/iteratorrange.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/iteratorrange.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

#include <dumux/microcirculation/redbloodcell.hh>
#include <dumux/microcirculation/rbciterator.hh>
#include <dumux/microcirculation/rbcmove.hh>

namespace Dumux {

/**
 * \file
 * \ingroup Microcirculation
 * \brief A network of red blood cells particles associated with a DUNE grid
 * \author Timo Koch
 * \tparam GridView a grid view on the dune grid
 *
 * The `RBCNetwork` class represents a network of red blood cell (RBC) particles associated with a DUNE grid.
 * It provides methods to initialize the network, compute the initial distribution of RBCs, and manage the RBC positions and radii.
 *
 * Usage:
 * - Create an instance of `RBCNetwork` by providing a grid view and a vector of radii for each element.
 * - Call the appropriate methods to initialize the network and compute the initial distribution of RBCs.
 * - Access and modify the RBC positions and radii using the provided methods.
 */
template<class GridView>
class RBCNetwork
{
    using GlobalPosition = RedBloodCell::GlobalPosition;

public:
    //! the tube network is given by a gridView and the discrete radius field (radius for each element)
    RBCNetwork(const GridView& gridView, const std::vector<double>& radii)
    : gridView_(gridView)
    , radii_(radii)
    {
        // compute an initial distribution of RBCs
        // make sure we roughly satisfy an initially given tube hematocrit
        // default: Schmid (2017) 10.2.1
        const double initTubeHematocrit = getParam<double>("RBCModel.InitTubeHematocrit", 0.2);
        if (initTubeHematocrit > 1.0)
            DUNE_THROW(ParameterException, "RBCModel.InitTubeHematocrit is larger than 1.0!");

        // create a random number generator
        std::mt19937 gen{std::random_device{}()};
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        // create id counter
        std::size_t rbcIdCounter = 0;

        // create orientation vectors and length
        orientations_.resize(gridView.size(0));
        vesselLength_.resize(gridView.size(0));

        // precompute node neighbors
        nodeNeighbors_.resize(gridView.size(1));

        // read thickness of the endothelial surface layer
        eslThickness_ = getParam<double>("VesselWall.EndothelialSurfaceLayerThickness", 0.0);

        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(element);
            const auto vesselRadius = radii_[eIdx];
            const auto rbcRadius = vesselRadius-eslThickness_;
            const auto geometry = element.geometry();
            const auto vesselLength = geometry.volume();
            const auto vesselVolume = M_PI*vesselRadius*vesselRadius*vesselLength;
            orientations_[eIdx] = geometry.corner(1)-geometry.corner(0);
            orientations_[eIdx] /= orientations_[eIdx].two_norm();
            vesselLength_[eIdx] = vesselLength;

            // compute number of RBCs given the vessel volume and the hematocrit
            // Hct = rbcVolume / vesselVolume * numRBC
            // if there is an endothelial surface layer the maximum hematocrit is smaller than 1.0
            const auto maxTubeHematocrit = rbcRadius*rbcRadius/(vesselRadius*vesselRadius);
            if (initTubeHematocrit > maxTubeHematocrit)
                std::cout << "-- Warning: initial tube hematocrit of "
                          << initTubeHematocrit << " is too high for segment " << eIdx
                          << " and endothelial surface layer thickness of " << eslThickness_ << " m. "
                          << "Using tube Hct = " << maxTubeHematocrit << " instead." << std::endl;

            const auto localTubeHematocrit = std::min(maxTubeHematocrit, initTubeHematocrit);
            const auto numRBC = static_cast<int>(std::floor(localTubeHematocrit*vesselVolume/RedBloodCell::volume()));
            const auto localSpacing = 1.0/static_cast<double>(numRBC);
            const auto localRBCLength = RedBloodCell::volume()/(M_PI*rbcRadius*rbcRadius)/vesselLength;

            if (localRBCLength > 1.0)
                DUNE_THROW(Dune::NotImplemented, "Vessel segments have to at least fit one RBC! At element: "
                                  << eIdx << ", length: " << vesselLength << ", RBC length: " << localRBCLength*vesselLength);

            // assume cylindrical RBCs with same diameter as vessel
            if (numRBC > 0)
            {
                // add the RBCs
                for (int i = 0; i < numRBC; ++i)
                {
                    // get random number between 0 and localSpacing
                    using LocalPosition = typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate;
                    const auto localPos = LocalPosition{dist(gen)*(localSpacing-localRBCLength) + i*localSpacing + 0.5*localRBCLength};
                    rbcs_.emplace_back(rbcIdCounter++, geometry.global(localPos));
                    rbcToElementIdx_.emplace_back(eIdx);
                }
            }

            const auto vIdx0 = gridView_.indexSet().subIndex(element, 0, 1);
            const auto vIdx1 = gridView_.indexSet().subIndex(element, 1, 1);
            nodeNeighbors_[vIdx0].push_back(eIdx);
            nodeNeighbors_[vIdx1].push_back(eIdx);
        }

        // add some inactive buffer RBCs
        const auto numBufferCells = getParam<std::size_t>("RBCModel.NumBufferCells", gridView.size(0)/2);
        rbcs_.reserve(rbcs_.size() + numBufferCells);
        rbcToElementIdx_.resize(rbcToElementIdx_.size() + numBufferCells, 0);

        // a pointer to the inactive RBCs is stored in a queue so we now which RBCs are inactive
        // these RBCs can then be activate by setting the flag and popping them out of the queue
        // deactivated RBCs can be added to the queue as shown here
        for (int i = 0; i < numBufferCells; ++i)
        {
            rbcs_.emplace_back(rbcIdCounter++);
            inactiveRBCs_.push_back(rbcs_.back().deactivate());
        }

        // now all rbcs have persistent locations in the vector
        // build the segment rbc map
        segmentRBCs_.resize(gridView.size(0));
        for (auto& rbc : activeRBCs_())
            segmentRBCs_[rbcToElementIdx_[rbc.id()]].push_back(&rbc);

        std::cout << "Created " << rbcs_.size() << " RBCs ("
                  << rbcs_.size()-numBufferCells << " active, " << numBufferCells << " inactive)" << std::endl;
    }

    //! the particle identifier
    /*!
     * \brief Return the number of rbcs
     * \param onlyActive if true only returns the number of active RBCs (that are inside the domain)
     */
    std::size_t size(bool onlyActive = true) const
    {
        if (onlyActive)
            return rbcs_.size()-inactiveRBCs_.size();
        else
            return rbcs_.size();
    }

    //! get a red blood cell from an id
    const RedBloodCell& rbc(std::size_t id) const
    { return rbcs_[id]; }

    /*!
     * \brief Const range generator for iterating over red blood cells
     * \note usage: for(const auto& rbc : particles(redBloodCellNetwork))
     * This is a free function found by means of ADL
     */
    friend inline Dune::IteratorRange<RBCConstIterator<std::vector<RedBloodCell>>>
    particles(const RBCNetwork& network)
    { return { {network.rbcs_.begin(), network.rbcs_.end()}, {network.rbcs_.end(), network.rbcs_.end()} }; }

    //! get the current orientation of an RBC in the network
    const GlobalPosition& orientation(const RedBloodCell& rbc) const
    { return orientations_[rbcToElementIdx_[rbc.id()]]; }

    //! get the current radius of an RBC in the network
    double radius(const RedBloodCell& rbc) const
    { return radii_[rbcToElementIdx_[rbc.id()]]-eslThickness_; }

    //! get the current length of an RBC in the network
    double rbcLength(const RedBloodCell& rbc) const
    {
        const auto r = radius(rbc);
        return RedBloodCell::volume()/(M_PI*r*r);
    }

    //! get the current element index the RBC belongs to
    std::size_t elementIdx(const RedBloodCell& rbc) const
    { return rbcToElementIdx_[rbc.id()]; }

    //! get the radius of a segment in the network
    double vesselRadius(std::size_t eIdx) const
    { return radii_[eIdx]; }

    //! get the radius of a segment in the network
    double rbcRadius(std::size_t eIdx) const
    { return radii_[eIdx]-eslThickness_; }

    //! get the vessel length of a segment in the network
    double vesselLength(std::size_t eIdx) const
    { return vesselLength_[eIdx]; }

    //! get the orientation of a segment in the network
    const GlobalPosition& orientation(std::size_t eIdx) const
    { return orientations_[eIdx]; }

    //! get the current length of an RBC in segment eIdx in the network
    double rbcLength(std::size_t eIdx) const
    {
        const auto r = rbcRadius(eIdx);
        return RedBloodCell::volume()/(M_PI*r*r);
    }

    //! get pointers to all RBCs of a segment
    const std::deque<RedBloodCell*>& segmentRBCs(std::size_t eIdx) const
    { return segmentRBCs_[eIdx]; }

    //! get pointers to all RBCs of a segment
    std::deque<RedBloodCell*>& segmentRBCs(std::size_t eIdx)
    { return segmentRBCs_[eIdx]; }

    //! get all neighbor elements of a vertex
    const std::vector<std::size_t>& nodeNeighbors(std::size_t vIdx) const
    { return nodeNeighbors_[vIdx]; }

    /*!
     * \brief advance the network
     * \param dt the time step size
     * \param gg the grid geometry
     * \param volumeFluxes the volume flux per element (positive means from 0->1 in local coordinates)
     */
    template<class GridGeometry>
    void advance(double dt, const GridGeometry& gg, const std::vector<double>& volumeFluxes)
    {
        RBCMove move(gg, volumeFluxes, *this);
        move.advance(dt);
    }

    //! add an RBC to the list of inactive RBCs
    void pushInactive(RedBloodCell* rbc)
    { inactiveRBCs_.push_back(rbc); }

    /*!
     * \brief bring one dorming RBC to life! (remove an RBC from the inactive list)
     * \note the RBC is automatically activated
     * \note Take care to place it somewhere in the network afterwards!!! Otherwise it's hanging somewhere in space.
     */
    RedBloodCell* spawnRBC()
    {
        if (inactiveRBCs_.size() == 0)
            DUNE_THROW(Dune::RangeError, "No inactive RBCs left to spawn! Increase RBCModel.NumBufferCells");

        auto* rbc = inactiveRBCs_.front();
        rbc->activate();
        inactiveRBCs_.pop_front();
        return rbc;
    }

    //! if there are any inactive RBCs left in the pool
    bool hasInactiveRBCs() const
    { return !inactiveRBCs_.empty(); }

    //! set the element index of an rbc
    void setElementIndex(RedBloodCell& rbc, std::size_t eIdx)
    { rbcToElementIdx_[rbc.id()] = eIdx; }

    //! get the current inflow spacing for an inflow vertex (local pos and current spacing)
    const std::array<double, 2>& nextInflowRBC(std::size_t vIdx) const
    { return nextInflowRBC_.at(vIdx); }

    //! get the current inflow spacing for an inflow vertex (local pos and current spacing)
    void setNextInflowRBC(std::size_t vIdx, const std::array<double, 2>& rbcInfo)
    { nextInflowRBC_[vIdx] = rbcInfo; }

    //! if the inflow already has been initialized
    bool inflowIsInitialized(std::size_t vIdx) const
    { return nextInflowRBC_.count(vIdx); }

private:
    //! get a range of all active RBCs
    Dune::IteratorRange<RBCIterator<std::vector<RedBloodCell>>> activeRBCs_()
    { return { {rbcs_.begin(), rbcs_.end()}, {rbcs_.end(), rbcs_.end()} }; }

    GridView gridView_;

    std::vector<double> radii_;
    std::vector<GlobalPosition> orientations_;

    // red blood cells
    std::vector<RedBloodCell> rbcs_;
    std::vector<typename GridView::IndexSet::IndexType> rbcToElementIdx_;
    std::deque<RedBloodCell*> inactiveRBCs_;
    std::vector<std::deque<RedBloodCell*>> segmentRBCs_;
    double eslThickness_;

    // additional network information
    std::vector<std::vector<std::size_t>> nodeNeighbors_;
    std::vector<double> vesselLength_;
    std::unordered_map<std::size_t, std::array<double, 2>> nextInflowRBC_;
};

} // end namespace Dumux

#endif

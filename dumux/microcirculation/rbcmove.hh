// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup Microcirculation
 * \brief The algorithms handling the RBC movement
 * \author Timo Koch
 */

#ifndef DUMUX_MICROCIRCULATION_RBC_MOVE_HH
#define DUMUX_MICROCIRCULATION_RBC_MOVE_HH

#include <cmath>
#include <vector>
#include <random>
#include <deque>
#include <unordered_map>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/parameters.hh>

#include <dumux/microcirculation/redbloodcell.hh>
#include <dumux/microcirculation/effectiveresistance.hh>

namespace Dumux {

// forward declaration
template<class GridView> class RBCNetwork;

/**
 * \file
 * \ingroup Microcirculation
 * \brief A class handling the movement of red blood cells (RBCs) in a microcirculation network
 * \author Timo Koch
 * \tparam GridGeometry the grid geometry of the segment grid
 *
 * The `RBCMove` class is responsible for handling the movement of red blood cells (RBCs) in a microcirculation network.
 * It provides methods to calculate the movement of RBCs based on the volume fluxes and the grid geometry.
 *
 * Usage:
 * - Create an instance of `RBCMove` by providing the grid geometry, volume fluxes, and the RBC network.
 * - Call the appropriate methods to calculate and update the RBC movement.
 */
template<class GridGeometry>
class RBCMove
{
    using GridView = typename GridGeometry::GridView;
    using RBCNetwork = Dumux::RBCNetwork<GridView>;
    using GlobalPosition = RedBloodCell::GlobalPosition;

    enum class SegmentType { inner, outflow, none };
    enum class NodeType { outflow, inflow, inner, converging, diverging, none };

public:
    /*!
     * \brief construct a new RBC move
     * \param gg the grid geometry
     * \param volumeFluxes the volume flux per element (positive means from 0->1 in local coordinates)
     * \param the rbcNetwork we are working with
     */
    RBCMove(const GridGeometry& gg, const std::vector<double>& volumeFluxes, RBCNetwork& rbcNetwork)
    : gg_(gg)
    , volumeFluxes_(volumeFluxes)
    , rbcNetwork_(rbcNetwork)
    , marked_(gg.gridView().size(0), false)
    , rndGen_(std::random_device{}())
    {
        // inflow params
        inflowHctIsRandom_ = getParam<bool>("RBCModel.RandomInflowTubeHematocrit", true);
        meanInflowTubeHct_ = getParam<double>("RBCModel.MeanInflowTubeHematocrit", 0.2);
        if (meanInflowTubeHct_ > 1.0)
            DUNE_THROW(ParameterException, "RBCModel.MeanInflowTubeHematocrit is larger than 1.0!");

        // resize the network entity classifiers
        segmentType_.resize(gg.gridView().size(0), SegmentType::none);
        nodeType_.resize(gg.gridView().size(1), NodeType::none);
        // compute local RBC velocities once per move
        localRBCVelocities_.resize(gg.gridView().size(0));

        // set the classifiers (TODO: do I have to do this every time?)
        const auto gridView = gg.gridView();
        for (const auto& element : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(element);
            const auto& flux = volumeFluxes[eIdx];
            bool isBoundary = false;

            const auto& orientation = rbcNetwork_.orientation(eIdx);
            localRBCVelocities_[eIdx] = localRBCVelocity_(eIdx, rbcNetwork_.vesselLength(eIdx), rbcNetwork_.segmentRBCs(eIdx).size());

            auto fvGeometry = localView(gg);
            fvGeometry.bind(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                const auto indexInInside = (scvf.unitOuterNormal()*orientation > 0.0) ? 1 : 0;
                const auto vIdx = gridView.indexSet().subIndex(element, indexInInside, 1);

                if (scvf.boundary())
                {
                    isBoundary = true;
                    if ((indexInInside == 0 && flux < 0.0) || (indexInInside == 1 && flux >= 0.0))
                    {
                        segmentType_[eIdx] = SegmentType::outflow;
                        nodeType_[vIdx] = NodeType::outflow;
                    }
                    else
                    {
                        segmentType_[eIdx] = (segmentType_[eIdx] == SegmentType::none) ? SegmentType::inner : segmentType_[eIdx];
                        nodeType_[vIdx] = NodeType::inflow;
                        // initialize the inflow node
                        if (!rbcNetwork_.inflowIsInitialized(vIdx))
                            initializeInflowNode_(vIdx, computeSegmentInfo_(eIdx, indexInInside, element.geometry()));
                    }
                }

                // inner nodes only have one outside neighbor
                else if (scvf.numOutsideScvs() == 1)
                    nodeType_[vIdx] = NodeType::inner;

                // bifurcation nodes have more
                else if (scvf.numOutsideScvs() > 1)
                {
                    int countInflow = 0; int countOutflow = 0;
                    if ((indexInInside == 0 && flux <= 0.0) || (indexInInside == 1 && flux >= 0.0))
                        countInflow += 1;
                    else
                        countOutflow += 1;

                    for (int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto& flipScvf = fvGeometry.flipScvf(scvf.index(), i);
                        const auto nIdx = flipScvf.insideScvIdx();
                        const auto& nflux = volumeFluxes[nIdx];
                        const auto nIndexInInside = (flipScvf.unitOuterNormal()*rbcNetwork_.orientation(nIdx) > 0.0) ? 1 : 0;
                        if ((nIndexInInside == 0 && nflux <= 0.0) || (nIndexInInside == 1 && nflux >= 0.0))
                            countInflow += 1;
                        else
                            countOutflow += 1;
                    }

                    if (countInflow > 1 && countOutflow == 1)
                        nodeType_[vIdx] = NodeType::converging;
                    else if (countInflow > 0 && countOutflow == 2)
                        nodeType_[vIdx] = NodeType::diverging;
                    else if (countOutflow > 2)
                        DUNE_THROW(Dune::NotImplemented, "RBCMove: Node has more than two outflow segments!");
                    else
                        DUNE_THROW(Dune::InvalidStateException, "RBCMove: Wrong node configuration! "
                                    << "outflow: "  << countOutflow << ", inflow: " << countInflow << ".");
                }
            }

            if (!isBoundary)
                segmentType_[eIdx] = SegmentType::inner;
        }
    }

    //! do a move with time step size dt
    void advance(const double dt)
    {
        // cache the time step size
        timeStepSize_ = dt;

        // (2) Move RBCs starting from outflow segments
        for (std::size_t eIdx = 0; eIdx < gg_.gridView().size(0); ++eIdx)
        {
            if (segmentType_[eIdx] != SegmentType::outflow)
                continue;

            const auto element = gg_.element(eIdx);
            const auto geometry = element.geometry();
            const auto startVIdx = gg_.gridView().indexSet().subIndex(element, 0, 1);
            const auto outflowIndex = nodeType_[startVIdx] == NodeType::outflow ? 0 : 1;

            const auto& orientation = rbcNetwork_.orientation(eIdx);
            const auto vesselLength = rbcNetwork_.vesselLength(eIdx);
            const auto localRBCVelocity = localRBCVelocities_[eIdx];
            const auto& segmentRBCs = rbcNetwork_.segmentRBCs(eIdx);

            // move all RBCs and deactivate those outside the segment
            // here the order doesn't matter
            for (auto* rbc : segmentRBCs)
                rbc->move(localRBCVelocity*dt*vesselLength, orientation);

            // an RBC leaves the domain when more than half of it is outside the domain
            const auto nextRBCLeaves = [&,this]()
            {
                if (segmentRBCs.empty())
                    return false;
                else
                {
                    auto localPos = localPos_(geometry, getFirstRBC_(eIdx, outflowIndex)->position(), /*clamp=*/false);
                    return (outflowIndex == 0) ? (localPos < 0.0) : (localPos > 1.0);
                }
            };

            while (nextRBCLeaves())
                deactivateFirstRBC_(eIdx, outflowIndex);

            // mark this segment as done and move on to next segment
            marked_[eIdx] = true;

            // collect the data the next segement needs from us
            auto segmentInfo = computeSegmentInfo_(eIdx, 1-outflowIndex, geometry);

            // now depending on the neighbor we use different treatments
            const auto otherNodeIdx = gg_.gridView().indexSet().subIndex(element, 1-outflowIndex, 1);
            advance_(otherNodeIdx, std::move(segmentInfo));
        }
    }

private:
    //! continue the movement at the node with index vIdx into the target segment with targetInfo
    template<class SegInfo>
    void advance_(std::size_t vIdx, SegInfo&& targetInfo)
    {
        switch (nodeType_[vIdx])
        {
            // for inner nodes continue moving RBCs
            case NodeType::inner:
                advanceAtInnerNode_(vIdx, std::move(targetInfo));
                break;

            // for converging bifurcations move RBCs one-by-one whoever arrives first
            case NodeType::converging:
                advanceAtConvergingNode_(vIdx, std::move(targetInfo));
                break;

            // for inflow bifurcations spawn elements acoording to inflow hematocrit distribution
            case NodeType::inflow:
                advanceAtInflowNode_(vIdx, std::move(targetInfo));
                break;

            // for diverging bifurcations continue if all branches are ready
            // we have two targets the passed here is only one of them but the info is already computed
            case NodeType::diverging:
                advanceAtDivergingNode_(vIdx, std::move(targetInfo));
                break;

            default:
                DUNE_THROW(Dune::InvalidStateException, "Found unmarked node: " << vIdx);
        }
    }

    //! advance RBCs within a segment (move at most til the end of the segment or the next RBC)
    template<class SegInfo>
    void advanceWithinSegment_(std::size_t vIdx, const SegInfo& segmentInfo)
    {
        // if this element has been already handled do nothing
        if (marked_[segmentInfo.eIdx])
            return;

        const auto& segmentRBCs = rbcNetwork_.segmentRBCs(segmentInfo.eIdx);
        if (segmentRBCs.size() > 0)
        {
            // the local (head) space tells us how far the first RBC can move out of the segment
            const auto localHeadSpace = segmentInfo.localSpace;
            const auto localRBCVelocity = localRBCVelocities_[segmentInfo.eIdx];

            // this algorithm depends on the local direction of movement
            if (segmentInfo.localVIdx == 0)
            {
                assert(localRBCVelocity <= 0.0);
                auto rbcIt = segmentRBCs.begin();
                const auto rbcEndIt = segmentRBCs.end();
                using std::max;
                auto endlocalPos = max(0.0, 0.5*segmentInfo.localRBCLength - localHeadSpace);
                for (; rbcIt != rbcEndIt; ++rbcIt)
                {
                    // find the target position (make sure not to run into the next RBC)
                    const auto localPos = localPos_(segmentInfo.geo, (*rbcIt)->position());
                    auto localTargetPos = localPos + localRBCVelocity*timeStepSize_;
                    localTargetPos = max(localTargetPos, endlocalPos);
                    if(!(localTargetPos > 0.0 - 1e-7 && localTargetPos < 1.0 + 1e-7))
                        DUNE_THROW(Dune::InvalidStateException, "Invalid local position: " << localTargetPos
                                        << "! This vessel should have been blocked!" << " Local RBC length: " << segmentInfo.localRBCLength);

                    // set RBC position
                    (*rbcIt)->setPosition(segmentInfo.geo.global({localTargetPos}));

                    // udpate
                    endlocalPos = localTargetPos + segmentInfo.localRBCLength;
                }
            }
            else
            {
                assert(localRBCVelocity >= 0.0);
                auto rbcIt = segmentRBCs.rbegin();
                const auto rbcEndIt = segmentRBCs.rend();
                using std::min;
                auto endlocalPos = min(1.0, 1.0 - (0.5*segmentInfo.localRBCLength - localHeadSpace));
                for (; rbcIt != rbcEndIt; ++rbcIt)
                {
                    // find the target position (make sure not to run into the next RBC)
                    const auto localPos = localPos_(segmentInfo.geo, (*rbcIt)->position());
                    auto localTargetPos = localPos + localRBCVelocity*timeStepSize_;
                    localTargetPos = min(localTargetPos, endlocalPos);
                    if(!(localTargetPos > 0.0 - 1e-7 && localTargetPos < 1.0 + 1e-7))
                        DUNE_THROW(Dune::InvalidStateException, "Invalid local position: " << localTargetPos
                                        << "! This vessel should have been blocked!" << " Half local RBC length: " << 0.5*segmentInfo.localRBCLength);

                    // set RBC position
                    (*rbcIt)->setPosition(segmentInfo.geo.global({localTargetPos}));

                    // udpate
                    endlocalPos = localTargetPos - segmentInfo.localRBCLength;
                }
            }
        }

        // mark this segment as done and move on to next node
        marked_[segmentInfo.eIdx] = true;

        // collect the data the next segement needs from us (we are the new target)
        auto nextNodeSegmentInfo = computeSegmentInfo_(segmentInfo.eIdx, 1-segmentInfo.localVIdx, segmentInfo.geo);

        // now depending on the node type we use different treatments
        const auto element = gg_.element(segmentInfo.eIdx);
        const auto otherNodeIdx = gg_.gridView().indexSet().subIndex(element, 1-segmentInfo.localVIdx, 1);
        advance_(otherNodeIdx, std::move(nextNodeSegmentInfo));
    }

    //! advance RBCs at an inner node
    template<class SegInfo>
    void advanceAtInnerNode_(std::size_t vIdx, SegInfo&& targetInfo)
    {
        // there is only one source element
        assert(rbcNetwork_.nodeNeighbors(vIdx).size() == 2);
        const auto eIdx = (rbcNetwork_.nodeNeighbors(vIdx)[0] == targetInfo.eIdx) ? rbcNetwork_.nodeNeighbors(vIdx)[1] : rbcNetwork_.nodeNeighbors(vIdx)[0];
        const auto element = gg_.element(eIdx);
        const auto geo = element.geometry();
        const auto localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;

        const auto localRBCVelocity = localRBCVelocities_[eIdx];
        const auto& segmentRBCs = rbcNetwork_.segmentRBCs(eIdx);
        // whether there is traffic jam in the target
        bool targetBlocked = targetInfo.localSpace < 0.501*targetInfo.localRBCLength;
        if (!targetBlocked && !segmentRBCs.empty())
        {
            using std::abs; using std::max;
            const auto localMoveDist = timeStepSize_*abs(localRBCVelocity);
            double localDist = max(0.0, abs(localPos_(geo, getFirstRBC_(eIdx, localVIdx)->position()) - 1.0*localVIdx));

            // as long as the move is larger than the local space and the target has space the first RBC moves to the target segment
            while (localMoveDist > localDist && !targetBlocked && !segmentRBCs.empty())
            {
                using std::min; using std::max;
                const auto dtTarget = timeStepSize_ - localDist/abs(localRBCVelocity);
                auto localTargetPos = min(max(dtTarget*localRBCVelocities_[targetInfo.eIdx], 0.0), targetInfo.localSpace);
                if (targetInfo.localVIdx == 1)
                    localTargetPos = 1.0 - localTargetPos;

                // move the position of the rbc
                auto* rbc = getFirstRBC_(eIdx, localVIdx);
                rbc->setPosition(targetInfo.geo.global({localTargetPos}));

                // move the rbc in the network
                popFirstRBC_(eIdx, localVIdx);
                pushFirstRBC_(rbc, targetInfo.eIdx, targetInfo.localVIdx);

                // update the local spaces
                updateTargetLocalSpace_(targetInfo);
                targetBlocked = targetInfo.localSpace < 0.501*targetInfo.localRBCLength;
                using std::max; using std::abs;
                if (!segmentRBCs.empty())
                    localDist = max(0.0, abs(localPos_(geo, getFirstRBC_(eIdx, localVIdx)->position()) - 1.0*localVIdx));
            }
        }

        // the rest of the cells stay within the segment and can be moved internally
        // however, there might be traffic jam (localMoveDist > localDist && targetBlocked)
        // so we need to take care to move at most the maximum allowed distance til the end of the segment or the next RBC
        // the local (head) space tells us how far the first RBC can move out of the segment
        auto segmentInfo = computeSegmentInfo_(eIdx, localVIdx, geo);
        segmentInfo.localSpace = targetInfo.localSpace*rbcNetwork_.vesselLength(targetInfo.eIdx)/rbcNetwork_.vesselLength(eIdx);
        advanceWithinSegment_(vIdx, segmentInfo);
    }

    //! advance RBCs at an inner node
    template<class SegInfo>
    void advanceAtConvergingNode_(std::size_t vIdx, SegInfo&& targetInfo)
    {
        // there are multiple sources, assume a maxium of 4
        Dune::ReservedVector<std::size_t, 4> sourceElems;
        for (const auto nIdx : rbcNetwork_.nodeNeighbors(vIdx))
            if (nIdx != targetInfo.eIdx)
                sourceElems.push_back(nIdx);

        // whether there is traffic jam in the target
        bool targetBlocked = targetInfo.localSpace < 0.501*targetInfo.localRBCLength;
        if (!targetBlocked)
        {
            // find the closest rbc (time-wise)
            const auto vertexPos = gg_.element(targetInfo.eIdx).template subEntity<1>(targetInfo.localVIdx).geometry().corner(0);
            auto [rbc, eIdx, localDist] = closestRBC_(vertexPos, sourceElems);
            if (rbc)
            {
                using std::abs;
                auto localMoveDist = timeStepSize_*abs(localRBCVelocities_[eIdx]);
                auto localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;

                // as long as the move is larger than the local space and the target has space the first RBC moves to the target segment (as long as upstream RBCs are found)
                while (localMoveDist > localDist && !targetBlocked && rbc)
                {
                    using std::min; using std::max;
                    const auto dtTarget = timeStepSize_ - localDist/abs(localRBCVelocities_[eIdx]);
                    auto localTargetPos = min(max(dtTarget*localRBCVelocities_[targetInfo.eIdx], 0.0), targetInfo.localSpace);
                    if (targetInfo.localVIdx == 1)
                        localTargetPos = 1.0 - localTargetPos;

                    // move the position of the rbc
                    rbc->setPosition(targetInfo.geo.global({localTargetPos}));

                    // move the rbc in the network
                    popFirstRBC_(eIdx, localVIdx);
                    pushFirstRBC_(rbc, targetInfo.eIdx, targetInfo.localVIdx);

                    // update the local spaces
                    updateTargetLocalSpace_(targetInfo);
                    targetBlocked = targetInfo.localSpace < 0.501*targetInfo.localRBCLength;

                    // update the closest RBC and the source data
                    std::tie(rbc, eIdx, localDist) = closestRBC_(vertexPos, sourceElems);
                    if (rbc)
                    {
                        localMoveDist = timeStepSize_*abs(localRBCVelocities_[eIdx]);
                        localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;
                    }
                }
            }
        }

        // traffic jam or done, advance all upstream branches locally
        for (const auto eIdx : sourceElems)
        {
            const auto localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;
            auto segmentInfo = computeSegmentInfo_(eIdx, localVIdx, gg_.element(eIdx).geometry());
            segmentInfo.localSpace = targetInfo.localSpace*rbcNetwork_.vesselLength(targetInfo.eIdx)/rbcNetwork_.vesselLength(eIdx);
            advanceWithinSegment_(vIdx, segmentInfo);
        }
    }

    //! initialize local information at inflow node
    template<class SegInfo>
    void initializeInflowNode_(std::size_t vIdx, SegInfo&& targetInfo)
    {
        // compute the local position of the next rbc
        const auto vesselLength = rbcNetwork_.vesselLength(targetInfo.eIdx);
        const auto currentLocalSpacing = generateSpacing_(rbcNetwork_.rbcLength(targetInfo.eIdx), rbcNetwork_.rbcRadius(targetInfo.eIdx), rbcNetwork_.vesselRadius(targetInfo.eIdx))/vesselLength;
        auto nextLocalPos = targetInfo.localSpace - currentLocalSpacing - 0.501*targetInfo.localRBCLength;
        if (targetInfo.localVIdx == 1)
            nextLocalPos = 1.0 - nextLocalPos;
        rbcNetwork_.setNextInflowRBC(vIdx, {nextLocalPos, currentLocalSpacing});
    }

    //! advance RBCs at an inflow node
    template<class SegInfo>
    void advanceAtInflowNode_(std::size_t vIdx, SegInfo&& targetInfo)
    {
        // move the virtual inflow rbc
        auto localTargetPos = rbcNetwork_.nextInflowRBC(vIdx)[0] + timeStepSize_*localRBCVelocities_[targetInfo.eIdx];
        auto currentLocalSpacing = rbcNetwork_.nextInflowRBC(vIdx)[1];

        // if there is another RBC in the inflow segment make sure to keep distance
        using std::min; using std::max;
        if (targetInfo.localSpace < 0.99)
        {
            if (targetInfo.localVIdx == 1)
                localTargetPos = max(localTargetPos, 1.0 - (targetInfo.localSpace - currentLocalSpacing - 0.501*targetInfo.localRBCLength));
            else
                localTargetPos = min(localTargetPos, (targetInfo.localSpace - currentLocalSpacing - 0.501*targetInfo.localRBCLength));
        }

        // add new RBCs as long as their local position is inside the domain
        while (localTargetPos > 0.0 - 1e-5 && localTargetPos < 1.0 + 1e-5)
        {
            // spawn new RBC with this position
            auto* rbc = rbcNetwork_.spawnRBC();
            rbc->setPosition(targetInfo.geo.global({localTargetPos}));
            pushFirstRBC_(rbc, targetInfo.eIdx, targetInfo.localVIdx); // add to network

            // update local spaces
            using std::max; using std::abs;
            updateTargetLocalSpace_(targetInfo);
            currentLocalSpacing = generateSpacing_(rbcNetwork_.rbcLength(targetInfo.eIdx), rbcNetwork_.rbcRadius(targetInfo.eIdx), rbcNetwork_.vesselRadius(targetInfo.eIdx))/rbcNetwork_.vesselLength(targetInfo.eIdx);
            localTargetPos = targetInfo.localSpace - currentLocalSpacing - 0.501*targetInfo.localRBCLength;
            if (targetInfo.localVIdx == 1)
                localTargetPos = 1.0 - localTargetPos;
        }

        // update persistent data for next time step
        rbcNetwork_.setNextInflowRBC(vIdx, {localTargetPos, currentLocalSpacing});
    }

    //! advance RBCs at a diverging node (bifurcation rule)
    template<class SegInfo>
    void advanceAtDivergingNode_(std::size_t vIdx, SegInfo&& firstTargetInfo)
    {
        // find other target and maybe multiple source elements
        Dune::ReservedVector<std::size_t, 4> sourceElems;
        std::size_t secondTargetEIdx{}, secondTargetLocalVIdx{};
        for (const auto nIdx : rbcNetwork_.nodeNeighbors(vIdx))
        {
            if (nIdx == firstTargetInfo.eIdx)
                continue;

            const auto element = gg_.element(nIdx);
            const auto testVIdx = gg_.gridView().indexSet().subIndex(element, 0, 1);
            const auto localVIdx = testVIdx == vIdx ? 0 : 1;
            if ((localVIdx == 0 && volumeFluxes_[nIdx] < 0.0) || (localVIdx == 1 && volumeFluxes_[nIdx] >= 0.0))
                sourceElems.push_back(nIdx);
            else // this is the other target, TODO this is ugly
            {
                secondTargetEIdx = nIdx;
                secondTargetLocalVIdx = localVIdx;
            }
        }

        // we only deal with this node if both targets have been already handled
        if (!marked_[secondTargetEIdx] || !marked_[firstTargetInfo.eIdx])
            return;

        std::array<SegInfo, 2> targetInfos{{std::move(firstTargetInfo),
                                            computeSegmentInfo_(secondTargetEIdx, secondTargetLocalVIdx,
                                                                gg_.element(secondTargetEIdx).geometry())}};

        bool bothTargetsBlocked = targetInfos[0].localSpace < 0.501*targetInfos[0].localRBCLength && targetInfos[1].localSpace < 0.501*targetInfos[1].localRBCLength;
        if (!bothTargetsBlocked)
        {
            using std::abs;
            const auto vertexPos = gg_.element(targetInfos[0].eIdx).template subEntity<1>(targetInfos[0].localVIdx).geometry().corner(0);
            auto [rbc, eIdx, localDist] = closestRBC_(vertexPos, sourceElems);
            auto localMoveDist = timeStepSize_*abs(localRBCVelocities_[eIdx]);

            // do we have a bifurcation event?
            if (rbc && localMoveDist > localDist)
            {
                // compute local index of this source element
                auto localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;

                // as long as the move is larger than the local space and the target has space the first RBC moves to one of the target segments (as long as upstream RBCs are found)
                while (localMoveDist > localDist && !bothTargetsBlocked && rbc)
                {
                    using std::min; using std::max;
                    const auto dtTarget = timeStepSize_ - localDist/abs(localRBCVelocities_[eIdx]);
                    auto& curTargetInfo = targetInfos[pickBifurcationTarget_(targetInfos)];
                    if (curTargetInfo.localSpace < 0.501*curTargetInfo.localRBCLength)
                        DUNE_THROW(Dune::InvalidStateException, "A blocked target has been chosen! Check the bifurcation rule!");

                    auto localTargetPos = min(max(dtTarget*localRBCVelocities_[curTargetInfo.eIdx], 0.0), curTargetInfo.localSpace);
                    if (curTargetInfo.localVIdx == 1)
                        localTargetPos = 1.0 - localTargetPos;

                    // move the position of the rbc
                    rbc->setPosition(curTargetInfo.geo.global({localTargetPos}));

                    // move the rbc in the network
                    popFirstRBC_(eIdx, localVIdx);
                    pushFirstRBC_(rbc, curTargetInfo.eIdx, curTargetInfo.localVIdx);

                    // update the local spaces
                    updateTargetLocalSpace_(curTargetInfo);
                    bothTargetsBlocked = targetInfos[0].localSpace < 0.501*targetInfos[0].localRBCLength && targetInfos[1].localSpace < 0.501*targetInfos[1].localRBCLength;

                    // update the closest RBC and the source data
                    std::tie(rbc, eIdx, localDist) = closestRBC_(vertexPos, sourceElems);
                    if (rbc)
                    {
                        localMoveDist = timeStepSize_*abs(localRBCVelocities_[eIdx]);
                        localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;
                    }
                }
            }
        }

        // traffic jam or done, advance all upstream branches locally
        for (const auto eIdx : sourceElems)
        {
            const auto localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;
            auto segmentInfo = computeSegmentInfo_(eIdx, localVIdx, gg_.element(eIdx).geometry());
            const auto targetIt = std::max_element(targetInfos.begin(), targetInfos.end(), [](const auto& a, const auto& b){ return a.localSpace < b.localSpace; });
            segmentInfo.localSpace = targetIt->localSpace*rbcNetwork_.vesselLength(targetIt->eIdx)/rbcNetwork_.vesselLength(eIdx);
            advanceWithinSegment_(vIdx, segmentInfo);
        }
    }

    //! compute the local RBC velocity
    double localRBCVelocity_(const std::size_t eIdx, const double vesselLength, const std::size_t numRBC)
    {
        using std::max;
        const auto r = rbcNetwork_.vesselRadius(eIdx);
        const auto volume = M_PI*r*r*vesselLength;
        const auto bulkVelocity = volumeFluxes_[eIdx]/(M_PI*r*r);
        // regularize on the lower end
        const auto tubeHematocrit = RedBloodCell::volume()*numRBC/volume;
        const auto dischargeHematocrit = Dumux::dischargeHematocritInVitro(r, tubeHematocrit);
        const auto rbcVelocity = Dumux::rbcVelocity(bulkVelocity, dischargeHematocrit, tubeHematocrit);
        return rbcVelocity/vesselLength;
    }

    //! get the current front runner
    RedBloodCell* getFirstRBC_(std::size_t eIdx, std::size_t localVIdx)
    {
        if (rbcNetwork_.segmentRBCs(eIdx).empty())
            DUNE_THROW(Dune::InvalidStateException, "Cannot get first RBC from a segment without RBCs!");
        return localVIdx == 1 ? rbcNetwork_.segmentRBCs(eIdx).back() : rbcNetwork_.segmentRBCs(eIdx).front();
    }

    //! get the current front runner
    const RedBloodCell* getFirstRBC_(std::size_t eIdx, std::size_t localVIdx) const
    {
        if (rbcNetwork_.segmentRBCs(eIdx).empty())
            DUNE_THROW(Dune::InvalidStateException, "Cannot get first RBC from a segment without RBCs!");
        return localVIdx == 1 ? rbcNetwork_.segmentRBCs(eIdx).back() : rbcNetwork_.segmentRBCs(eIdx).front();
    }

    //! deactivate the current front runner (because it left the domain)
    void deactivateFirstRBC_(std::size_t eIdx, std::size_t localVIdx)
    {
        auto& segmentRBCs = rbcNetwork_.segmentRBCs(eIdx);
        if (localVIdx == 1)
        {
            rbcNetwork_.pushInactive(segmentRBCs.back()->deactivate());
            segmentRBCs.pop_back();
        }
        else
        {
            rbcNetwork_.pushInactive(segmentRBCs.front()->deactivate());
            segmentRBCs.pop_front();
        }
    }

    //! pop the current front runner (because it jumped to another segment)
    void popFirstRBC_(std::size_t eIdx, std::size_t localVIdx)
    {
        (localVIdx == 1) ? rbcNetwork_.segmentRBCs(eIdx).pop_back() : rbcNetwork_.segmentRBCs(eIdx).pop_front();
    }

    //! push a new element at the beginning
    void pushFirstRBC_(RedBloodCell* rbc, std::size_t eIdx, std::size_t localVIdx)
    {
        rbcNetwork_.setElementIndex(*rbc, eIdx); // reset the index to the new element
        (localVIdx == 1) ? rbcNetwork_.segmentRBCs(eIdx).push_back(rbc) : rbcNetwork_.segmentRBCs(eIdx).push_front(rbc);
    }

    //! segment information
    template<class Geometry>
    struct SegmentInfo { const std::size_t eIdx, localVIdx; const double localRBCLength; const Geometry geo; double localSpace; };

    //! compute local information about a segment
    template<class Geometry>
    SegmentInfo<Geometry> computeSegmentInfo_(std::size_t eIdx, std::size_t localVIdx, const Geometry& geo) const
    {
        using std::abs; using std::max;
        const double localRBCLength = rbcNetwork_.rbcLength(eIdx)/rbcNetwork_.vesselLength(eIdx);
        const double localSpace = localTargetSpace_(eIdx, localVIdx, geo, localRBCLength);
        return { eIdx, localVIdx, localRBCLength, geo, localSpace };
    }

    //! compute local space in a target
    template<class SegInfo>
    void updateTargetLocalSpace_(SegInfo& targetInfo) const
    { targetInfo.localSpace = localTargetSpace_(targetInfo.eIdx, targetInfo.localVIdx, targetInfo.geo, targetInfo.localRBCLength); }

    //! compute local space in a target
    template<class Geometry>
    double localTargetSpace_(std::size_t eIdx, std::size_t localVIdx, const Geometry& geo, const double localRBCLength) const
    {
        const auto& segmentRBCs = rbcNetwork_.segmentRBCs(eIdx);
        if (segmentRBCs.empty())
            return 1.0;

        using std::abs; using std::max; using std::min;
        const auto spaceToNextRBC = max(0.0, abs(localPos_(geo, getFirstRBC_(eIdx, localVIdx)->position()) - 1.0*localVIdx) - 0.501*localRBCLength);
        return max(0.0, min(1.0-localRBCLength*segmentRBCs.size(), spaceToNextRBC));
    }

    //! get the local coordinate from a global coordinate (orthogonal projection)
    template<class Geometry>
    double localPos_(const Geometry& geometry, const GlobalPosition& point, bool clamp = true) const
    {
        const auto& p = geometry.corner(0);
        const auto& q = geometry.corner(1);
        const auto pq = q-p;
        using std::min; using std::max;
        return clamp ? max(0.0, min(1.0, (pq)*(point-p)/(pq*pq))) : (pq)*(point-p)/(pq*pq);
    }

    //! find the closest RBC coming towards a node (when multiple source branches)
    template<class SourceElems>
    std::tuple<RedBloodCell*, std::size_t, double> closestRBC_(const GlobalPosition& vPos, const SourceElems& sourceElems)
    {
        double minDt = 1e100;
        auto closestRBC = std::tuple<RedBloodCell*, std::size_t, double>{nullptr, 0, 0.0};
        for (const auto eIdx : sourceElems)
        {
            const auto localVIdx = volumeFluxes_[eIdx] > 0.0 ? 1 : 0;
            const auto& segmentRBCs = rbcNetwork_.segmentRBCs(eIdx);
            if (segmentRBCs.empty())
                continue;
            else
            {
                using std::max; using std::abs;
                auto* rbc = getFirstRBC_(eIdx, localVIdx);
                const auto localDist = (rbc->position()-vPos).two_norm()/rbcNetwork_.vesselLength(eIdx);
                const auto dt = localDist/abs(localRBCVelocities_[eIdx]);
                if (dt < minDt)
                {
                    minDt = dt;
                    closestRBC = std::tuple<RedBloodCell*, std::size_t, double>{rbc, eIdx, localDist};
                }
            }
        }
        return closestRBC;
    }

    //! generate a random spacing for the boundary condition at inflow nodes (see Schmid (2017) Section 11.1.4)
    double generateSpacing_(const double rbcLength, const double rbcRadius, const double vesselRadius)
    {
        using std::min; using std::log;
        const auto maxTubeHematocrit = rbcRadius*rbcRadius/(vesselRadius*vesselRadius);
        const auto localTubeHematocrit = min(maxTubeHematocrit, meanInflowTubeHct_);
        if (meanInflowTubeHct_ > maxTubeHematocrit)
            std::cout << "-- Warning: mean inflow tube hematocrit of "
                      << meanInflowTubeHct_ << " is too high."
                      << " Using tube Hct = " << maxTubeHematocrit << " instead." << std::endl;

        if (inflowHctIsRandom_)
        {
            // TODO: check this (discharge hematocrit? standard deviation too low?)
            std::lognormal_distribution<double> dist(log(rbcLength*(maxTubeHematocrit - localTubeHematocrit)/localTubeHematocrit), 0.1);
            return dist(rndGen_);
        }
        else
            return rbcLength*(maxTubeHematocrit - localTubeHematocrit)/localTubeHematocrit;
    }

    //! pick one target branch at bifurcation
    //!  This assumes that at least on of the targets is not blocked!
    template<class SegInfo>
    int pickBifurcationTarget_(const std::array<SegInfo, 2>& infos)
    {
        // is the first branch is blocked take the second branch (we know know that not both are blocked)
        if (infos[0].localSpace < 0.501*infos[0].localRBCLength)
            return 1;
        else if (infos[1].localSpace < 0.501*infos[1].localRBCLength)
            return 0;

        // otherwise return the one with the highest bulk velocity
        // TODO: Bifucation rule for larger vessels
        else
        {
            const std::array<double, 2> r = {{ rbcNetwork_.vesselRadius(infos[0].eIdx), rbcNetwork_.vesselRadius(infos[1].eIdx) }};
            using std::abs;
            if (abs(volumeFluxes_[infos[0].eIdx]/(M_PI*r[0]*r[0])) > abs(volumeFluxes_[infos[1].eIdx]/(M_PI*r[1]*r[1])))
                return 0;
            else
                return 1;
        }
    }

    const GridGeometry& gg_;
    const std::vector<double>& volumeFluxes_;
    RBCNetwork& rbcNetwork_;
    std::vector<bool> marked_; //!< flags which segment has already been handled
    std::vector<double> localRBCVelocities_;

    // inflow
    double meanInflowTubeHct_;
    bool inflowHctIsRandom_;

    double timeStepSize_ = 0.0;

    // inflow boundaries
    std::mt19937 rndGen_;

    // additional network information
    std::vector<SegmentType> segmentType_;
    std::vector<NodeType> nodeType_;
};

} // end namespace Dumux

#endif

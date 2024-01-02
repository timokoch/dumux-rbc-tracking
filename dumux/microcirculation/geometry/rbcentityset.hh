// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \author Timo Koch
 * \ingroup Microcirculation
 * \brief An interface for a set of rbc geometries
 * \note This can be used e.g. to contruct a bounding box volume hierarchy of an rbc network
 */
#ifndef DUMUX_MICROCIRCULATION_RBC_ENTITY_SET_HH
#define DUMUX_MICROCIRCULATION_RBC_ENTITY_SET_HH

#include <memory>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/microcirculation/redbloodcell.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief An interface for a set of geometric entities
 * \note This can be used e.g. to contruct a bounding box volume hierarchy of a grid
 * It defines the minimum requirement for such a set
 */
template <class RBCNetwork>
class RBCEntitySet
{
    /*!
     * \brief The current geometrical configuration of an RBC
     */
    class RBCGeometry
    {
        using GlobalPosition = RedBloodCell::GlobalPosition;
        using LocalPosition = Dune::FieldVector<GlobalPosition::value_type, 1>;
    public:

        static constexpr int mydimension = 1;
        static constexpr int coorddimension = 3;

        using ctype = GlobalPosition::value_type;
        using GlobalCoordinate = GlobalPosition;

        RBCGeometry(const RedBloodCell& rbc, const RBCNetwork& rbcNetwork)
        {
            const auto& orientation = rbcNetwork.orientation(rbc);
            const auto length = rbcNetwork.rbcLength(rbc);

            corners_[0] = rbc.position(); corners_[0].axpy(-0.5*length, orientation);
            corners_[1] = rbc.position(); corners_[1].axpy(0.5*length, orientation);
        }

        //! the global position given a local position
        GlobalPosition global(const LocalPosition& localPos) const
        {
            auto globalPos = corners_[0];
            globalPos.axpy(localPos[0], corners_[1]-corners_[0]);
            return globalPos;
        }

        //! the number of corners
        static constexpr unsigned int corners()
        { return 2; }

        //! get corner by local index
        const GlobalPosition& corner(unsigned int i) const
        { return corners_[i]; }

    private:
        std::array<GlobalPosition, 2> corners_;
    };

    /*!
     * \brief Wrapper class around an RBC turning it into a geometric entity
     */
    class RBCEntity
    {
    public:
        using Geometry = RBCGeometry;

        //! wrap rbc into a geometric entity
        RBCEntity(const RedBloodCell& rbc, const RBCNetwork& rbcNetwork, std::size_t index, std::size_t id)
        : geo_(rbc, rbcNetwork)
        , index_(index)
        , id_(id)
        {}

        const Geometry& geometry() const
        { return geo_; }

        std::size_t index() const
        { return index_; }

        std::size_t id() const
        { return id_; }

    private:
        Geometry geo_;
        std::size_t index_, id_;
    };

public:
    using Entity = RBCEntity;

    RBCEntitySet(std::shared_ptr<const RBCNetwork> rbcNetwork)
    : rbcNetwork_(rbcNetwork)
    {
        std::size_t counter = 0;
        entities_.reserve(rbcNetwork->size());
        for (const auto& rbc : particles(*rbcNetwork))
            entities_.emplace_back(rbc, *rbcNetwork, counter++, rbc.id());
    }

    /*!
     * \brief The world dimension of the entity set
     */
    static constexpr int dimensionworld = 3;

    /*!
     * \brief the coordinate type
     */
    using ctype = double;

    /*!
     * \brief the number of entities in this set
     */
    decltype(auto) size() const
    { return rbcNetwork_->size(); }

    /*!
     * \brief begin iterator to enable range-based for iteration
     */
    decltype(auto) begin() const
    { return entities_.begin(); }

    /*!
     * \brief end iterator to enable range-based for iteration
     */
    decltype(auto) end() const
    { return entities_.end(); }

    /*!
     * \brief get an entities index
     */
    std::size_t index(const Entity& e) const
    { return e.index(); }

    /*!
     * \brief get an entity from an index
     */
    Entity entity(std::size_t index) const
    { return entities_[index]; }

    /*!
     * \brief get an rbc from an index
     */
    const RedBloodCell& rbc(std::size_t index) const
    { return rbcNetwork_->rbc(entities_[index].id()); }

private:
    std::shared_ptr<const RBCNetwork> rbcNetwork_;
    std::vector<RBCEntity> entities_;
};

} // end namespace Dumux

#endif

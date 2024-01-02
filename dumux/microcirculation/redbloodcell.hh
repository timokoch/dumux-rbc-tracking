// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup Microcirculation
 * \brief A red blood cell particle
 * \author Timo Koch
 */

#ifndef DUMUX_MICROCIRCULATION_REDBLOODCELL_HH
#define DUMUX_MICROCIRCULATION_REDBLOODCELL_HH

#include <dune/common/fvector.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/**
 * \ingroup Microcirculation
 * \brief A red blood cell particle
 * \author Timo Koch
 */
class RedBloodCell
{
public:
    // we always consider three-dimensional RBCs and tube networks
    using GlobalPosition = Dune::FieldVector<double, 3>;

    /**
     * \brief create a new red blood cell (RBC)
     * \param id a unique identifier
     * \param pos the global position of the RBC
     * \param active if the RBC is currently in active state
     */
    RedBloodCell(std::size_t id, const GlobalPosition& pos, bool active = true)
    : id_(id)
    , position_(pos)
    , active_(active)
    {}

    /**
     * \brief Construct inactive RBC
     */
    RedBloodCell(std::size_t id)
    : id_(id)
    , position_(0)
    , active_(false)
    {}

    //! the current global position of the RBC
    const GlobalPosition& position() const
    { return position_; }

    //! set a new global position
    void setPosition(const GlobalPosition& pos)
    { position_ = pos; }

    //! move by length in direction
    void move(double step, const GlobalPosition& dir)
    { position_.axpy(step, dir); }

    /**
     * \brief the RBC volume in cubic meter
     * \note Pries et al (1990): rat 55fl, human 92fl
     */
    static double volume()
    {
        static const double volume = getParam<double>("RBC.Volume", 92e-18);
        return volume;
    }

    //! the particle identifier
    std::size_t id() const
    { return id_; }

    /**
     * \brief deactivate this rbc (e.g. leaves model domain)
     * \return a pointer to this RBC
     */
    RedBloodCell* deactivate()
    { active_ = false; return this; }

    /**
     * \brief activate this rbc (e.g. enters model domain)
     */
    void activate()
    { active_ = true; }

    /**
     * \brief is the RBC currently in active mode?
     */
    bool isActive() const
    { return active_; }

private:
    std::size_t id_;
    GlobalPosition position_;
    bool active_;
};

} // end namespace Dumux

#endif

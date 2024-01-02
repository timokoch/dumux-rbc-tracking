// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup Microcirculation
 * \brief Iterators over red blood cells
 * \author Timo Koch
 */

#ifndef DUMUX_MICROCIRCULATION_RBC_ITERATOR_HH
#define DUMUX_MICROCIRCULATION_RBC_ITERATOR_HH

#include <vector>
#include <iterator>

#include <dune/common/iteratorrange.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/iteratorrange.hh>

#include <dumux/microcirculation/redbloodcell.hh>

namespace Dumux {

/*!
 * \ingroup Microcirculation
 * \brief Iterators over red blood cells
 * \note usage: for(auto&& rbc : particles(redBloodCellNetwork))
 */
template<class Container>
class RBCIterator
: public Dune::ForwardIteratorFacade<RBCIterator<Container>, RedBloodCell>
{
    using Iterator = typename Container::iterator;
public:
    //! default construtor
    RBCIterator() = default;

    //! create from iterator
    RBCIterator(const Iterator& it, const Iterator& endIt)
    : it_{it}, endIt_{endIt}
    {
        while (it_ != endIt_ && !it_->isActive())
            ++it_;
    }

    //! dereferencing yields a red blood cell
    RedBloodCell& dereference() const
    {
        return *it_;
    }

    //! test for equality
    bool equals(const RBCIterator& other) const
    {
        return it_ == other.it_;
    }

    //! increment until the next active RBC
    void increment()
    {
        ++it_;
        while (it_ != endIt_ && !it_->isActive())
            ++it_;
    }

private:
    Iterator it_;
    const Iterator endIt_;
};

/*!
 * \ingroup Microcirculation
 * \brief Iterators over red blood cells
 * \note usage: for(auto&& rbc : particles(redBloodCellNetwork))
 */
template<class Container>
class RBCConstIterator
: public Dune::ForwardIteratorFacade<RBCConstIterator<Container>, const RedBloodCell>
{
    using Iterator = typename Container::const_iterator;
public:
    //! default construtor
    RBCConstIterator() = default;

    //! create from iterator
    RBCConstIterator(const Iterator& it, const Iterator& endIt)
    : it_{it}, endIt_{endIt}
    {
        while (it_ != endIt_ && !it_->isActive())
            ++it_;
    }

    //! dereferencing yields a red blood cell
    const RedBloodCell& dereference() const
    {
        return *it_;
    }

    //! test for equality
    bool equals(const RBCConstIterator& other) const
    {
        return it_ == other.it_;
    }

    //! increment until the next active RBC
    void increment()
    {
        ++it_;
        while (it_ != endIt_ && !it_->isActive())
            ++it_;
    }

private:
    Iterator it_;
    const Iterator endIt_;
};

} // end namespace Dumux

#endif

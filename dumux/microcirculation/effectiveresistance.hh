// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup Microcirculation
 * \brief Functional relationships for effective vessel resistance computation
 * \author Timo Koch
 */

#ifndef DUMUX_MICROCIRCULATION_EFFECTIVE_RESISTANCE_HH
#define DUMUX_MICROCIRCULATION_EFFECTIVE_RESISTANCE_HH

#include <cmath>

namespace Dumux {

/**
 * \file
 * \ingroup Microcirculation
 * \brief Compute the discharge hematocrit according to Pries et al (1990)
 */
inline double dischargeHematocritInVitro(const double radius, const double tubeHematocrit)
{
    // this is the quadratic formula from Pries et al (1992) Eq. (10)
    // with the constants from Pries et al (1990) Eq. (1)
    // also see Schmid (2017) Eq. (8.9)
    using std::exp; using std::sqrt;
    const auto d = 2.0*radius*1e6;
    const auto x = 1.0 + 1.7*exp(-0.415*d) - 0.6*exp(-0.011*d);
    const auto xFactor = -x/(2.0 - 2.0*x);
    return xFactor + sqrt(xFactor*xFactor + tubeHematocrit/(1.0-x));
}

/**
 * \file
 * \ingroup Microcirculation
 * \brief Compute rbc velocity, see Schmid (2017) Eq. (8.10)
 */
inline double rbcVelocity(const double bulkVelocity, const double dischargeHematocrit, const double tubeHematocrit)
{
    if (tubeHematocrit < 1e-5)
        return bulkVelocity;
    else
        return bulkVelocity*dischargeHematocrit/tubeHematocrit;
}

/**
 * \file
 * \ingroup Microcirculation
 * \brief Compute the relative effective viscosity (in vitro), see Schmid (2017) Eq. (8.13)
 */
inline double relEffViscosityInVitro(const double radius, const double dischargeHematocrit)
{
    if (dischargeHematocrit < 1e-5)
        return 1.0;
    else
    {
        using std::exp; using std::sqrt; using std::pow;
        const auto d = 2.0*radius*1e6;
        const auto dFactor = 1.0/(1.0 + 1e-11*pow(d, 12));
        const auto c = (0.8 + exp(-0.075*d))*(dFactor - 1.0) + dFactor;
        const auto mu45 = 220.0*exp(-1.3*d) + 3.2 - 2.44*exp(-0.06*pow(d, 0.645));
        return 1.0 + (mu45 - 1.0)*(pow(1.0 - dischargeHematocrit, c) - 1.0)/(pow(1.0 - 0.45, c) - 1.0);
    }
}

} // end namespace Dumux

#endif

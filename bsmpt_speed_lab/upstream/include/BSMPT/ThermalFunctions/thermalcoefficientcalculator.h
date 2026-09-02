// Copyright (C) 2021  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <functional>
#include <vector>

namespace BSMPT
{
namespace ThermalFunctions
{
/**
 * @brief The ThermalCoefficientCalculator class is a thread-safe wrapper around
 * the calculation of the coefficients for the thermal expansions used in
 * ThermalFunctions
 */
class ThermalCoefficientCalculator
{
public:
  /**
   * @brief ThermalCoefficientCalculator
   * @param func Describes the function to calculate the coefficents
   * @param maxOrderToPrecalc defines until which order the coefficents are
   * precalculated
   */
  ThermalCoefficientCalculator(std::function<double(int)> func,
                               int maxOrderToPrecalc);
  /**
   * @brief GetCoefficentAtOrder
   * @param n describes the order at which the coefficient should be calculated
   * @return The coefficent at the given order. If it was precalculated the
   * result stored in the map will be returned, otherwise the result will be
   * calculated.
   */
  // This accessor is used inside the innermost thermal-function expansions.
  // Keep the implementation in the header so optimized callers can remove
  // the call boundary (and, for the fixed orders used there, the redundant
  // bounds check) without changing coefficient values or summation order.
  double GetCoefficentAtOrder(int n) const
  {
    if (n <= MaxOrderToSave)
    {
      // All callers pass a non-negative order.  The upper bound is already
      // established by the branch above; operator[] avoids repeating the
      // vector's lower/upper bounds checks in the thermal-function loop.
      return PreCalculatedCoefficents[static_cast<std::size_t>(n)];
    }
    else
    {
      return Calculater(n);
    }
  }

  // Coefficients are populated once in the constructor and are immutable
  // afterwards.  Hot callers with a fixed expansion order can cache this
  // pointer outside their summation loop.
  const double *GetPreCalculatedCoefficentsData() const noexcept
  {
    return PreCalculatedCoefficents.data();
  }

private:
  /**
   * @brief Calculater stores the function to calculate the coefficient
   */
  std::function<double(int)> Calculater;
  /**
   * @brief MaxOrderToSave defines the highest order until which the coefficents
   * are precalculated
   */
  const int MaxOrderToSave;

  /**
   * @brief PreCalculatedCoefficents Stores the precalculated coefficents. This
   * will only be changed during the constructor to remain thread-safe!
   */
  // Orders 0..MaxOrderToSave are dense.  Contiguous storage avoids a tree
  // lookup in every term of the very hot thermal-function expansions while
  // retaining the exact precomputed double values.
  std::vector<double> PreCalculatedCoefficents;
};

} // namespace ThermalFunctions
} // namespace BSMPT

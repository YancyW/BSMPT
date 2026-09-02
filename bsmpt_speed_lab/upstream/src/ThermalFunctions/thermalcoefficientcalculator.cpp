// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/ThermalFunctions/thermalcoefficientcalculator.h>

namespace BSMPT
{
namespace ThermalFunctions
{

ThermalCoefficientCalculator::ThermalCoefficientCalculator(
    std::function<double(int)> func,
    int maxOrderToPrecalc)
    : Calculater{func}
    , MaxOrderToSave{maxOrderToPrecalc}
{
  PreCalculatedCoefficents.reserve(MaxOrderToSave + 1);
  for (int i{0}; i <= MaxOrderToSave; ++i)
  {
    PreCalculatedCoefficents.push_back(Calculater(i));
  }
}

} // namespace ThermalFunctions
} // namespace BSMPT

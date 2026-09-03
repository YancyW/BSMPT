// Copyright (C) 2020  Philipp Basler, Margarete Mühlleitner and Jonas Müller
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

/**
 * @file
 */

#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/ThermalFunctions/thermalcoefficientcalculator.h>
#include <BSMPT/bounce_solution/calcgw_profiler.h>
#include <BSMPT/models/SMparam.h>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <map>

#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>

namespace BSMPT
{
namespace ThermalFunctions
{

namespace
{
void ProfileThermalRepeat(bool fermion, double x, int diff) noexcept
{
  if (!CalcGWProfiler::thermal_repeat_enabled()) return;
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(x));
  std::memcpy(&bits, &x, sizeof(bits));
  struct LastKey
  {
    std::uint64_t bits = 0;
    int diff           = 0;
    bool fermion       = false;
    bool valid         = false;
  };
  thread_local LastKey last;
  const bool hit = last.valid && last.bits == bits && last.diff == diff &&
                   last.fermion == fermion;
  CalcGWProfiler::thermal_repeat_call(fermion, diff, hit);
  last = {bits, diff, fermion, true};
}

ThermalCoefficientCalculator FermionInterpolatedLowCoefficientCalculator(
    [](int l) -> double
    {
      if (l < 2)
      {
        return 0;
      }
      else
      {
        return gsl_sf_doublefact(2 * l - 3) * gsl_sf_zeta(2 * l - 1) /
               (gsl_sf_doublefact(2 * l) * (l + 1)) * (pow(2, 2 * l - 1) - 1);
      }
    },
    5);

ThermalCoefficientCalculator BosonInterpolatedLowCoefficientCalculator(
    [](int l) -> double
    {
      if (l < 2)
      {
        return 0;
      }
      else
      {
        return gsl_sf_doublefact(2 * l - 3) * gsl_sf_zeta(2 * l - 1) /
               (gsl_sf_doublefact(2 * l) * (l + 1));
      }
    },
    5);

ThermalCoefficientCalculator JInterpolatedHighCoefficientCalculator(
    [](int l) -> double
    {
      return 1 / (std::pow(2, l) * gsl_sf_fact(l)) * gsl_sf_gamma(2.5 + l) /
             gsl_sf_gamma(2.5 - l);
    },
    5);

bool UseFastThermalPowers()
{
  static const bool enabled = []
  {
    const char *env = std::getenv("BSMPT_USE_THERMAL_FAST_POWERS");
    return env != nullptr && env[0] == '1';
  }();
  return enabled;
}

bool UseFastFermionPowers()
{
  static const bool enabled = []
  {
    const char *env = std::getenv("BSMPT_USE_FERMION_FAST_POWERS");
    return env != nullptr && env[0] == '1';
  }();
  return UseFastThermalPowers() || enabled;
}

bool UseFastBosonPowers()
{
  static const bool enabled = []
  {
    const char *env = std::getenv("BSMPT_USE_BOSON_FAST_POWERS");
    return env != nullptr && env[0] == '1';
  }();
  return UseFastThermalPowers() || enabled;
}

bool UseFastThermalHighPowers()
{
  static const bool enabled = []
  {
    const char *env = std::getenv("BSMPT_USE_THERMAL_FAST_HIGH_POWERS");
    return env != nullptr && env[0] == '1';
  }();
  return enabled;
}

bool UseThermalCoefficientData()
{
  static const bool enabled = []
  {
    const char *env = std::getenv("BSMPT_USE_THERMAL_COEFF_DATA");
    return env != nullptr && env[0] == '1';
  }();
  return enabled;
}

bool UseThermalExactLastCache()
{
  static const bool enabled = []
  {
    const char *env = std::getenv("BSMPT_USE_THERMAL_EXACT_LAST_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  return enabled;
}

struct ThermalLastValue
{
  std::uint64_t bits = 0;
  int diff           = 0;
  double value       = 0;
  bool valid         = false;
};

std::uint64_t DoubleBits(double value) noexcept
{
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}
} // namespace

double JfermionInterpolatedLow(const double &x, const int &n, int diff)
{
  if (x == 0 and diff == 0)
  {
    return -7 * pow(M_PI, 4) / 360.0;
  }
  else if (x == 0 and diff == 1)
  {
    return pow(M_PI, 2) / 24;
  }
  using std::log;
  using std::pow;
  double res = 0;
  double cf  = 1.5 + 2 * log(4 * M_PI) - 2 * C_euler_gamma - 2 * log(4);
  if (diff == 0)
  {
    res = -7 * pow(M_PI, 4) / 360.0;
    res += pow(M_PI, 2) / 24 * x;
    res += 1 / 32.0 * pow(x, 2) * (log(x) - cf);
    double sum = 0;
    if (n <= 5)
    {
      const double *coefficients =
          FermionInterpolatedLowCoefficientCalculator
              .GetPreCalculatedCoefficentsData();
      for (int l = 2; l <= n; l++)
      {
        double Kl = coefficients[l];
        sum += pow(-x / (4 * pow(M_PI, 2)), l) * Kl;
      }
    }
    else
    {
      for (int l = 2; l <= n; l++)
      {
        double Kl = FermionInterpolatedLowCoefficientCalculator
                        .GetCoefficentAtOrder(l);
        sum += pow(-x / (4 * pow(M_PI, 2)), l) * Kl;
      }
    }
    res += -pow(M_PI, 2) * x * sum;
  }
  else if (diff == 1)
  {
    res = pow(M_PI, 2) / 24.0;
    res += x * (-6 * cf + 3) / 96;
    res += x * log(x) / 16;
    double sum = 0;
    if (n <= 5)
    {
      const double *coefficients =
          FermionInterpolatedLowCoefficientCalculator
              .GetPreCalculatedCoefficentsData();
      for (int l = 2; l <= n; l++)
      {
        double Kl = coefficients[l];
        sum += -Kl * pow(-x / 4.0, l) * (l + 1) * pow(M_PI, 2 - 2 * l);
      }
    }
    else
    {
      for (int l = 2; l <= n; l++)
      {
        double Kl = FermionInterpolatedLowCoefficientCalculator
                        .GetCoefficentAtOrder(l);
        sum += -Kl * pow(-x / 4.0, l) * (l + 1) * pow(M_PI, 2 - 2 * l);
      }
    }
    res += sum;
  }

  return res;
}

namespace
{
double JfermionInterpolatedLow4Exact(const double &x, int diff)
{
  if (x == 0 and diff == 0)
    return -7 * pow(M_PI, 4) / 360.0;
  if (x == 0 and diff == 1)
    return pow(M_PI, 2) / 24;

  using std::log;
  using std::pow;
  double res = 0;
  const double cf =
      1.5 + 2 * log(4 * M_PI) - 2 * C_euler_gamma - 2 * log(4);
  const double *coefficients =
      FermionInterpolatedLowCoefficientCalculator
          .GetPreCalculatedCoefficentsData();
  if (UseFastFermionPowers())
  {
    const double x2 = x * x;
    if (diff == 0)
    {
      const double ratio  = -x / (4 * M_PI * M_PI);
      const double ratio2 = ratio * ratio;
      const double ratio3 = ratio2 * ratio;
      const double ratio4 = ratio3 * ratio;
      res                 = -7 * pow(M_PI, 4) / 360.0;
      res += M_PI * M_PI / 24 * x;
      res += x2 * (log(x) - cf) / 32.0;
      double sum = ratio2 * coefficients[2];
      sum += ratio3 * coefficients[3];
      sum += ratio4 * coefficients[4];
      res += -M_PI * M_PI * x * sum;
    }
    else if (diff == 1)
    {
      const double ratio  = -x / 4.0;
      const double ratio2 = ratio * ratio;
      const double ratio3 = ratio2 * ratio;
      const double ratio4 = ratio3 * ratio;
      res                 = M_PI * M_PI / 24.0;
      res += x * (-6 * cf + 3) / 96;
      res += x * log(x) / 16;
      double sum = -coefficients[2] * ratio2 * 3 / (M_PI * M_PI);
      sum += -coefficients[3] * ratio3 * 4 / pow(M_PI, 4);
      sum += -coefficients[4] * ratio4 * 5 / pow(M_PI, 6);
      res += sum;
    }
    return res;
  }
  if (diff == 0)
  {
    res = -7 * pow(M_PI, 4) / 360.0;
    res += pow(M_PI, 2) / 24 * x;
    res += 1 / 32.0 * pow(x, 2) * (log(x) - cf);
    double sum = 0;
    sum += pow(-x / (4 * pow(M_PI, 2)), 2) * coefficients[2];
    sum += pow(-x / (4 * pow(M_PI, 2)), 3) * coefficients[3];
    sum += pow(-x / (4 * pow(M_PI, 2)), 4) * coefficients[4];
    res += -pow(M_PI, 2) * x * sum;
  }
  else if (diff == 1)
  {
    res = pow(M_PI, 2) / 24.0;
    res += x * (-6 * cf + 3) / 96;
    res += x * log(x) / 16;
    double sum = 0;
    sum += -coefficients[2] * pow(-x / 4.0, 2) * 3 * pow(M_PI, -2);
    sum += -coefficients[3] * pow(-x / 4.0, 3) * 4 * pow(M_PI, -4);
    sum += -coefficients[4] * pow(-x / 4.0, 4) * 5 * pow(M_PI, -6);
    res += sum;
  }
  return res;
}
} // namespace

double JbosonInterpolatedLow(const double &x, const int &n, int diff)
{
  if (x == 0 and diff == 0)
  {
    return -pow(M_PI, 4) / 45.0;
  }
  else if (x == 0 and diff == 1)
  {
    return pow(M_PI, 2) / 12.0;
  }
  using std::log;
  using std::pow;
  using std::sqrt;
  double cb  = 1.5 + 2 * std::log(4 * M_PI) - 2 * C_euler_gamma;
  double res = 0;
  if (UseFastBosonPowers() && n == 3)
  {
    const double x2      = x * x;
    const double sqrt_x  = sqrt(x);
    const double ratio   = -x / (4 * M_PI * M_PI);
    const double ratio2  = ratio * ratio;
    const double ratio3  = ratio2 * ratio;
    double coeff2, coeff3;
    if (UseThermalCoefficientData())
    {
      const double *coefficients =
          BosonInterpolatedLowCoefficientCalculator
              .GetPreCalculatedCoefficentsData();
      coeff2 = coefficients[2];
      coeff3 = coefficients[3];
    }
    else
    {
      coeff2 =
          BosonInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(2);
      coeff3 =
          BosonInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(3);
    }
    if (diff == 0)
    {
      res = -pow(M_PI, 4) / 45.0;
      res += M_PI * M_PI * x / 12.0;
      res += -M_PI * x * sqrt_x / 6;
      res += -x2 * (log(x) - cb) / 32.0;
      const double sum = ratio2 * coeff2 + ratio3 * coeff3;
      res += M_PI * M_PI * x * sum;
    }
    else if (diff == 1)
    {
      const double derivative_ratio = -x / 4.0;
      const double derivative_ratio2 = derivative_ratio * derivative_ratio;
      const double derivative_ratio3 =
          derivative_ratio2 * derivative_ratio;
      res = M_PI * M_PI / 12.0;
      res += x * (6 * cb - 3) / 96.0;
      res += -x * log(x) / 16.0;
      res += -M_PI * sqrt_x / 4.0;
      res += coeff2 * derivative_ratio2 * 3 / (M_PI * M_PI);
      res += coeff3 * derivative_ratio3 * 4 / pow(M_PI, 4);
    }
    return res;
  }
  if (diff == 0)
  {
    res = -pow(M_PI, 4) / 45.0;
    res += pow(M_PI, 2) * x / 12.0;
    res += -M_PI * pow(x, 1.5) / 6;
    res += -pow(x, 2) * (log(x) - cb) / 32.0;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      double Kl =
          BosonInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += pow(-x / (4 * pow(M_PI, 2)), l) * Kl;
    }
    res += pow(M_PI, 2) * x * sum;
  }
  else if (diff == 1)
  {
    res = pow(M_PI, 2) / 12.0;
    res += x * (6 * cb - 3) / 96.0;
    res += -x * log(x) / 16.0;
    res += -M_PI * sqrt(x) / 4.0;
    double sum = 0;
    for (int l = 2; l <= n; l++)
    {
      double Kl =
          BosonInterpolatedLowCoefficientCalculator.GetCoefficentAtOrder(l);
      sum += Kl * pow(-x / 4.0, l) * (l + 1) * pow(M_PI, 2 - 2 * l);
    }
    res += sum;
  }
  return res;
}

double JInterpolatedHigh(const double &x, const int &n, int diff)
{
  using std::exp;
  using std::pow;
  using std::sqrt;

  double res = 0;
  if (UseFastThermalHighPowers() && n == 3)
  {
    const double sqrt_x = sqrt(x);
    const double inv_sqrt_x = 1.0 / sqrt_x;
    const double inv_x = inv_sqrt_x * inv_sqrt_x;
    const double *coefficients =
        JInterpolatedHighCoefficientCalculator
            .GetPreCalculatedCoefficentsData();
    if (diff == 0)
    {
      double sum = coefficients[0];
      sum += coefficients[1] * inv_sqrt_x;
      sum += coefficients[2] * inv_x;
      sum += coefficients[3] * inv_x * inv_sqrt_x;
      res = -exp(-sqrt_x) * sqrt(M_PI / 2 * x * sqrt_x) * sum;
    }
    else if (diff == 1)
    {
      double sum = coefficients[0] * sqrt_x * (2 * sqrt_x - 3);
      sum += coefficients[1] * (2 + 2 * sqrt_x - 3);
      sum += coefficients[2] * inv_sqrt_x * (4 + 2 * sqrt_x - 3);
      sum += coefficients[3] * inv_x * (6 + 2 * sqrt_x - 3);
      const double x_three_quarters = sqrt(x * sqrt_x);
      res = exp(-sqrt_x) * sqrt(2 * M_PI) /
            (8 * x_three_quarters) * sum;
    }
    return res;
  }
  if (diff == 0)
  {
    double sum = 0;
    if (n <= 5)
    {
      const double *coefficients =
          JInterpolatedHighCoefficientCalculator
              .GetPreCalculatedCoefficentsData();
      for (int l = 0; l <= n; l++)
      {
        double Kl = coefficients[l];
        sum += Kl * pow(x, -l / 2.0);
      }
    }
    else
    {
      for (int l = 0; l <= n; l++)
      {
        double Kl = JInterpolatedHighCoefficientCalculator
                        .GetCoefficentAtOrder(l);
        sum += Kl * pow(x, -l / 2.0);
      }
    }
    res = -exp(-sqrt(x)) * sqrt(M_PI / 2 * pow(x, 1.5)) * sum;
  }
  else if (diff == 1)
  {
    double sum = 0;
    if (n <= 5)
    {
      const double *coefficients =
          JInterpolatedHighCoefficientCalculator
              .GetPreCalculatedCoefficentsData();
      for (int l = 0; l <= n; l++)
      {
        double Kl = coefficients[l];
        sum += Kl * pow(x, (1 - l) / 2) * (2 * l + 2 * sqrt(x) - 3);
      }
    }
    else
    {
      for (int l = 0; l <= n; l++)
      {
        double Kl = JInterpolatedHighCoefficientCalculator
                        .GetCoefficentAtOrder(l);
        sum += Kl * pow(x, (1 - l) / 2) * (2 * l + 2 * sqrt(x) - 3);
      }
    }
    res = exp(-sqrt(x)) * sqrt(2 * M_PI) / (8 * pow(x, 3.0 / 4.0)) * sum;
  }
  return res;
}

double JfermionInterpolated(const double &x, int diff)
{
  ProfileThermalRepeat(true, x, diff);
  thread_local ThermalLastValue cache;
  const bool use_cache = UseThermalExactLastCache();
  const std::uint64_t bits = use_cache ? DoubleBits(x) : 0;
  if (use_cache && cache.valid && cache.bits == bits && cache.diff == diff)
    return cache.value;
  static const bool UseExactLow4 = []
  {
    const char *env = std::getenv("BSMPT_USE_FERMION_LOW4_EXACT");
    return env != nullptr && env[0] == '1';
  }();
  double res = 0;
  if (x >= C_FermionTheta)
  {
    res = -JInterpolatedHigh(x, 3, diff);
  }
  else
  {
    res = -(UseExactLow4 ? JfermionInterpolatedLow4Exact(x, diff)
                         : JfermionInterpolatedLow(x, 4, diff));
    if (diff == 0) res += -C_FermionShift;
  }
  if (use_cache) cache = {bits, diff, res, true};
  return res;
}

double JbosonInterpolated(const double &x, int diff)
{
  ProfileThermalRepeat(false, x, diff);
  thread_local ThermalLastValue cache;
  const bool use_cache = UseThermalExactLastCache();
  const std::uint64_t bits = use_cache ? DoubleBits(x) : 0;
  if (use_cache && cache.valid && cache.bits == bits && cache.diff == diff)
    return cache.value;
  double res = 0;
  if (x >= C_BosonTheta)
  {
    res = JInterpolatedHigh(x, 3, diff);
    if (diff == 0) res -= C_BosonShift;
  }
  else if (x >= 0)
  {
    res = JbosonInterpolatedLow(x, 3, diff);
  }
  else if (x < 0)
  {
    res = JbosonInterpolatedNegative(x, diff);
  }
  if (use_cache) cache = {bits, diff, res, true};
  return res;
}

double JbosonInterpolatedNegative(const double &x, int diff)
{
  if (x >= 0) return 0;
  double PotVal = 0;

  if (diff == 0)
  {
    PotVal = JbosonNegativeSpline(-x) + 3.533375127e-06;
  }
  else if (diff == 1)
  {
    PotVal = -JbosonNegativeSpline.deriv(1, -x);
  }

  return PotVal;
}

} // namespace ThermalFunctions
} // namespace BSMPT

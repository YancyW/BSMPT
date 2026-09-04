// Copyright (C) 2018  Philipp Basler and Margarete Mühlleitner
// SPDX-FileCopyrightText: 2021 Philipp Basler, Margarete Mühlleitner and Jonas
// Müller
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gsl/gsl_sf_gamma.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <random>

#include "Eigen/Dense"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_zeta.h>
#include <cstdint>
#include <cstring>
#include <unordered_map>

#include <BSMPT/ThermalFunctions/NegativeBosonSpline.h>
#include <BSMPT/ThermalFunctions/ThermalFunctions.h>
#include <BSMPT/bounce_solution/calcgw_profiler.h>
#include <BSMPT/models/ClassPotentialOrigin.h>
#include <BSMPT/models/IncludeAllModels.h>
#include <BSMPT/utility/Logger.h>
#include <BSMPT/utility/utility.h>
using namespace Eigen;

namespace BSMPT
{
namespace
{
uint64_t ExactDoubleBits(double value) noexcept
{
  uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

void ProfileExactModelRepeat(CalcGWProfiler::ExactRepeatMetric metric,
                             const void *model,
                             const std::vector<double> &v,
                             double temperature,
                             int diff,
                             int order)
{
  if (!CalcGWProfiler::exact_repeat_enabled()) return;
  struct LastKey
  {
    const void *model = nullptr;
    std::vector<std::uint64_t> fields;
    std::uint64_t temperature = 0;
    int diff                 = 0;
    int order                = 0;
    bool valid               = false;
  };
  thread_local std::array<LastKey, 3> last;
  const auto index = static_cast<std::size_t>(metric);
  auto &key        = last.at(index);
  bool hit = key.valid && key.model == model && key.diff == diff &&
             key.order == order &&
             key.temperature == ExactDoubleBits(temperature) &&
             key.fields.size() == v.size();
  if (hit)
  {
    for (std::size_t i = 0; i < v.size(); ++i)
      if (key.fields[i] != ExactDoubleBits(v[i]))
      {
        hit = false;
        break;
      }
  }
  key.model       = model;
  key.temperature = ExactDoubleBits(temperature);
  key.diff        = diff;
  key.order       = order;
  key.valid       = true;
  key.fields.resize(v.size());
  for (std::size_t i = 0; i < v.size(); ++i)
    key.fields[i] = ExactDoubleBits(v[i]);
  CalcGWProfiler::exact_repeat_call(metric, hit);
}

struct BosonContributionCacheKey
{
  uint64_t massSquared = 0;
  uint64_t temperature = 0;
  uint64_t cb          = 0;
  uint64_t devMassSquared = 0;
  int derivative       = 0;

  bool operator==(const BosonContributionCacheKey &other) const noexcept
  {
    return massSquared == other.massSquared &&
           temperature == other.temperature && cb == other.cb &&
           devMassSquared == other.devMassSquared &&
           derivative == other.derivative;
  }
};

struct BosonContributionCacheKeyHash
{
  std::size_t operator()(const BosonContributionCacheKey &key) const noexcept
  {
    std::size_t hash = static_cast<std::size_t>(key.massSquared);
    hash ^= static_cast<std::size_t>(key.temperature) +
            static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) +
            (hash >> 2);
    hash ^= static_cast<std::size_t>(key.cb) +
            static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) +
            (hash >> 2);
    hash ^= static_cast<std::size_t>(key.devMassSquared) +
            static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) +
            (hash >> 2);
    hash ^= static_cast<std::size_t>(key.derivative) +
            static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) +
            (hash >> 2);
    return hash;
  }
};

struct FermionContributionCacheKey
{
  uint64_t massSquared = 0;
  uint64_t temperature  = 0;
  int derivative        = 0;

  bool operator==(const FermionContributionCacheKey &other) const noexcept
  {
    return massSquared == other.massSquared &&
           temperature == other.temperature && derivative == other.derivative;
  }
};

struct FermionContributionCacheKeyHash
{
  std::size_t operator()(const FermionContributionCacheKey &key) const noexcept
  {
    std::size_t hash = static_cast<std::size_t>(key.massSquared);
    hash ^= static_cast<std::size_t>(key.temperature) +
            static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) +
            (hash >> 2);
    hash ^= static_cast<std::size_t>(key.derivative) +
            static_cast<std::size_t>(0x9e3779b97f4a7c15ULL) + (hash << 6) +
            (hash >> 2);
    return hash;
  }
};

// Real symmetric counterpart of FirstDerivativeOfEigenvalues.  Higgs and
// gauge mass matrices are real, so an opt-in caller can avoid converting them
// to MatrixXcd while retaining the same 1e-6 degeneracy prescription.
std::vector<double>
FirstDerivativeOfRealEigenvalues(const MatrixXd &M, const MatrixXd &MDiff)
{
  const Eigen::Index nRows = M.rows();
  const Eigen::Index nCols = M.cols();
  const double EVThres    = std::pow(10, -6);
  if (nRows != nCols || MDiff.rows() != nRows || MDiff.cols() != nCols)
    throw std::runtime_error(
        "ERROR ! M and MDiff need to be equally sized quadratic matrices !\n");
  const std::size_t nSize = static_cast<std::size_t>(nRows);

  SelfAdjointEigenSolver<MatrixXd> es(M);
  if (es.info() != Success)
    return {};

  const auto eigenvalues = es.eigenvalues();
  const auto eigenvectors = es.eigenvectors();
  std::vector<double> values(nSize);
  std::vector<double> derivatives(nSize);
  for (std::size_t i = 0; i < nSize; ++i)
  {
    values[i] = eigenvalues[i];
    if (std::abs(values[i]) < EVThres) values[i] = 0;
  }

  // Eigen sorts the spectrum, so equal-within-threshold eigenvalues form
  // contiguous groups.  This is equivalent to the pairwise Mapping used by
  // Class_Potential_Origin::FirstDerivativeOfEigenvalues.
  for (std::size_t first = 0; first < nSize;)
  {
    std::size_t last = first + 1;
    while (last < nSize &&
           std::abs(values[first] - values[last]) <= EVThres)
      ++last;

    const std::size_t groupSize = last - first;
    if (groupSize == 1)
    {
      derivatives[first] =
          (eigenvectors.col(first).transpose() * MDiff *
           eigenvectors.col(first))(0, 0);
    }
    else
    {
      const MatrixXd phi = eigenvectors.block(0, first, nRows, groupSize);
      const MatrixXd projected = phi.transpose() * MDiff * phi;
      SelfAdjointEigenSolver<MatrixXd> groupSolver(projected,
                                                    EigenvaluesOnly);
      if (groupSolver.info() != Success) return {};
      for (std::size_t i = 0; i < groupSize; ++i)
        derivatives[first + i] = groupSolver.eigenvalues()[i];
    }
    first = last;
  }

  std::vector<double> result;
  result.reserve(2 * nRows);
  for (double value : values) result.push_back(value);
  for (double derivative : derivatives)
  {
    if (std::abs(derivative) < EVThres) derivative = 0;
    result.push_back(derivative);
  }
  return result;
}

} // namespace

namespace
{
// Per-V1Loop constants.  The powers and denominator are deliberately formed
// with the same std::pow expressions as boson()/fermion(); the cached callers
// below only replace repeated evaluation across particles.
struct V1LoopThermalConstants
{
  double temperatureSquared;
  double temperatureCubed;
  double temperatureFourth;
  double piDenominator;
};

bool UseV1LoopThermalContext() noexcept
{
  static const bool enabled = []
  {
    const char *env = std::getenv("BSMPT_USE_V1LOOP_THERMAL_CONTEXT");
    return env != nullptr && env[0] == '1';
  }();
  return enabled;
}

double BosonWithV1LoopThermalContext(
    const Class_Potential_Origin &model,
    double MassSquared,
    double Temp,
    double cb,
    int diff,
    double DevMassSquared,
    const V1LoopThermalConstants &constants)
{
  double res = 0;
  if (diff >= 0) res = model.CWTerm(MassSquared, cb, diff);
  const double Ratio = MassSquared / constants.temperatureSquared;
  if (diff == 0)
  {
    res += constants.temperatureFourth / constants.piDenominator *
           ThermalFunctions::JbosonInterpolated(Ratio);
  }
  else if (diff > 0)
  {
    res += constants.temperatureSquared / constants.piDenominator *
           ThermalFunctions::JbosonInterpolated(Ratio, 1);
  }
  else if (diff == -1)
  {
    res += 1.0 / constants.piDenominator *
           (4 * constants.temperatureCubed *
                ThermalFunctions::JbosonInterpolated(Ratio, 0) -
            (2 * Temp * MassSquared - Temp * Temp * DevMassSquared) *
                ThermalFunctions::JbosonInterpolated(Ratio, 1));
  }
  return res;
}

double FermionWithV1LoopThermalContext(
    const Class_Potential_Origin &model,
    double MassSquared,
    double Temp,
    int diff,
    const V1LoopThermalConstants &constants)
{
  double res = 0;
  if (diff >= 0) res = model.CWTerm(std::abs(MassSquared), C_CWcbFermion, diff);
  const double Ratio = MassSquared / constants.temperatureSquared;
  if (diff == 0)
  {
    res += constants.temperatureFourth / constants.piDenominator *
           ThermalFunctions::JfermionInterpolated(Ratio);
  }
  else if (diff > 0)
  {
    res += constants.temperatureSquared / constants.piDenominator *
           ThermalFunctions::JfermionInterpolated(Ratio, 1);
  }
  else if (diff == -1)
  {
    res += 1.0 / constants.piDenominator *
           (4 * constants.temperatureCubed *
                ThermalFunctions::JfermionInterpolated(Ratio, 0) -
            2 * Temp * MassSquared *
                ThermalFunctions::JfermionInterpolated(Ratio, 1));
  }
  return res;
}
} // namespace

Class_Potential_Origin::Class_Potential_Origin()
    : Class_Potential_Origin(GetSMConstants())
{
}

Class_Potential_Origin::Class_Potential_Origin(const ISMConstants &smConstants)
    : SMConstants{smConstants}
    , scale{SMConstants.C_vev0}

{
  // TODO Auto-generated constructor stub
}

Class_Potential_Origin::~Class_Potential_Origin()
{
  // TODO Auto-generated destructor stub
}

/**
 * This will call set_gen(par), SetCurvatureArrays, set_CT_Pot_Par(parCT),
 * CalculateDebye() as well as CalculateDebyeGauge()
 */
void Class_Potential_Origin::set_All(const std::vector<double> &par,
                                     const std::vector<double> &parCT)
{

  set_gen(par);
  if (!SetCurvatureDone) SetCurvatureArrays();
  set_CT_Pot_Par(parCT);
  CalculateDebye();
  CalculateDebyeGauge();
}

void Class_Potential_Origin::Prepare_Triple()
{
  for (std::size_t a = 0; a < NHiggs; a++)
  {
    for (std::size_t b = 0; b < NHiggs; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        LambdaHiggs_3_CT[a][b][i] = Curvature_Higgs_CT_L3[a][b][i];
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          LambdaHiggs_3_CT[a][b][i] +=
              Curvature_Higgs_CT_L4[a][b][i][j] * HiggsVev[j];
        }
      }
    }
  }
}

double Class_Potential_Origin::FCW(double MassSquared) const
{
  double res = 0;
  double x;
  double Boarder = std::pow(10, -200);
  if (std::isnan(MassSquared))
    x = Boarder;
  else if (std::abs(MassSquared) < Boarder)
    x = Boarder;
  else
    x = std::abs(MassSquared);
  res = std::log(x) - 2 * std::log(scale);
  return res;
}

double
Class_Potential_Origin::CWTerm(double MassSquared, double cb, int diff) const
{

  if (std::abs(MassSquared) < C_threshold) return 0;
  double LogTerm = 0, PotVal = 0;
  LogTerm = FCW(MassSquared);
  if (diff == 0)
    PotVal =
        1.0 / (64 * M_PI * M_PI) * MassSquared * MassSquared * (LogTerm - cb);
  else if (diff > 0)
  {
    PotVal = 1.0 / (32 * M_PI * M_PI) * MassSquared * (LogTerm - cb + 0.5);
  }
  return PotVal;
}

double Class_Potential_Origin::boson(double MassSquared,
                                     double Temp,
                                     double cb,
                                     int diff,
                                     double DevMassSquared) const
{
  double res = 0;
  if (diff >= 0) res = CWTerm(MassSquared, cb, diff);
  if (Temp == 0) return res;
  double Ratio = MassSquared / std::pow(Temp, 2);
  if (diff == 0)
  {
    res += std::pow(Temp, 4) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JbosonInterpolated(Ratio);
  }
  else if (diff > 0)
  {
    res += std::pow(Temp, 2) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JbosonInterpolated(Ratio, 1);
  }
  else if (diff == -1)
  {
    res += 1.0 / (2 * std::pow(M_PI, 2)) *
           (4 * std::pow(Temp, 3) *
                ThermalFunctions::JbosonInterpolated(Ratio, 0) -
            (2 * Temp * MassSquared - Temp * Temp * DevMassSquared) *
                ThermalFunctions::JbosonInterpolated(Ratio, 1));
  }
  return res;
}

double
Class_Potential_Origin::fermion(double MassSquared, double Temp, int diff) const
{
  double res = 0;
  if (diff >= 0) res = CWTerm(std::abs(MassSquared), C_CWcbFermion, diff);
  double Ratio = MassSquared / std::pow(Temp, 2);
  if (Temp == 0) return res;
  if (diff == 0)
  {
    res += std::pow(Temp, 4) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JfermionInterpolated(Ratio);
  }

  else if (diff > 0)
  {
    res += std::pow(Temp, 2) / (2 * std::pow(M_PI, 2)) *
           ThermalFunctions::JfermionInterpolated(Ratio, 1);
  }
  else if (diff == -1)
  {
    res += 1.0 / (2 * std::pow(M_PI, 2)) *
           (4 * std::pow(Temp, 3) *
                ThermalFunctions::JfermionInterpolated(Ratio, 0) -
            2 * Temp * MassSquared *
                ThermalFunctions::JfermionInterpolated(Ratio, 1));
  }
  return res;
}

std::vector<double> Class_Potential_Origin::FirstDerivativeOfEigenvalues(
    const Ref<MatrixXcd> M,
    const Ref<MatrixXcd> MDiff) const
{
  std::vector<double> res;
  const std::size_t nRows = M.rows();
  const std::size_t nCols = M.cols();

  // Every successful path appends n eigenvalues followed by n derivatives.
  // Reserve once to avoid the geometric-growth reallocations in this hot
  // helper; this does not change the values or their summation order.
  res.reserve(2 * nRows);

  const double EVThres = std::pow(10, -6);

  if (nCols != nRows)
  {
    throw std::runtime_error("ERROR ! M needs to be an quadratic Matrix for "
                             "calculating the derivatives !\n");
  }

  const std::size_t nSize = nRows;

  SelfAdjointEigenSolver<MatrixXcd> es;
  es.compute(M);

  std::vector<std::complex<double>> Eigenvalues(nSize);
  std::vector<std::complex<double>> Derivatives(nSize);
  std::vector<double> AlreadyCalculated(
      nSize); // Array to check which EVs already been calculated.
  for (std::size_t i = 0; i < nSize; i++)
  {
    Eigenvalues[i] = es.eigenvalues()[i];
    if (std::abs(Eigenvalues[i]) < EVThres)
    {
      Eigenvalues[i] = 0;
    }
  }

  std::vector<std::vector<double>> Mapping(nSize, std::vector<double>(nSize));

  for (std::size_t i = 0; i < nSize; i++)
  {
    AlreadyCalculated[i] = -1;
    for (std::size_t j = i; j < nSize; j++)
    {
      if (std::abs(Eigenvalues[i] - Eigenvalues[j]) > EVThres)
      {
        Mapping[i][j] = 0;
      }
      else
      {
        Mapping[i][j] = 1;
      }
    }
  }
  for (std::size_t i = 1; i < nSize; i++)
  {
    for (std::size_t j = 0; j < i; j++)
      Mapping[i][j] = Mapping[j][i];
  }

  for (std::size_t p = 0; p < nSize; p++)
  {
    if (AlreadyCalculated[p] == -1)
    {
      std::size_t NumOfReps = 0;
      for (std::size_t i = p + 1; i < nSize; i++)
      {
        NumOfReps += Mapping[p][i];
      }
      if (NumOfReps == 0)
      {
        VectorXcd v(nSize);
        v = es.eigenvectors().col(p);
        // Derivatives[p] = (v.transpose()*MDiff*v).value();
        Derivatives[p]       = (v.adjoint() * MDiff * v).value();
        AlreadyCalculated[p] = 1;
      }
      else
      {
        MatrixXcd Phi(nSize, NumOfReps + 1);
        std::size_t helpCol = 0;
        MatrixXcd MXWork(NumOfReps + 1, NumOfReps + 1);

        for (std::size_t i = p; i < nSize; i++)
        {
          if (Mapping[p][i] == 1)
          {
            Phi.col(helpCol) = es.eigenvectors().col(i);
            helpCol++;
          }
        }
        MXWork = Phi.adjoint() * MDiff * Phi;
        SelfAdjointEigenSolver<MatrixXcd> esWork;
        esWork.compute(MXWork);
        helpCol = 0;
        for (std::size_t i = p; i < nSize; i++)
        {
          if (Mapping[p][i] == 1)
          {
            AlreadyCalculated[i] = 1;
            Derivatives[i]       = esWork.eigenvalues()[helpCol];
            helpCol++;
          }
        }
      }
    }
  }
  for (std::size_t i = 0; i < nSize; i++)
  {
    if (std::abs(Derivatives[i]) < EVThres) Derivatives[i] = 0;
  }

  for (std::size_t i = 0; i < nSize; i++)
  {
    res.push_back(Eigenvalues[i].real());
  }
  for (std::size_t i = 0; i < nSize; i++)
    res.push_back(Derivatives[i].real());
  return res;
}

double Class_Potential_Origin::fbaseTri(double MassSquaredA,
                                        double MassSquaredB,
                                        double MassSquaredC) const
{
  double res  = 0;
  double mas  = MassSquaredA;
  double mbs  = MassSquaredB;
  double mcs  = MassSquaredC;
  double LogA = 0, LogB = 0, LogC = 0;
  double Thres = 1e-8;
  if (std::abs(mas) < Thres) mas = 0;
  if (std::abs(mbs) < Thres) mbs = 0;
  if (std::abs(mcs) < Thres) mcs = 0;
  if (std::abs(mas - mbs) < Thres) mas = mbs;
  if (std::abs(mas - mcs) < Thres) mas = mcs;
  if (std::abs(mbs - mcs) < Thres) mbs = mcs;

  if (mas != 0) LogA = std::log(mas) - 2 * std::log(scale);
  if (mbs != 0) LogB = std::log(mbs) - 2 * std::log(scale);
  if (mcs != 0) LogC = std::log(mcs) - 2 * std::log(scale);

  std::size_t C = 1;
  if (mas == 0 and mbs == 0 and mcs == 0)
    res = 0;
  else if (mas != 0 and mbs == 0 and mcs == 0)
  {
    C   = 2;
    res = 1.0 / mas * (LogA - 1);
  }
  else if (mas == 0 and mbs != 0 and mcs == 0)
  {
    C   = 3;
    res = (LogB - 1) / mbs;
  }
  else if (mas == 0 and mbs == 0 and mcs != 0)
  {
    C   = 4;
    res = (LogC - 1) / mcs;
  }
  else if (mas == mbs and mas != 0 and mas != mcs and mcs != 0)
  {
    C   = 6;
    res = (mbs - mcs + mcs * std::log(mcs / mbs)) / std::pow(mbs - mcs, 2);
  }
  else if (mas == mcs and mas != 0 and mas != mbs and mbs != 0)
  {
    C   = 7;
    res = (mbs * log(mbs / mcs) - mbs + mcs) / std::pow(mbs - mcs, 2);
  }
  else if (mbs == mcs and mas != 0 and mbs != mas and mbs != 0)
  {
    C   = 8;
    res = (mas * std::log(mas / mcs) - mas + mcs) / std::pow(mas - mcs, 2);
  }
  else if (mas == mbs and mas == mcs and mas != 0)
  {
    C   = 9;
    res = 1.0 / (2 * mcs);
  }
  else if (mas == mbs and mas != mcs and mas != 0 and mcs == 0)
  {
    C   = 10;
    res = 1.0 / mbs;
  }
  else if (mas == mcs and mas != mbs and mas != 0 and mbs == 0)
  {
    C   = 11;
    res = 1.0 / mas;
  }
  else if (mbs == mcs and mbs != 0 and mbs != mas and mas == 0)
  {
    C   = 12;
    res = 1.0 / mbs;
  }
  else
  {
    C   = 5;
    res = mas * LogA / ((mas - mbs) * (mas - mcs)) +
          mbs * LogB / ((mbs - mas) * (mbs - mcs));
    res += mcs * LogC / ((mcs - mas) * (mcs - mbs));
  }

  if (std::isnan(res) or std::isinf(res))
  {
    std::string throwstring = "Found nan at line = ";
    throwstring += std::to_string(InputLineNumber);
    throwstring += " in function ";
    throwstring += __func__;
    throwstring += "\n";
    std::stringstream ss;
    ss << "Found nan at line = " << InputLineNumber << " in function "
       << __func__ << std::endl;
    ss << mas << sep << mbs << sep << mcs << sep << res << sep << C
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str(), __FILE__, __LINE__);
    throw std::runtime_error(throwstring.c_str());
  }

  return res;
}

double Class_Potential_Origin::fbaseFour(double MassSquaredA,
                                         double MassSquaredB,
                                         double MassSquaredC,
                                         double MassSquaredD) const
{

  double res  = 0;
  double mas  = MassSquaredA;
  double mbs  = MassSquaredB;
  double mcs  = MassSquaredC;
  double mds  = MassSquaredD;
  double LogA = 0, LogB = 0, LogC = 0, LogD = 0;
  double thresZero = 1e-6;
  double thresDeg  = 1e-6;
  if (std::abs(mas) < thresZero) mas = 0;
  if (std::abs(mbs) < thresZero) mbs = 0;
  if (std::abs(mcs) < thresZero) mcs = 0;
  if (std::abs(mds) < thresZero) mds = 0;
  if (std::abs(mas - mbs) < thresDeg) mbs = mas;
  if (std::abs(mas - mcs) < thresDeg) mcs = mas;
  if (std::abs(mas - mds) < thresDeg) mds = mas;
  if (std::abs(mbs - mcs) < thresDeg) mcs = mbs;
  if (std::abs(mbs - mds) < thresDeg) mds = mbs;
  if (std::abs(mcs - mds) < thresDeg) mds = mcs;

  if (mas != 0) LogA = std::log(mas) - 2 * std::log(scale);
  if (mbs != 0) LogB = std::log(mbs) - 2 * std::log(scale);
  if (mcs != 0) LogC = std::log(mcs) - 2 * std::log(scale);
  if (mds != 0) LogD = std::log(mds) - 2 * std::log(scale);

  std::size_t C = 0;

  // all masses are zero
  if (mas == 0 and mbs == 0 and mcs == 0 and mds == 0)
  {
    C = 1;
    // f0000
    res = 0;
  }

  // one mass is non-zero, the other ones are zero
  else if (mas != 0 and mbs == 0 and mcs == 0 and mds == 0)
  {
    C = 2;
    // fa000
    res = LogA / (mas * mas);
  }
  else if (mas == 0 and mbs != 0 and mcs == 0 and mds == 0)
  {
    C = 3;
    // f0b00
    res = LogB / (mbs * mbs);
  }
  else if (mas == 0 and mbs == 0 and mcs != 0 and mds == 0)
  {
    C = 4;
    // f00c0
    res = LogC / (mcs * mcs);
  }
  else if (mas == 0 and mbs == 0 and mcs == 0 and mds != 0)
  {
    C = 5;
    // f000d
    res = LogD / (mds * mds);
  }

  // two masses are non-zero, the other masses are zero
  // 1) they are equal
  else if (mas == mbs and mas != 0 and mcs == 0 and mds == 0)
  {
    C = 6;
    // faa00
    res = (2. - LogA) / (mas * mas);
  }
  else if (mas == mcs and mas != 0 and mbs == 0 and mds == 0)
  {
    C = 7;
    // fa0a0
    res = (2. - LogA) / (mas * mas);
  }
  else if (mas == mds and mas != 0 and mbs == 0 and mcs == 0)
  {
    C = 8;
    // fa00a
    res = (2. - LogA) / (mas * mas);
  }
  else if (mbs == mcs and mbs != 0 and mas == 0 and mds == 0)
  {
    C = 9;
    // f0bb0
    res = (2. - LogB) / (mbs * mbs);
  }
  else if (mbs == mds and mbs != 0 and mas == 0 and mcs == 0)
  {
    C = 10;
    // f0b0b
    res = (2. - LogB) / (mbs * mbs);
  }
  else if (mcs == mds and mcs != 0 and mas == 0 and mbs == 0)
  {
    C = 11;
    // f00cc
    res = (2. - LogC) / (mcs * mcs);
  }

  // 2) they are non-equal
  else if (mas != 0 and mbs != 0 and mas != mbs and mcs == 0 and mds == 0)
  {
    C = 12;
    // fab00
    res = 1. / (mas * mbs) +
          (mbs * LogA - mas * LogB) / (mas * mbs * (mas - mbs));
  }
  else if (mas != 0 and mcs != 0 and mas != mcs and mbs == 0 and mds == 0)
  {
    C = 13;
    // fa0c0
    res = 1. / (mas * mcs) +
          (mcs * LogA - mas * LogC) / (mas * mcs * (mas - mcs));
  }
  else if (mas != 0 and mds != 0 and mas != mds and mbs == 0 and mcs == 0)
  {
    C = 14;
    // fa00d
    res = 1. / (mas * mds) +
          (mds * LogA - mas * LogD) / (mas * mds * (mas - mds));
  }
  else if (mbs != 0 and mcs != 0 and mbs != mcs and mas == 0 and mds == 0)
  {
    C = 15;
    // f0bc0
    res = 1. / (mbs * mcs) +
          (mcs * LogB - mbs * LogC) / (mbs * mcs * (mbs - mcs));
  }
  else if (mbs != 0 and mds != 0 and mbs != mds and mas == 0 and mcs == 0)
  {
    C = 16;
    // f0b0d
    res = 1. / (mbs * mds) +
          (mds * LogB - mbs * LogD) / (mbs * mds * (mbs - mds));
  }
  else if (mcs != 0 and mds != 0 and mcs != mds and mas == 0 and mbs == 0)
  {
    C = 17;
    // f00cd
    res = 1. / (mcs * mds) +
          (mds * LogC - mcs * LogD) / (mcs * mds * (mcs - mds));
  }

  // three masses are non-zero, the remaining mass is zero
  // 1) the three masses are equal
  else if (mas == 0 and mbs == mcs and mbs == mds and mbs != 0)
  {
    C = 18;
    // f0bbb
    res = -1. / (2. * mbs * mbs);
  }
  else if (mbs == 0 and mas == mcs and mas == mds and mas != 0)
  {
    C = 19;
    // fa0aa
    res = -1. / (2. * mas * mas);
  }
  else if (mcs == 0 and mas == mbs and mas == mds and mas != 0)
  {
    C = 20;
    // faa0a
    res = -1. / (2. * mas * mas);
  }
  else if (mds == 0 and mas == mbs and mas == mcs and mas != 0)
  {
    C = 21;
    // faaa0
    res = -1. / (2. * mas * mas);
  }

  // 2) two of the three non-zero masses are equal
  else if (mas == 0 and mbs != 0 and mcs != 0 and mds != 0)
  {
    if (mbs == mcs and mds != mbs)
    {
      C = 22;
      // f0bbd
      res = (mbs - mds + mbs * LogD - mbs * LogB) /
            (mbs * (mbs - mds) * (mbs - mds));
    }
    else if (mbs == mds and mcs != mbs)
    {
      C = 23;
      // f0bcb
      res = (mbs - mcs + mbs * LogC - mbs * LogB) /
            (mbs * (mbs - mcs) * (mbs - mcs));
    }
    else if (mcs == mds and mbs != mcs)
    {
      C = 24;
      // f0bcc
      res = (mcs - mbs + mcs * LogB - mcs * LogC) /
            (mcs * (mbs - mcs) * (mbs - mcs));
    }
  }
  else if (mas != 0 and mbs == 0 and mcs != 0 and mds != 0)
  {
    if (mas == mcs and mds != mas)
    {
      C = 25;
      // fa0ad
      res = (mas - mds + mas * LogD - mas * LogA) /
            (mas * (mas - mds) * (mas - mds));
    }
    else if (mas == mds and mcs != mas)
    {
      C = 26;
      // fa0ca
      res = (mas - mcs + mas * LogC - mas * LogA) /
            (mas * (mas - mcs) * (mas - mcs));
    }
    else if (mcs == mds and mas != mcs)
    {
      C = 27;
      // fa0cc
      res = (mcs - mas + mcs * LogA - mcs * LogC) /
            (mcs * (mcs - mas) * (mcs - mas));
    }
  }
  else if (mas != 0 and mbs != 0 and mcs == 0 and mds != 0)
  {
    if (mas == mbs and mds != mas)
    {
      C = 28;
      // faa0d
      res = (mas - mds + mas * LogD - mas * LogA) /
            (mas * (mas - mds) * (mas - mds));
    }
    else if (mas == mds and mbs != mas)
    {
      C = 29;
      // fab0a
      res = (mas - mbs + mas * LogB - mas * LogA) /
            (mas * (mas - mbs) * (mas - mbs));
    }
    else if (mbs == mds and mas != mbs)
    {
      C = 30;
      // fab0b
      res = (mbs - mas + mbs * LogA - mbs * LogB) /
            (mbs * (mas - mbs) * (mas - mbs));
    }
  }
  else if (mas != 0 and mbs != 0 and mcs != 0 and mds == 0)
  {
    if (mas == mbs and mcs != mas)
    {
      C = 31;
      // faac0
      res = (mas - mcs + mas * LogC - mas * LogA) /
            (mas * (mas - mcs) * (mas - mcs));
    }
    else if (mas == mcs and mbs != mas)
    {
      C = 32;
      // faba0
      res = (mas - mbs + mas * LogB - mas * LogA) /
            (mas * (mas - mbs) * (mas - mbs));
    }
    else if (mbs == mcs and mas != mbs)
    {
      C = 33;
      // fabb0
      res = (mbs - mas + mbs * LogA - mbs * LogB) /
            (mbs * (mas - mbs) * (mas - mbs));
    }
  }

  // 3) all three non-zero masses are different
  else if (mas == 0 and mbs != mcs and mbs != mds and mcs != mds and
           mbs != 0 and mcs != 0 and mds != 0)
  {
    C = 34;
    // f0bcd
    res = LogB / (mbs - mcs) / (mbs - mds) + LogC / (mcs - mbs) / (mcs - mds) +
          LogD / (mds - mbs) / (mds - mcs);
  }
  else if (mbs == 0 and mas != mcs and mas != mds and mcs != mds and
           mas != 0 and mcs != 0 and mds != 0)
  {
    C = 35;
    // fa0cd
    res = LogA / (mas - mcs) / (mas - mds) + LogC / (mcs - mas) / (mcs - mds) +
          LogD / (mds - mas) / (mds - mcs);
  }
  else if (mcs == 0 and mas != mbs and mas != mds and mbs != mds and
           mas != 0 and mbs != 0 and mds != 0)
  {
    C = 36;
    // fab0d
    res = LogA / (mas - mbs) / (mas - mds) + LogB / (mbs - mas) / (mbs - mds) +
          LogD / (mds - mas) / (mds - mbs);
  }
  else if (mds == 0 and mas != mbs and mas != mcs and mbs != mcs and
           mas != 0 and mbs != 0 and mcs != 0)
  {
    C = 37;
    // fabc0
    res = LogA / (mas - mbs) / (mas - mcs) + LogB / (mbs - mas) / (mbs - mcs) +
          LogC / (mcs - mas) / (mcs - mbs);
  }

  // all four masses are non-zero
  // 1) all four masses are equal
  else if (mas == mbs and mbs == mcs and mcs == mds and mas != 0)
  {
    C = 38;
    // faaaa
    res = -1 / (6 * mas * mas);
  }
  // 2) only two of the masses are equal
  // 2.1) remaining two are not equal
  else if (mas == mbs and mcs != mds and mas != mcs and mas != mds and
           mas != 0 and mcs != 0 and mds != 0)
  {
    C = 39;
    // faacd
    res = (mcs * (mas - mds) * (mas - mds) * LogC -
           mds * (mas - mcs) * (mas - mcs) * LogD +
           (mcs - mds) *
               ((mas - mcs) * (mas - mds) + (-mas * mas + mcs * mds) * LogA)) /
          ((mas - mcs) * (mas - mcs) * (mas - mds) * (mas - mds) * (mcs - mds));
  }
  else if (mas == mcs and mbs != mds and mas != mbs and mas != mds and
           mas != 0 and mbs != 0 and mds != 0)
  {
    C = 40;
    // fabad
    res = (mbs * (mas - mds) * (mas - mds) * LogB -
           mds * (mas - mbs) * (mas - mbs) * LogD +
           (mbs - mds) *
               ((mas - mbs) * (mas - mds) + (-mas * mas + mbs * mds) * LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mds) * (mas - mds) * (mbs - mds));
  }
  else if (mas == mds and mbs != mcs and mas != mbs and mas != mcs and
           mas != 0 and mbs != 0 and mcs != 0)
  {
    C = 41;
    // fabca
    res = (mbs * (mas - mcs) * (mas - mcs) * LogB -
           mcs * (mas - mbs) * (mas - mbs) * LogC +
           (mbs - mcs) *
               ((mas - mbs) * (mas - mcs) + (-mas * mas + mbs * mcs) * LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mcs) * (mas - mcs) * (mbs - mcs));
  }
  else if (mbs == mcs and mas != mds and mbs != mas and mbs != mds and
           mbs != 0 and mas != 0 and mds != 0)
  {
    C = 42;
    // fabbd
    res = (mas * (mbs - mds) * (mbs - mds) * LogA -
           mds * (mas - mbs) * (mas - mbs) * LogD +
           (mas - mds) *
               (-(mas - mbs) * (mbs - mds) + (-mbs * mbs + mas * mds) * LogB)) /
          ((mas - mbs) * (mas - mbs) * (mbs - mds) * (mbs - mds) * (mas - mds));
  }
  else if (mbs == mds and mas != mcs and mbs != mas and mbs != mcs and
           mbs != 0 and mas != 0 and mcs != 0)
  {
    C = 43;
    // fabcb
    res = (mas * (mbs - mcs) * (mbs - mcs) * LogA -
           mcs * (mas - mbs) * (mas - mbs) * LogC +
           (mas - mcs) *
               (-(mas - mbs) * (mbs - mcs) + (-mbs * mbs + mas * mcs) * LogB)) /
          ((mas - mbs) * (mas - mbs) * (mbs - mcs) * (mbs - mcs) * (mas - mcs));
  }
  else if (mcs == mds and mas != mbs and mcs != mas and mcs != mbs and
           mcs != 0 and mcs != 0 and mds != 0)
  {
    C = 44;
    // fabcc
    res = (mas * (mbs - mcs) * (mbs - mcs) * LogA -
           mbs * (mas - mcs) * (mas - mcs) * LogB +
           (mas - mbs) *
               ((mas - mcs) * (mbs - mcs) + (-mcs * mcs + mas * mbs) * LogC)) /
          ((mas - mcs) * (mas - mcs) * (mbs - mcs) * (mbs - mcs) * (mas - mbs));
  }
  // 2.2) remaining two are also equal
  else if (mas == mbs and mcs == mds and mas != mcs and mas != 0 and mcs != 0)
  {
    C = 45;
    // faacc
    res = (2. * (mas - mcs) + (mas + mcs) * (LogC - LogA)) /
          ((mas - mcs) * (mas - mcs) * (mas - mcs));
  }
  else if (mas == mcs and mbs == mds and mas != mbs and mas != 0 and mbs != 0)
  {
    C = 46;
    // fabab
    res = (2. * (mas - mbs) + (mas + mbs) * (LogB - LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mbs));
  }
  else if (mas == mds and mbs == mcs and mas != mbs and mas != 0 and mbs != 0)
  {
    C = 47;
    // fabba
    res = (2. * (mas - mbs) + (mas + mbs) * (LogB - LogA)) /
          ((mas - mbs) * (mas - mbs) * (mas - mbs));
  }

  // 3) three of the masses are equal
  else if (mas == mbs and mas == mcs and mas != 0 and mas != mds and mds != 0)
  {
    C = 48;
    // faaad
    res = (-mas * mas + mds * mds + 2. * mas * mds * (LogA - LogD)) /
          (2. * mas * (mas - mds) * (mas - mds) * (mas - mds));
  }
  else if (mas == mbs and mas == mds and mas != 0 and mas != mcs and mcs != 0)
  {
    C = 49;
    // faaca
    res = (-mas * mas + mcs * mcs + 2. * mas * mcs * (LogA - LogC)) /
          (2. * mas * (mas - mcs) * (mas - mcs) * (mas - mcs));
  }
  else if (mas == mcs and mas == mds and mas != 0 and mas != mbs and mbs != 0)
  {
    C = 50;
    // fabaa
    res = (-mas * mas + mbs * mbs + 2. * mas * mbs * (LogA - LogB)) /
          (2. * mas * (mas - mbs) * (mas - mbs) * (mas - mbs));
  }
  else if (mbs == mcs and mbs == mds and mbs != 0 and mas != mbs and mas != 0)
  {
    C = 51;
    // fabbb
    res = (-mbs * mbs + mas * mas + 2. * mbs * mas * (LogB - LogA)) /
          (2. * mbs * (mbs - mas) * (mbs - mas) * (mbs - mas));
  }

  // 4) all four masses are non-equal
  else
  {
    C = 52;
    // fabcd
    res = mas * LogA / ((mas - mbs) * (mas - mcs) * (mas - mds));
    res += mbs * LogB / ((mbs - mas) * (mbs - mcs) * (mbs - mds));
    res += mcs * LogC / ((mcs - mas) * (mcs - mbs) * (mcs - mds));
    res += mds * LogD / ((mds - mas) * (mds - mbs) * (mds - mcs));
  }

  if (std::isnan(res) or std::isinf(res))
  {
    std::string throwstring = "Found nan at line = ";
    throwstring += std::to_string(InputLineNumber);
    throwstring += " in function ";
    throwstring += __func__;
    throwstring += "\n";
    std::stringstream ss;
    // maximize the numerical precision
    typedef std::numeric_limits<double> dbl;
    ss << std::setprecision(dbl::max_digits10);
    ss << "Found nan at line = " << InputLineNumber << " in function "
       << __func__ << std::endl;
    ss << mas << sep << mbs << sep << mcs << sep << mds << sep << C << sep
       << res << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str(), __FILE__, __LINE__);
    throw std::runtime_error(throwstring.c_str());
  }

  return res;
}

double Class_Potential_Origin::fbase(double MassSquaredA,
                                     double MassSquaredB) const
{
  double res  = 0;
  double LogA = 0;
  if (MassSquaredA == 0 and MassSquaredB == 0) return 1;
  double ZB = std::pow(10, -5);
  if (MassSquaredA != 0) LogA = std::log(MassSquaredA) - 2 * std::log(scale);
  if (std::abs(MassSquaredA - MassSquaredB) > ZB)
  {
    double LogB = 0;
    if (MassSquaredB != 0) LogB = std::log(MassSquaredB) - 2 * std::log(scale);
    if (MassSquaredA == 0)
      res = LogB;
    else if (MassSquaredB == 0)
      res = LogA;
    else
      res = (LogA * MassSquaredA - LogB * MassSquaredB) /
            (MassSquaredA - MassSquaredB);
  }
  else
  {
    res = 1 + LogA;
  }
  return res;
}

std::vector<double>
Class_Potential_Origin::SecondDerivativeOfEigenvaluesNonRepeated(
    const Eigen::Ref<Eigen::MatrixXd> M,
    const Eigen::Ref<Eigen::MatrixXd> MDiffX,
    const Eigen::Ref<Eigen::MatrixXd> MDiffY,
    const Eigen::Ref<Eigen::MatrixXd> MDiffXY) const
{
  std::vector<double> res;
  const std::size_t nRows = M.rows();
  const std::size_t nCols = M.cols();

  const std::size_t EVThres = std::pow(10, -6);

  if (nCols != nRows)
  {
    throw std::runtime_error("ERROR ! M needs to be an quadratic Matrix for "
                             "calculating the derivatives !\n");
  }

  const std::size_t nSize = nRows;

  SelfAdjointEigenSolver<MatrixXd> es;
  es.compute(M);

  std::vector<double> Eigenvalues(nSize);
  for (std::size_t i = 0; i < nSize; i++)
    Eigenvalues[i] = es.eigenvalues()[i];
  for (std::size_t i = 0; i < nSize - 1; i++)
  {
    if (std::abs(Eigenvalues[i] - Eigenvalues[i + 1]) < EVThres)
    {
      Logger::Write(LoggingLevel::Default, "ERROR ! repeated eigenvalues.");
    }
  }

  std::vector<std::vector<double>> Deriv(nSize, std::vector<double>(4));
  VectorXd v(nSize);
  MatrixXd C(nSize, nSize), E(nSize, nSize), Identity(nSize, nSize);
  Identity = MatrixXd::Identity(nSize, nSize);
  VectorXd vDiffX(nSize), vDiffY(nSize);

  for (std::size_t i = 0; i < nSize; i++)
  {
    Deriv[i][0] = Eigenvalues[i];
    v           = es.eigenvectors().col(i);
    Deriv[i][1] = v.transpose() * MDiffX * v;
    Deriv[i][2] = v.transpose() * MDiffY * v;

    C = (M - Deriv[i][0] * Identity).transpose() *
            (M - Deriv[i][0] * Identity) +
        v * v.transpose();
    E = (M - Deriv[i][0] * Identity).transpose() *
        (MDiffX - Deriv[i][1] * Identity);

    vDiffX = C.colPivHouseholderQr().solve(-E * v);

    E = (M - Deriv[i][0] * Identity).transpose() *
        (MDiffY - Deriv[i][2] * Identity);
    vDiffY = C.colPivHouseholderQr().solve(-E * v);

    Deriv[i][3] = v.transpose() * MDiffXY * v;
    Deriv[i][3] += v.transpose() * (MDiffX - Deriv[i][1] * Identity) * vDiffY;
    Deriv[i][3] += v.transpose() * (MDiffY - Deriv[i][2] * Identity) * vDiffX;
  }

  for (const auto &x : Eigenvalues)
    res.push_back(x);

  for (std::size_t i = 0; i < nSize; i++)
  {
    for (std::size_t j = 0; j < 4; j++)
      res.push_back(Deriv[i][j]);
  }
  return res;
}

// Sanity check to make sure HiggsRotationMatrix is a proper rotation
// matrix, i.e. its inverse should correspond to its transpose, and its
// determinant should be +1 or -1
bool Class_Potential_Origin::CheckRotationMatrix()
{
  MatrixXd mat(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      mat(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  const double det_precision = 1e-10; // det precision
  const double el_precision =
      1e-8; // element precision. different precision to prevent numerical
            // instabilities for small elements

  if (!almost_the_same(std::abs(mat.determinant()), 1., det_precision))
  {
    return false;
  }

  auto inv    = mat.inverse();
  auto transp = mat.transpose();

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (!almost_the_same(inv(i, j), transp(i, j), el_precision))
      {
        return false;
      }
    }
  }

  return true;
}

void Class_Potential_Origin::CalculatePhysicalCouplings()
{
  if (!SetCurvatureDone) SetCurvatureArrays();
  const double ZeroMass = std::pow(10, -5);
  MatrixXd MassHiggs(NHiggs, NHiggs), MassGauge(NGauge, NGauge);
  MatrixXcd MassQuark(NQuarks, NQuarks), MassLepton(NLepton, NLepton);
  MassHiggs  = MatrixXd::Zero(NHiggs, NHiggs);
  MassGauge  = MatrixXd::Zero(NGauge, NGauge);
  MassQuark  = MatrixXcd::Zero(NQuarks, NQuarks);
  MassLepton = MatrixXcd::Zero(NLepton, NLepton);

  MassSquaredGauge.resize(NGauge);
  MassSquaredHiggs.resize(NHiggs);
  MassSquaredQuark.resize(NQuarks);
  MassSquaredLepton.resize(NLepton);
  HiggsRotationMatrix.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
    HiggsRotationMatrix[i].resize(NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      MassHiggs(i, j) += Curvature_Higgs_L2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        MassHiggs(i, j) += Curvature_Higgs_L3[i][j][k] * HiggsVev[k];
        for (std::size_t l = 0; l < NHiggs; l++)
          MassHiggs(i, j) +=
              0.5 * Curvature_Higgs_L4[i][j][k][l] * HiggsVev[k] * HiggsVev[l];
      }
    }
  }

  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          MassGauge(a, b) += 0.5 * Curvature_Gauge_G2H2[a][b][i][j] *
                             HiggsVev[i] * HiggsVev[j];
        }
      }
    }
  }

  MatrixXcd MIJQuarks = QuarkMassMatrix(HiggsVev);

  MassQuark = MIJQuarks.conjugate() * MIJQuarks;

  MatrixXcd MIJLeptons = LeptonMassMatrix(HiggsVev);

  MassLepton = MIJLeptons.conjugate() * MIJLeptons;

  MatrixXd HiggsRot(NHiggs, NHiggs), GaugeRot(NGauge, NGauge),
      QuarkRot(NQuarks, NQuarks), LepRot(NLepton, NLepton);
  HiggsRot = MatrixXd::Identity(NHiggs, NHiggs);
  GaugeRot = MatrixXd::Identity(NGauge, NGauge);
  QuarkRot = MatrixXd::Identity(NQuarks, NQuarks);
  LepRot   = MatrixXd::Identity(NLepton, NLepton);

  SelfAdjointEigenSolver<MatrixXd> es;

  es.compute(MassHiggs);
  HiggsRot = es.eigenvectors().transpose();
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(HiggsRot(i, j)) < std::pow(10, -10)) HiggsRot(i, j) = 0;
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    MassSquaredHiggs[i] = es.eigenvalues()[i];
    if (std::abs(MassSquaredHiggs[i]) < ZeroMass) MassSquaredHiggs[i] = 0;
  }

  es.compute(MassGauge);
  GaugeRot = es.eigenvectors().transpose();

  for (std::size_t i = 0; i < NGauge; i++)
  {
    MassSquaredGauge[i] = es.eigenvalues()[i];
    if (std::abs(MassSquaredGauge[i]) < ZeroMass) MassSquaredGauge[i] = 0;
  }

  SelfAdjointEigenSolver<MatrixXcd> esQuark(MassQuark);
  QuarkRot = esQuark.eigenvectors().transpose().real();
  for (std::size_t i = 0; i < NQuarks; i++)
    MassSquaredQuark[i] = esQuark.eigenvalues().real()[i];

  SelfAdjointEigenSolver<MatrixXcd> esLepton(MassLepton);
  LepRot = esLepton.eigenvectors().transpose().real();
  for (std::size_t i = 0; i < NLepton; i++)
    MassSquaredLepton[i] = esLepton.eigenvalues().real()[i];

  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        LambdaGauge_3[a][b][i] = 0;
        for (std::size_t j = 0; j < NHiggs; j++)
          LambdaGauge_3[a][b][i] +=
              Curvature_Gauge_G2H2[a][b][i][j] * HiggsVev[j];
      }
    }
  }
  for (std::size_t a = 0; a < NHiggs; a++)
  {
    for (std::size_t b = 0; b < NHiggs; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        LambdaHiggs_3[a][b][i] = Curvature_Higgs_L3[a][b][i];

        for (std::size_t j = 0; j < NHiggs; j++)
        {
          LambdaHiggs_3[a][b][i] +=
              Curvature_Higgs_L4[a][b][i][j] * HiggsVev[j];
        }
      }
    }
  }

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < NQuarks; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        LambdaQuark_3[i][j][k] = 0;
        for (std::size_t l = 0; l < NQuarks; l++)
        {
          LambdaQuark_3[i][j][k] +=
              conj(Curvature_Quark_F2H1[i][l][k]) * MIJQuarks(l, j);
          LambdaQuark_3[i][j][k] +=
              conj(MIJQuarks(i, l)) * Curvature_Quark_F2H1[l][j][k];
        }
        for (std::size_t m = 0; m < NHiggs; m++)
        {
          LambdaQuark_4[i][j][k][m] = 0;
          for (std::size_t l = 0; l < NQuarks; l++)
          {
            LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][k]) *
                                         Curvature_Quark_F2H1[l][j][m];
            LambdaQuark_4[i][j][k][m] += conj(Curvature_Quark_F2H1[i][l][m]) *
                                         Curvature_Quark_F2H1[l][j][k];
          }
        }
      }
    }
  }

  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < NLepton; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        LambdaLepton_3[i][j][k] = 0;
        for (std::size_t l = 0; l < NLepton; l++)
        {
          LambdaLepton_3[i][j][k] +=
              conj(Curvature_Lepton_F2H1[i][l][k]) * MIJLeptons(l, j);
          LambdaLepton_3[i][j][k] +=
              conj(MIJLeptons(i, l)) * Curvature_Lepton_F2H1[l][j][k];
        }
        for (std::size_t m = 0; m < NHiggs; m++)
        {
          LambdaLepton_4[i][j][k][m] = 0;
          for (std::size_t l = 0; l < NLepton; l++)
          {
            LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][k]) *
                                          Curvature_Lepton_F2H1[l][j][m];
            LambdaLepton_4[i][j][k][m] += conj(Curvature_Lepton_F2H1[i][l][m]) *
                                          Curvature_Lepton_F2H1[l][j][k];
          }
        }
      }
    }
  }

  // Rotate and save std::size_to corresponding vectors

  Couplings_Higgs_Quartic.resize(NHiggs);
  Couplings_Higgs_Triple.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    Couplings_Higgs_Quartic[i].resize(NHiggs);
    Couplings_Higgs_Triple[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      Couplings_Higgs_Quartic[i][j].resize(NHiggs);
      Couplings_Higgs_Triple[i][j].resize(NHiggs);
      for (std::size_t k = 0; k < NHiggs; k++)
        Couplings_Higgs_Quartic[i][j][k].resize(NHiggs);
    }
  }

  Couplings_Gauge_Higgs_22.resize(NGauge);
  Couplings_Gauge_Higgs_21.resize(NGauge);
  for (std::size_t a = 0; a < NGauge; a++)
  {
    Couplings_Gauge_Higgs_22[a].resize(NGauge);
    Couplings_Gauge_Higgs_21[a].resize(NGauge);
    for (std::size_t b = 0; b < NGauge; b++)
    {
      Couplings_Gauge_Higgs_22[a][b].resize(NHiggs);
      Couplings_Gauge_Higgs_21[a][b].resize(NHiggs);
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Gauge_Higgs_22[a][b][i].resize(NHiggs);
      }
    }
  }

  Couplings_Quark_Higgs_22.resize(NQuarks);
  Couplings_Quark_Higgs_21.resize(NQuarks);
  for (std::size_t a = 0; a < NQuarks; a++)
  {
    Couplings_Quark_Higgs_22[a].resize(NQuarks);
    Couplings_Quark_Higgs_21[a].resize(NQuarks);
    for (std::size_t b = 0; b < NQuarks; b++)
    {
      Couplings_Quark_Higgs_22[a][b].resize(NHiggs);
      Couplings_Quark_Higgs_21[a][b].resize(NHiggs);
      for (std::size_t i = 0; i < NHiggs; i++)
        Couplings_Quark_Higgs_22[a][b][i].resize(NHiggs);
    }
  }

  Couplings_Lepton_Higgs_22.resize(NLepton);
  Couplings_Lepton_Higgs_21.resize(NLepton);
  for (std::size_t a = 0; a < NLepton; a++)
  {
    Couplings_Lepton_Higgs_22[a].resize(NLepton);
    Couplings_Lepton_Higgs_21[a].resize(NLepton);
    for (std::size_t b = 0; b < NLepton; b++)
    {
      Couplings_Lepton_Higgs_22[a][b].resize(NHiggs);
      Couplings_Lepton_Higgs_21[a][b].resize(NHiggs);
      for (std::size_t i = 0; i < NHiggs; i++)
        Couplings_Lepton_Higgs_22[a][b][i].resize(NHiggs);
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        Couplings_Higgs_Triple[i][j][k] = 0;
        for (std::size_t is = 0; is < NHiggs; is++)
        {
          for (std::size_t js = 0; js < NHiggs; js++)
          {
            for (std::size_t ks = 0; ks < NHiggs; ks++)
            {
              Couplings_Higgs_Triple[i][j][k] +=
                  HiggsRot(i, is) * HiggsRot(j, js) * HiggsRot(k, ks) *
                  LambdaHiggs_3[is][js][ks];
            }
          }
        }
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          Couplings_Higgs_Quartic[i][j][k][l] = 0;
          for (std::size_t is = 0; is < NHiggs; is++)
          {
            for (std::size_t js = 0; js < NHiggs; js++)
            {
              for (std::size_t ks = 0; ks < NHiggs; ks++)
              {
                for (std::size_t ls = 0; ls < NHiggs; ls++)
                {
                  Couplings_Higgs_Quartic[i][j][k][l] +=
                      HiggsRot(i, is) * HiggsRot(j, js) * HiggsRot(k, ks) *
                      HiggsRot(l, ls) * Curvature_Higgs_L4[is][js][ks][ls];
                }
              }
            }
          }
        }
      }
    }
  }

  // Gauge Rot
  for (std::size_t a = 0; a < NGauge; a++)
  {
    for (std::size_t b = 0; b < NGauge; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Gauge_Higgs_21[a][b][i] = 0;
        for (std::size_t as = 0; as < NGauge; as++)
        {
          for (std::size_t bs = 0; bs < NGauge; bs++)
          {
            for (std::size_t is = 0; is < NHiggs; is++)
              Couplings_Gauge_Higgs_21[a][b][i] +=
                  GaugeRot(a, as) * GaugeRot(b, bs) * HiggsRot(i, is) *
                  LambdaGauge_3[as][bs][is];
          }
        }
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          Couplings_Gauge_Higgs_22[a][b][i][j] = 0;
          for (std::size_t as = 0; as < NGauge; as++)
          {
            for (std::size_t bs = 0; bs < NGauge; bs++)
            {
              for (std::size_t is = 0; is < NHiggs; is++)
              {
                for (std::size_t js = 0; js < NHiggs; js++)
                {
                  double RotFac = GaugeRot(a, as) * GaugeRot(b, bs) *
                                  HiggsRot(i, is) * HiggsRot(j, js);
                  Couplings_Gauge_Higgs_22[a][b][i][j] +=
                      RotFac * Curvature_Gauge_G2H2[as][bs][is][js];
                }
              }
            }
          }
        }
      }
    }
  }

  // Quark

  for (std::size_t a = 0; a < NQuarks; a++)
  {
    for (std::size_t b = 0; b < NQuarks; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Quark_Higgs_21[a][b][i] = 0;
        for (std::size_t as = 0; as < NQuarks; as++)
        {
          for (std::size_t bs = 0; bs < NQuarks; bs++)
          {
            for (std::size_t is = 0; is < NHiggs; is++)
            {
              double RotFac =
                  QuarkRot(a, as) * QuarkRot(b, bs) * HiggsRot(i, is);
              Couplings_Quark_Higgs_21[a][b][i] +=
                  RotFac * LambdaQuark_3[as][bs][is];
            }
          }
        }
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          Couplings_Quark_Higgs_22[a][b][i][j] = 0;
          for (std::size_t as = 0; as < NQuarks; as++)
          {
            for (std::size_t bs = 0; bs < NQuarks; bs++)
            {
              for (std::size_t is = 0; is < NHiggs; is++)
              {
                for (std::size_t js = 0; js < NHiggs; js++)
                {
                  double RotFac = QuarkRot(a, as) * QuarkRot(b, bs) *
                                  HiggsRot(i, is) * HiggsRot(j, js);
                  Couplings_Quark_Higgs_22[a][b][i][j] +=
                      RotFac * LambdaQuark_4[as][bs][is][js];
                }
              }
            }
          }
        }
      }
    }
  }

  // Lepton

  for (std::size_t a = 0; a < NLepton; a++)
  {
    for (std::size_t b = 0; b < NLepton; b++)
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        Couplings_Lepton_Higgs_21[a][b][i] = 0;
        for (std::size_t as = 0; as < NLepton; as++)
        {
          for (std::size_t bs = 0; bs < NLepton; bs++)
          {
            for (std::size_t is = 0; is < NHiggs; is++)
            {
              double RotFac = LepRot(a, as) * LepRot(b, bs) * HiggsRot(i, is);
              Couplings_Lepton_Higgs_21[a][b][i] +=
                  RotFac * LambdaLepton_3[as][bs][is];
            }
          }
        }
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          Couplings_Lepton_Higgs_22[a][b][i][j] = 0;
          for (std::size_t as = 0; as < NLepton; as++)
          {
            for (std::size_t bs = 0; bs < NLepton; bs++)
            {
              for (std::size_t is = 0; is < NHiggs; is++)
              {
                for (std::size_t js = 0; js < NHiggs; js++)
                {
                  double RotFac = LepRot(a, as) * LepRot(b, bs) *
                                  HiggsRot(i, is) * HiggsRot(j, js);
                  Couplings_Lepton_Higgs_22[a][b][i][j] +=
                      RotFac * LambdaLepton_4[as][bs][is][js];
                }
              }
            }
          }
        }
      }
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRotationMatrix[i][j] = HiggsRot(i, j);
    }
  }

  CalcCouplingsDone = true;

  return;
}

std::vector<double> Class_Potential_Origin::WeinbergFirstDerivative() const
{
  std::vector<double> res;
  if (!CalcCouplingsDone)
  {
    //        CalculatePhysicalCouplings();
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }
  const double NumZero = std::pow(10, -10);
  VectorXd FirstDeriv(NHiggs), FirstDerivGauge(NHiggs), FirstDerivHiggs(NHiggs),
      FirstDerivQuark(NHiggs), FirstDerivLepton(NHiggs);
  FirstDeriv       = VectorXd::Zero(NHiggs);
  FirstDerivGauge  = VectorXd::Zero(NHiggs);
  FirstDerivHiggs  = VectorXd::Zero(NHiggs);
  FirstDerivQuark  = VectorXd::Zero(NHiggs);
  FirstDerivLepton = VectorXd::Zero(NHiggs);
  double epsilon   = 1.0 / (16 * M_PI * M_PI);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t a = 0; a < NGauge; a++)
    {
      if (MassSquaredGauge[a] != 0)
      {
        FirstDerivGauge[i] +=
            MassSquaredGauge[a] * Couplings_Gauge_Higgs_21[a][a][i] *
            (std::log(MassSquaredGauge[a] / std::pow(scale, 2)) - C_CWcbGB +
             0.5);
      }
    }

    for (std::size_t a = 0; a < NHiggs; a++)
    {
      if (MassSquaredHiggs[a] != 0)
      {

        FirstDerivHiggs[i] +=
            MassSquaredHiggs[a] * Couplings_Higgs_Triple[a][a][i] *
            (std::log(MassSquaredHiggs[a] / std::pow(scale, 2)) - C_CWcbHiggs +
             0.5);
      }
    }
    for (std::size_t a = 0; a < NQuarks; a++)
    {
      if (MassSquaredQuark[a] != 0)
      {
        double Coup = Couplings_Quark_Higgs_21[a][a][i].real();
        FirstDerivQuark[i] +=
            MassSquaredQuark[a] * Coup *
            (std::log(MassSquaredQuark[a] / std::pow(scale, 2)) -
             C_CWcbFermion + 0.5);
      }
    }
    for (std::size_t a = 0; a < NLepton; a++)
    {
      if (MassSquaredLepton[a] != 0)
      {
        double Coup = Couplings_Lepton_Higgs_21[a][a][i].real();
        FirstDerivLepton[i] +=
            MassSquaredLepton[a] * Coup *
            (std::log(MassSquaredLepton[a] / std::pow(scale, 2)) -
             C_CWcbFermion + 0.5);
      }
    }
  }
  FirstDerivGauge *= 1.5;
  FirstDerivHiggs *= 0.5;
  FirstDerivQuark *= -3;
  FirstDerivLepton *= -1;

  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
  }

  FirstDeriv = HiggsRot.transpose() * (FirstDerivGauge + FirstDerivHiggs +
                                       FirstDerivQuark + FirstDerivLepton);
  FirstDeriv *= epsilon;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (std::abs(FirstDeriv[i]) < NumZero) FirstDeriv[i] = 0;
  }

  for (std::size_t i = 0; i < NHiggs; i++)
    res.push_back(FirstDeriv[i]);

  return res;
}

Eigen::MatrixXd
Class_Potential_Origin::WeinbergSecondDerivativeAsMatrixXd() const
{
  if (!CalcCouplingsDone)
  {
    //        CalculatePhysicalCouplings();
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }
  const double NumZero = std::pow(10, -10);
  MatrixXd GaugePart(NHiggs, NHiggs), HiggsPart(NHiggs, NHiggs),
      QuarkPart(NHiggs, NHiggs), LeptonPart(NHiggs, NHiggs);
  GaugePart  = MatrixXd::Zero(NHiggs, NHiggs);
  HiggsPart  = MatrixXd::Zero(NHiggs, NHiggs);
  QuarkPart  = MatrixXd::Zero(NHiggs, NHiggs);
  LeptonPart = MatrixXd::Zero(NHiggs, NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t a = 0; a < NGauge; a++)
      {
        for (std::size_t b = 0; b < NGauge; b++)
        {
          double Coup1 = Couplings_Gauge_Higgs_21[a][b][i];
          double Coup2 = Couplings_Gauge_Higgs_21[b][a][j];
          double Br =
              fbase(MassSquaredGauge[a], MassSquaredGauge[b]) - C_CWcbGB + 0.5;
          GaugePart(i, j) += Coup1 * Coup2 * Br;
        }
        if (MassSquaredGauge[a] != 0)
        {
          GaugePart(i, j) +=
              MassSquaredGauge[a] * Couplings_Gauge_Higgs_22[a][a][i][j] *
              (std::log(MassSquaredGauge[a] / std::pow(scale, 2)) - C_CWcbGB +
               0.5);
        }
      }

      for (std::size_t a = 0; a < NHiggs; a++)
      {
        for (std::size_t b = 0; b < NHiggs; b++)
        {
          double Coup1 = Couplings_Higgs_Triple[a][b][i];
          double Coup2 = Couplings_Higgs_Triple[b][a][j];
          double Br    = fbase(MassSquaredHiggs[a], MassSquaredHiggs[b]) -
                      C_CWcbHiggs + 0.5;
          HiggsPart(i, j) += Coup1 * Coup2 * Br;
        }
        if (MassSquaredHiggs[a] != 0)
        {
          HiggsPart(i, j) +=
              MassSquaredHiggs[a] * Couplings_Higgs_Quartic[a][a][i][j] *
              (std::log(MassSquaredHiggs[a] / std::pow(scale, 2)) -
               C_CWcbHiggs + 0.5);
        }
      }

      for (std::size_t a = 0; a < NQuarks; a++)
      {
        for (std::size_t b = 0; b < NQuarks; b++)
        {
          double Coup = (Couplings_Quark_Higgs_21[a][b][i] *
                         Couplings_Quark_Higgs_21[b][a][j])
                            .real();
          double Br = fbase(MassSquaredQuark[a], MassSquaredQuark[b]) -
                      C_CWcbFermion + 0.5;
          QuarkPart(i, j) += Coup * Br;
        }
        if (MassSquaredQuark[a] != 0)
        {
          double Coup = Couplings_Quark_Higgs_22[a][a][i][j].real();
          QuarkPart(i, j) +=
              MassSquaredQuark[a] * Coup *
              (std::log(MassSquaredQuark[a] / std::pow(scale, 2)) -
               C_CWcbFermion + 0.5);
        }
      }

      for (std::size_t a = 0; a < NLepton; a++)
      {
        for (std::size_t b = 0; b < NLepton; b++)
        {
          double Coup = (Couplings_Lepton_Higgs_21[a][b][i] *
                         Couplings_Lepton_Higgs_21[b][a][j])
                            .real();
          double Br = fbase(MassSquaredLepton[a], MassSquaredLepton[b]) -
                      C_CWcbFermion + 0.5;
          LeptonPart(i, j) += Coup * Br;
        }
        if (MassSquaredLepton[a] != 0)
        {
          double Coup = Couplings_Lepton_Higgs_22[a][a][i][j].real();
          LeptonPart(i, j) +=
              Coup * MassSquaredLepton[a] *
              (std::log(MassSquaredLepton[a] / std::pow(scale, 2)) -
               C_CWcbFermion + 0.5);
        }
      }
    }
  }

  HiggsPart *= 0.5;
  GaugePart *= 1.5;
  QuarkPart *= -3;
  LeptonPart *= -1;

  MatrixXd Storage(NHiggs, NHiggs);
  Storage = HiggsPart + GaugePart + QuarkPart + LeptonPart;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(Storage(i, j)) < NumZero) Storage(i, j) = 0;
    }
  }
  MatrixXd ResMatrix;
  MatrixXd HiggsRot(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      HiggsRot(i, j) = HiggsRotationMatrix[i][j];
    }
  }

  ResMatrix =
      0.5 * HiggsRot.transpose() * (Storage + Storage.transpose()) * HiggsRot;
  double epsilon = 1.0 / (16.0 * M_PI * M_PI);
  ResMatrix *= epsilon;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      if (std::abs(ResMatrix(i, j)) < NumZero) ResMatrix(i, j) = 0;
    }
  }
  return ResMatrix;
}
std::vector<double> Class_Potential_Origin::WeinbergSecondDerivative() const
{

  auto ResMatrix = WeinbergSecondDerivativeAsMatrixXd();
  std::vector<double> res;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      res.push_back(ResMatrix(j, i));
    }
  }

  return res;
}

std::vector<double> Class_Potential_Origin::WeinbergThirdDerivative() const
{

  if (not CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }
  const double NumZero = std::pow(10, -10);
  double epsilon       = 1.0 / (16.0 * M_PI * M_PI);

  std::vector<double> res;

  std::vector<std::vector<std::vector<std::complex<double>>>> restmp;
  std::vector<std::vector<std::vector<std::complex<double>>>> QuarkPart;
  std::vector<std::vector<std::vector<std::complex<double>>>> LeptonPart;
  std::vector<std::vector<std::vector<std::complex<double>>>> QuarkPartSym;
  std::vector<std::vector<std::vector<std::complex<double>>>> LeptonPartSym;
  restmp.resize(NHiggs);
  QuarkPart.resize(NHiggs);
  LeptonPart.resize(NHiggs);
  QuarkPartSym.resize(NHiggs);
  LeptonPartSym.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    restmp[i].resize(NHiggs);
    QuarkPart[i].resize(NHiggs);
    LeptonPart[i].resize(NHiggs);
    QuarkPartSym[i].resize(NHiggs);
    LeptonPartSym[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      restmp[i][j].resize(NHiggs);
      QuarkPart[i][j].resize(NHiggs);
      LeptonPart[i][j].resize(NHiggs);
      QuarkPartSym[i][j].resize(NHiggs);
      LeptonPartSym[i][j].resize(NHiggs);
    }
  }

  std::vector<std::vector<std::vector<double>>> resGaugeBase(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));
  std::vector<std::vector<std::vector<double>>> Higgspart(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  std::vector<std::vector<std::vector<double>>> GaugePart(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  std::vector<std::vector<std::vector<double>>> HiggspartSym(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  std::vector<std::vector<std::vector<double>>> GaugePartSym(
      NHiggs,
      std::vector<std::vector<double>>(NHiggs, std::vector<double>(NHiggs)));

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        Higgspart[i][j][k] = 0;
        for (std::size_t a = 0; a < NHiggs; a++)
        {
          for (std::size_t b = 0; b < NHiggs; b++)
          {
            for (std::size_t c = 0; c < NHiggs; c++)
            {
              double f1 = fbaseTri(MassSquaredHiggs[a],
                                   MassSquaredHiggs[b],
                                   MassSquaredHiggs[c]);
              double f2 = Couplings_Higgs_Triple[a][b][i];
              double f3 = Couplings_Higgs_Triple[b][c][j];
              double f4 = Couplings_Higgs_Triple[c][a][k];
              Higgspart[i][j][k] += 2 * f1 * f2 * f3 * f4;
            }
            double f1 = Couplings_Higgs_Quartic[a][b][i][j];
            double f2 = Couplings_Higgs_Triple[b][a][k];
            double f3 = fbase(MassSquaredHiggs[a], MassSquaredHiggs[b]) -
                        C_CWcbHiggs + 0.5;
            Higgspart[i][j][k] += 3 * f1 * f2 * f3;
          }
        }

        GaugePart[i][j][k] = 0;
        for (std::size_t a = 0; a < NGauge; a++)
        {
          for (std::size_t b = 0; b < NGauge; b++)
          {
            for (std::size_t c = 0; c < NGauge; c++)
            {
              double f1 = fbaseTri(MassSquaredGauge[a],
                                   MassSquaredGauge[b],
                                   MassSquaredGauge[c]);
              double f2 = Couplings_Gauge_Higgs_21[a][b][i];
              double f3 = Couplings_Gauge_Higgs_21[b][c][j];
              double f4 = Couplings_Gauge_Higgs_21[c][a][k];
              GaugePart[i][j][k] += 2 * f1 * f2 * f3 * f4;
            }
            double f1 = Couplings_Gauge_Higgs_22[a][b][i][j];
            double f2 = Couplings_Gauge_Higgs_21[b][a][k];
            double f3 = fbase(MassSquaredGauge[a], MassSquaredGauge[b]) -
                        C_CWcbGB + 0.5;
            GaugePart[i][j][k] += 3 * f1 * f2 * f3;
          }
        }

        QuarkPart[i][j][k] = 0;
        for (std::size_t a = 0; a < NQuarks; a++)
        {
          for (std::size_t b = 0; b < NQuarks; b++)
          {
            for (std::size_t c = 0; c < NQuarks; c++)
            {
              std::complex<double> f1 = fbaseTri(MassSquaredQuark[a],
                                                 MassSquaredQuark[b],
                                                 MassSquaredQuark[c]);
              std::complex<double> f2 = Couplings_Quark_Higgs_21[a][b][i];
              std::complex<double> f3 = Couplings_Quark_Higgs_21[b][c][j];
              std::complex<double> f4 = Couplings_Quark_Higgs_21[c][a][k];
              QuarkPart[i][j][k] += 2.0 * f1 * f2 * f3 * f4;
            }
            std::complex<double> f1 = Couplings_Quark_Higgs_22[a][b][i][j];
            std::complex<double> f2 = Couplings_Quark_Higgs_21[b][a][k];
            std::complex<double> f3 =
                fbase(MassSquaredQuark[a], MassSquaredQuark[b]) -
                C_CWcbFermion + 0.5;
            QuarkPart[i][j][k] += 3.0 * f1 * f2 * f3;
          }
        }
        LeptonPart[i][j][k] = 0;
        for (std::size_t a = 0; a < NLepton; a++)
        {
          for (std::size_t b = 0; b < NLepton; b++)
          {
            for (std::size_t c = 0; c < NLepton; c++)
            {
              std::complex<double> f1 = fbaseTri(MassSquaredLepton[a],
                                                 MassSquaredLepton[b],
                                                 MassSquaredLepton[c]);
              std::complex<double> f2 = Couplings_Lepton_Higgs_21[a][b][i];
              std::complex<double> f3 = Couplings_Lepton_Higgs_21[b][c][j];
              std::complex<double> f4 = Couplings_Lepton_Higgs_21[c][a][k];
              LeptonPart[i][j][k] += 2.0 * f1 * f2 * f3 * f4;
            }
            std::complex<double> f1 = Couplings_Lepton_Higgs_22[a][b][i][j];
            std::complex<double> f2 = Couplings_Lepton_Higgs_21[b][a][k];
            std::complex<double> f3 =
                fbase(MassSquaredLepton[a], MassSquaredLepton[b]) -
                C_CWcbFermion + 0.5;
            LeptonPart[i][j][k] += 3.0 * f1 * f2 * f3;
          }
        }
      }
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        HiggspartSym[i][j][k] = Higgspart[i][j][k] + Higgspart[i][k][j];
        HiggspartSym[i][j][k] += Higgspart[j][i][k] + Higgspart[j][k][i];
        HiggspartSym[i][j][k] += Higgspart[k][i][j] + Higgspart[k][j][i];
        HiggspartSym[i][j][k] *= 1.0 / 6.0;

        GaugePartSym[i][j][k] = GaugePart[i][j][k] + GaugePart[i][k][j];
        GaugePartSym[i][j][k] += GaugePart[j][i][k] + GaugePart[j][k][i];
        GaugePartSym[i][j][k] += GaugePart[k][i][j] + GaugePart[k][j][i];
        GaugePartSym[i][j][k] *= 1.0 / 6.0;

        QuarkPartSym[i][j][k] = QuarkPart[i][j][k] + QuarkPart[i][k][j];
        QuarkPartSym[i][j][k] += QuarkPart[j][i][k] + QuarkPart[j][k][i];
        QuarkPartSym[i][j][k] += QuarkPart[k][i][j] + QuarkPart[k][j][i];
        QuarkPartSym[i][j][k] *= 1.0 / 6.0;

        LeptonPartSym[i][j][k] = LeptonPart[i][j][k] + LeptonPart[i][k][j];
        LeptonPartSym[i][j][k] += LeptonPart[j][i][k] + LeptonPart[j][k][i];
        LeptonPartSym[i][j][k] += LeptonPart[k][i][j] + LeptonPart[k][j][i];
        LeptonPartSym[i][j][k] *= 1.0 / 6.0;
      }
    }
  }

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        restmp[i][j][k] = 0.5 * HiggspartSym[i][j][k];
        restmp[i][j][k] += 1.5 * GaugePartSym[i][j][k];
        restmp[i][j][k] += -1.0 * LeptonPartSym[i][j][k];
        restmp[i][j][k] += -3.0 * QuarkPartSym[i][j][k];
      }
    }
  }

  for (std::size_t l = 0; l < NHiggs; l++)
  {
    for (std::size_t m = 0; m < NHiggs; m++)
    {
      for (std::size_t n = 0; n < NHiggs; n++)
      {
        resGaugeBase[l][m][n] = 0;
        for (std::size_t i = 0; i < NHiggs; i++)
        {
          for (std::size_t j = 0; j < NHiggs; j++)
          {
            for (std::size_t k = 0; k < NHiggs; k++)
            {
              double RotFac = HiggsRotationMatrix[i][l] *
                              HiggsRotationMatrix[j][m] *
                              HiggsRotationMatrix[k][n];
              resGaugeBase[l][m][n] += RotFac * restmp[i][j][k].real();
            }
          }
        }
        resGaugeBase[l][m][n] *= epsilon;
        if (std::abs(resGaugeBase[l][m][n]) < NumZero)
          resGaugeBase[l][m][n] = 0;
      }
    }
  }

  for (std::size_t l = 0; l < NHiggs; l++)
  {
    for (std::size_t m = 0; m < NHiggs; m++)
    {
      for (std::size_t n = 0; n < NHiggs; n++)
      {
        res.push_back(resGaugeBase[l][m][n]);
      }
    }
  }

  return res;
}

std::vector<double> Class_Potential_Origin::WeinbergForthDerivative() const
{

  if (not CalcCouplingsDone)
  {
    std::string retmes = __func__;
    retmes += " tries to use Physical couplings but they are not initialised.";
    throw std::runtime_error(retmes);
  }

  const double NumZero = std::pow(10, -10);
  double epsilon       = 1.0 / (16.0 * M_PI * M_PI);

  std::vector<double> res;

  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      restmp;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      QuarkPart;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      LeptonPart;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      QuarkPartSym;
  std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>
      LeptonPartSym;
  std::vector<std::vector<std::vector<std::vector<double>>>> resGaugeBase;
  std::vector<std::vector<std::vector<std::vector<double>>>> HiggsPart;
  std::vector<std::vector<std::vector<std::vector<double>>>> GaugePart;
  std::vector<std::vector<std::vector<std::vector<double>>>> HiggsPartSym;
  std::vector<std::vector<std::vector<std::vector<double>>>> GaugePartSym;

  restmp.resize(NHiggs);
  QuarkPart.resize(NHiggs);
  LeptonPart.resize(NHiggs);
  QuarkPartSym.resize(NHiggs);
  LeptonPartSym.resize(NHiggs);
  resGaugeBase.resize(NHiggs);
  HiggsPart.resize(NHiggs);
  GaugePart.resize(NHiggs);
  HiggsPartSym.resize(NHiggs);
  GaugePartSym.resize(NHiggs);

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    restmp[i].resize(NHiggs);
    QuarkPart[i].resize(NHiggs);
    LeptonPart[i].resize(NHiggs);
    QuarkPartSym[i].resize(NHiggs);
    LeptonPartSym[i].resize(NHiggs);
    resGaugeBase[i].resize(NHiggs);
    HiggsPart[i].resize(NHiggs);
    GaugePart[i].resize(NHiggs);
    HiggsPartSym[i].resize(NHiggs);
    GaugePartSym[i].resize(NHiggs);
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      restmp[i][j].resize(NHiggs);
      QuarkPart[i][j].resize(NHiggs);
      LeptonPart[i][j].resize(NHiggs);
      QuarkPartSym[i][j].resize(NHiggs);
      LeptonPartSym[i][j].resize(NHiggs);
      resGaugeBase[i][j].resize(NHiggs);
      HiggsPart[i][j].resize(NHiggs);
      GaugePart[i][j].resize(NHiggs);
      HiggsPartSym[i][j].resize(NHiggs);
      GaugePartSym[i][j].resize(NHiggs);
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        restmp[i][j][k].resize(NHiggs);
        QuarkPart[i][j][k].resize(NHiggs);
        LeptonPart[i][j][k].resize(NHiggs);
        QuarkPartSym[i][j][k].resize(NHiggs);
        LeptonPartSym[i][j][k].resize(NHiggs);
        resGaugeBase[i][j][k].resize(NHiggs);
        HiggsPart[i][j][k].resize(NHiggs);
        GaugePart[i][j][k].resize(NHiggs);
        HiggsPartSym[i][j][k].resize(NHiggs);
        GaugePartSym[i][j][k].resize(NHiggs);
      }
    }
  }

  for (std::size_t i1 = 0; i1 < NHiggs; i1++)
  {
    for (std::size_t i2 = 0; i2 < NHiggs; i2++)
    {
      for (std::size_t i3 = 0; i3 < NHiggs; i3++)
      {
        for (std::size_t i4 = 0; i4 < NHiggs; i4++)
        {
          HiggsPart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NHiggs; a++)
          {
            for (std::size_t b = 0; b < NHiggs; b++)
            {
              double f1 = Couplings_Higgs_Quartic[a][b][i1][i4];
              double f2 = Couplings_Higgs_Quartic[b][a][i2][i3];
              double f3 = fbase(MassSquaredHiggs[a], MassSquaredHiggs[b]) -
                          C_CWcbHiggs + 0.5;

              HiggsPart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NHiggs; c++)
              {
                for (std::size_t d = 0; d < NHiggs; d++)
                {
                  double f11 = fbaseFour(MassSquaredHiggs[a],
                                         MassSquaredHiggs[b],
                                         MassSquaredHiggs[c],
                                         MassSquaredHiggs[d]);
                  double f12 = Couplings_Higgs_Triple[a][b][i4];
                  double f13 = Couplings_Higgs_Triple[b][c][i1];
                  double f14 = Couplings_Higgs_Triple[c][d][i2];
                  double f15 = Couplings_Higgs_Triple[d][a][i3];

                  HiggsPart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                double f21 = fbaseTri(MassSquaredHiggs[a],
                                      MassSquaredHiggs[b],
                                      MassSquaredHiggs[c]);
                double f22 = Couplings_Higgs_Triple[a][b][i4];
                double f23 = Couplings_Higgs_Quartic[b][c][i1][i2];
                double f24 = Couplings_Higgs_Triple[c][a][i3];

                HiggsPart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }

          GaugePart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NGauge; a++)
          {
            for (std::size_t b = 0; b < NGauge; b++)
            {
              double f1 = Couplings_Gauge_Higgs_22[a][b][i1][i4];
              double f2 = Couplings_Gauge_Higgs_22[b][a][i2][i3];
              double f3 = fbase(MassSquaredGauge[a], MassSquaredGauge[b]) -
                          C_CWcbGB + 0.5;

              GaugePart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NGauge; c++)
              {
                for (std::size_t d = 0; d < NGauge; d++)
                {
                  double f11 = fbaseFour(MassSquaredGauge[a],
                                         MassSquaredGauge[b],
                                         MassSquaredGauge[c],
                                         MassSquaredGauge[d]);
                  double f12 = Couplings_Gauge_Higgs_21[a][b][i4];
                  double f13 = Couplings_Gauge_Higgs_21[b][c][i1];
                  double f14 = Couplings_Gauge_Higgs_21[c][d][i2];
                  double f15 = Couplings_Gauge_Higgs_21[d][a][i3];

                  GaugePart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                double f21 = fbaseTri(MassSquaredGauge[a],
                                      MassSquaredGauge[b],
                                      MassSquaredGauge[c]);
                double f22 = Couplings_Gauge_Higgs_21[a][b][i4];
                double f23 = Couplings_Gauge_Higgs_22[b][c][i1][i2];
                double f24 = Couplings_Gauge_Higgs_21[c][a][i3];

                GaugePart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }

          QuarkPart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NQuarks; a++)
          {
            for (std::size_t b = 0; b < NQuarks; b++)
            {
              std::complex<double> f1 = Couplings_Quark_Higgs_22[a][b][i1][i4];
              std::complex<double> f2 = Couplings_Quark_Higgs_22[b][a][i2][i3];
              std::complex<double> f3 =
                  fbase(MassSquaredQuark[a], MassSquaredQuark[b]) -
                  C_CWcbFermion + 0.5;

              QuarkPart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NQuarks; c++)
              {
                for (std::size_t d = 0; d < NQuarks; d++)
                {
                  std::complex<double> f11 = fbaseFour(MassSquaredQuark[a],
                                                       MassSquaredQuark[b],
                                                       MassSquaredQuark[c],
                                                       MassSquaredQuark[d]);
                  std::complex<double> f12 = Couplings_Quark_Higgs_21[a][b][i4];
                  std::complex<double> f13 = Couplings_Quark_Higgs_21[b][c][i1];
                  std::complex<double> f14 = Couplings_Quark_Higgs_21[c][d][i2];
                  std::complex<double> f15 = Couplings_Quark_Higgs_21[d][a][i3];

                  QuarkPart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                std::complex<double> f21 = fbaseTri(MassSquaredQuark[a],
                                                    MassSquaredQuark[b],
                                                    MassSquaredQuark[c]);
                std::complex<double> f22 = Couplings_Quark_Higgs_21[a][b][i4];
                std::complex<double> f23 =
                    Couplings_Quark_Higgs_22[b][c][i1][i2];
                std::complex<double> f24 = Couplings_Quark_Higgs_21[c][a][i3];

                QuarkPart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }

          LeptonPart[i1][i2][i3][i4] = 0;
          for (std::size_t a = 0; a < NLepton; a++)
          {
            for (std::size_t b = 0; b < NLepton; b++)
            {
              std::complex<double> f1 = Couplings_Lepton_Higgs_22[a][b][i1][i4];
              std::complex<double> f2 = Couplings_Lepton_Higgs_22[b][a][i2][i3];
              std::complex<double> f3 =
                  fbase(MassSquaredLepton[a], MassSquaredLepton[b]) -
                  C_CWcbFermion + 0.5;

              LeptonPart[i1][i2][i3][i4] += f1 * f2 * f3;

              for (std::size_t c = 0; c < NLepton; c++)
              {
                for (std::size_t d = 0; d < NLepton; d++)
                {
                  std::complex<double> f11 = fbaseFour(MassSquaredLepton[a],
                                                       MassSquaredLepton[b],
                                                       MassSquaredLepton[c],
                                                       MassSquaredLepton[d]);
                  std::complex<double> f12 =
                      Couplings_Lepton_Higgs_21[a][b][i4];
                  std::complex<double> f13 =
                      Couplings_Lepton_Higgs_21[b][c][i1];
                  std::complex<double> f14 =
                      Couplings_Lepton_Higgs_21[c][d][i2];
                  std::complex<double> f15 =
                      Couplings_Lepton_Higgs_21[d][a][i3];

                  LeptonPart[i1][i2][i3][i4] +=
                      2.0 * f11 * f12 * f13 * f14 * f15;
                }

                std::complex<double> f21 = fbaseTri(MassSquaredLepton[a],
                                                    MassSquaredLepton[b],
                                                    MassSquaredLepton[c]);
                std::complex<double> f22 = Couplings_Lepton_Higgs_21[a][b][i4];
                std::complex<double> f23 =
                    Couplings_Lepton_Higgs_22[b][c][i1][i2];
                std::complex<double> f24 = Couplings_Lepton_Higgs_21[c][a][i3];

                LeptonPart[i1][i2][i3][i4] += 4.0 * f21 * f22 * f23 * f24;
              }
            }
          }
        }
      }
    }
  }

  for (std::size_t i1 = 0; i1 < NHiggs; i1++)
  {
    for (std::size_t i2 = 0; i2 < NHiggs; i2++)
    {
      for (std::size_t i3 = 0; i3 < NHiggs; i3++)
      {
        for (std::size_t i4 = 0; i4 < NHiggs; i4++)
        {
          HiggsPartSym[i1][i2][i3][i4] =
              HiggsPart[i1][i2][i3][i4] + HiggsPart[i1][i2][i4][i3] +
              HiggsPart[i1][i3][i2][i4] + HiggsPart[i1][i3][i4][i2] +
              HiggsPart[i1][i4][i2][i3] + HiggsPart[i1][i4][i3][i2];
          HiggsPartSym[i1][i2][i3][i4] +=
              HiggsPart[i2][i1][i3][i4] + HiggsPart[i2][i1][i4][i3] +
              HiggsPart[i2][i3][i1][i4] + HiggsPart[i2][i3][i4][i1] +
              HiggsPart[i2][i4][i1][i3] + HiggsPart[i2][i4][i3][i1];
          HiggsPartSym[i1][i2][i3][i4] +=
              HiggsPart[i3][i1][i2][i4] + HiggsPart[i3][i1][i4][i2] +
              HiggsPart[i3][i2][i1][i4] + HiggsPart[i3][i2][i4][i1] +
              HiggsPart[i3][i4][i1][i2] + HiggsPart[i3][i4][i2][i1];
          HiggsPartSym[i1][i2][i3][i4] +=
              HiggsPart[i4][i1][i2][i3] + HiggsPart[i4][i1][i3][i2] +
              HiggsPart[i4][i2][i1][i3] + HiggsPart[i4][i2][i3][i1] +
              HiggsPart[i4][i3][i1][i2] + HiggsPart[i4][i3][i2][i1];
          HiggsPartSym[i1][i2][i3][i4] *= 1.0 / 24.0;

          GaugePartSym[i1][i2][i3][i4] =
              GaugePart[i1][i2][i3][i4] + GaugePart[i1][i2][i4][i3] +
              GaugePart[i1][i3][i2][i4] + GaugePart[i1][i3][i4][i2] +
              GaugePart[i1][i4][i2][i3] + GaugePart[i1][i4][i3][i2];
          GaugePartSym[i1][i2][i3][i4] +=
              GaugePart[i2][i1][i3][i4] + GaugePart[i2][i1][i4][i3] +
              GaugePart[i2][i3][i1][i4] + GaugePart[i2][i3][i4][i1] +
              GaugePart[i2][i4][i1][i3] + GaugePart[i2][i4][i3][i1];
          GaugePartSym[i1][i2][i3][i4] +=
              GaugePart[i3][i1][i2][i4] + GaugePart[i3][i1][i4][i2] +
              GaugePart[i3][i2][i1][i4] + GaugePart[i3][i2][i4][i1] +
              GaugePart[i3][i4][i1][i2] + GaugePart[i3][i4][i2][i1];
          GaugePartSym[i1][i2][i3][i4] +=
              GaugePart[i4][i1][i2][i3] + GaugePart[i4][i1][i3][i2] +
              GaugePart[i4][i2][i1][i3] + GaugePart[i4][i2][i3][i1] +
              GaugePart[i4][i3][i1][i2] + GaugePart[i4][i3][i2][i1];
          GaugePartSym[i1][i2][i3][i4] *= 1.0 / 24.0;

          QuarkPartSym[i1][i2][i3][i4] =
              QuarkPart[i1][i2][i3][i4] + QuarkPart[i1][i2][i4][i3] +
              QuarkPart[i1][i3][i2][i4] + QuarkPart[i1][i3][i4][i2] +
              QuarkPart[i1][i4][i2][i3] + QuarkPart[i1][i4][i3][i2];
          QuarkPartSym[i1][i2][i3][i4] +=
              QuarkPart[i2][i1][i3][i4] + QuarkPart[i2][i1][i4][i3] +
              QuarkPart[i2][i3][i1][i4] + QuarkPart[i2][i3][i4][i1] +
              QuarkPart[i2][i4][i1][i3] + QuarkPart[i2][i4][i3][i1];
          QuarkPartSym[i1][i2][i3][i4] +=
              QuarkPart[i3][i1][i2][i4] + QuarkPart[i3][i1][i4][i2] +
              QuarkPart[i3][i2][i1][i4] + QuarkPart[i3][i2][i4][i1] +
              QuarkPart[i3][i4][i1][i2] + QuarkPart[i3][i4][i2][i1];
          QuarkPartSym[i1][i2][i3][i4] +=
              QuarkPart[i4][i1][i2][i3] + QuarkPart[i4][i1][i3][i2] +
              QuarkPart[i4][i2][i1][i3] + QuarkPart[i4][i2][i3][i1] +
              QuarkPart[i4][i3][i1][i2] + QuarkPart[i4][i3][i2][i1];
          QuarkPartSym[i1][i2][i3][i4] *= 1.0 / 24.0;

          LeptonPartSym[i1][i2][i3][i4] =
              LeptonPart[i1][i2][i3][i4] + LeptonPart[i1][i2][i4][i3] +
              LeptonPart[i1][i3][i2][i4] + LeptonPart[i1][i3][i4][i2] +
              LeptonPart[i1][i4][i2][i3] + LeptonPart[i1][i4][i3][i2];
          LeptonPartSym[i1][i2][i3][i4] +=
              LeptonPart[i2][i1][i3][i4] + LeptonPart[i2][i1][i4][i3] +
              LeptonPart[i2][i3][i1][i4] + LeptonPart[i2][i3][i4][i1] +
              LeptonPart[i2][i4][i1][i3] + LeptonPart[i2][i4][i3][i1];
          LeptonPartSym[i1][i2][i3][i4] +=
              LeptonPart[i3][i1][i2][i4] + LeptonPart[i3][i1][i4][i2] +
              LeptonPart[i3][i2][i1][i4] + LeptonPart[i3][i2][i4][i1] +
              LeptonPart[i3][i4][i1][i2] + LeptonPart[i3][i4][i2][i1];
          LeptonPartSym[i1][i2][i3][i4] +=
              LeptonPart[i4][i1][i2][i3] + LeptonPart[i4][i1][i3][i2] +
              LeptonPart[i4][i2][i1][i3] + LeptonPart[i4][i2][i3][i1] +
              LeptonPart[i4][i3][i1][i2] + LeptonPart[i4][i3][i2][i1];
          LeptonPartSym[i1][i2][i3][i4] *= 1.0 / 24.0;
        }
      }
    }
  }

  for (std::size_t i1 = 0; i1 < NHiggs; i1++)
  {
    for (std::size_t i2 = 0; i2 < NHiggs; i2++)
    {
      for (std::size_t i3 = 0; i3 < NHiggs; i3++)
      {
        for (std::size_t i4 = 0; i4 < NHiggs; i4++)
        {
          restmp[i1][i2][i3][i4] = 3.0 * 0.5 * HiggsPartSym[i1][i2][i3][i4];
          restmp[i1][i2][i3][i4] += 3.0 * 1.5 * GaugePartSym[i1][i2][i3][i4];
          restmp[i1][i2][i3][i4] +=
              3.0 * (-1.0) * LeptonPartSym[i1][i2][i3][i4];
          restmp[i1][i2][i3][i4] += 3.0 * (-3.0) * QuarkPartSym[i1][i2][i3][i4];
        }
      }
    }
  }

  for (std::size_t j1 = 0; j1 < NHiggs; j1++)
  {
    for (std::size_t j2 = 0; j2 < NHiggs; j2++)
    {
      for (std::size_t j3 = 0; j3 < NHiggs; j3++)
      {
        for (std::size_t j4 = 0; j4 < NHiggs; j4++)
        {
          resGaugeBase[j1][j2][j3][j4] = 0;

          for (std::size_t i1 = 0; i1 < NHiggs; i1++)
          {
            for (std::size_t i2 = 0; i2 < NHiggs; i2++)
            {
              for (std::size_t i3 = 0; i3 < NHiggs; i3++)
              {
                for (std::size_t i4 = 0; i4 < NHiggs; i4++)
                {
                  double RotFac = HiggsRotationMatrix[i1][j1] *
                                  HiggsRotationMatrix[i2][j2] *
                                  HiggsRotationMatrix[i3][j3] *
                                  HiggsRotationMatrix[i4][j4];
                  resGaugeBase[j1][j2][j3][j4] +=
                      RotFac * restmp[i1][i2][i3][i4].real();
                }
              }
            }
          }

          resGaugeBase[j1][j2][j3][j4] *= epsilon;
          if (std::abs(resGaugeBase[j1][j2][j3][j4]) < NumZero)
            resGaugeBase[j1][j2][j3][j4] = 0;
        }
      }
    }
  }

  for (std::size_t j1 = 0; j1 < NHiggs; j1++)
  {
    for (std::size_t j2 = 0; j2 < NHiggs; j2++)
    {
      for (std::size_t j3 = 0; j3 < NHiggs; j3++)
      {
        for (std::size_t j4 = 0; j4 < NHiggs; j4++)
        {
          res.push_back(resGaugeBase[j4][j3][j2][j1]);
        }
      }
    }
  }

  return res;
}

MatrixXd Class_Potential_Origin::HiggsMassMatrix(const std::vector<double> &v,
                                                 double Temp,
                                                 int diff) const
{
  MatrixXd res(NHiggs, NHiggs);
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    res          = HiggsMassMatrix(Transformedv, Temp, diff);
    return res;
  }
  if (!SetCurvatureDone)
  {
    //        SetCurvatureArrays();
    throw std::runtime_error(
        "SetCurvatureDone is not set. The Model is not initiliased correctly");
  }

  if (diff == 0)
  {
    static const bool SkipZeroFields = []
    {
      const char *env = std::getenv("BSMPT_SKIP_ZERO_HIGGS_FIELDS");
      return env != nullptr && env[0] == '1';
    }();
    static const bool SkipZeroR2HDMTensorTerms = []
    {
      const char *env = std::getenv("BSMPT_SKIP_ZERO_R2HDM_HIGGS_TERMS");
      return env != nullptr && env[0] == '1';
    }();
    const bool useSparseTerms = SkipZeroR2HDMTensorTerms &&
                                Model == ModelID::ModelIDs::R2HDM &&
                                NHiggs == 8;
    static const bool UseR2HDMIndexCache = []
    {
      const char *env = std::getenv("BSMPT_USE_R2HDM_HIGGS_INDEX_CACHE");
      return env != nullptr && env[0] == '1';
    }();
    const bool useIndexCache = UseR2HDMIndexCache &&
                               Model == ModelID::ModelIDs::R2HDM &&
                               NHiggs == 8 && HiggsTensorIndexCacheReady &&
                               HiggsTensorIndexCache.size() == NHiggs * NHiggs;
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = i; j < NHiggs; j++)
      {
        res(i, j) = Curvature_Higgs_L2[i][j];
        if (useIndexCache)
        {
          const auto &terms = HiggsTensorIndexCache[i * NHiggs + j];
          for (const auto &term : terms)
          {
            const std::size_t k = term.k;
            if (SkipZeroFields && v[k] == 0.0) continue;
            if (term.l3Nonzero)
              res(i, j) += Curvature_Higgs_L3[i][j][k] * v[k];
            for (const std::size_t l : term.l4Nonzero)
            {
              if (SkipZeroFields && v[l] == 0.0) continue;
              res(i, j) +=
                  0.5 * Curvature_Higgs_L4[i][j][k][l] * v[k] * v[l];
            }
          }
        }
        else
        {
          for (std::size_t k = 0; k < NHiggs; k++)
          {
            if (SkipZeroFields && v[k] == 0.0) continue;
            if (!useSparseTerms || Curvature_Higgs_L3[i][j][k] != 0.0)
              res(i, j) += Curvature_Higgs_L3[i][j][k] * v[k];
            for (std::size_t l = 0; l < NHiggs; l++)
            {
              if ((SkipZeroFields && v[l] == 0.0) ||
                  (useSparseTerms &&
                   Curvature_Higgs_L4[i][j][k][l] == 0.0))
                continue;
              res(i, j) +=
                  0.5 * Curvature_Higgs_L4[i][j][k][l] * v[k] * v[l];
            }
          }
        }

        if (Temp != 0)
        {
          res(i, j) += DebyeHiggs[i][j] * std::pow(Temp, 2);
        }
      }
    }
    for (std::size_t i{1}; i < NHiggs; ++i)
    {
      for (std::size_t j{0}; j < i; ++j)
      {
        res(i, j) = res(j, i);
      }
    }
  }
  else if (static_cast<size_t>(diff) <= NHiggs and diff > 0)
  {
    static const bool SkipZeroFields = []
    {
      const char *env = std::getenv("BSMPT_SKIP_ZERO_HIGGS_FIELDS");
      return env != nullptr && env[0] == '1';
    }();
    static const bool SkipZeroR2HDMTensorTerms = []
    {
      const char *env = std::getenv("BSMPT_SKIP_ZERO_R2HDM_HIGGS_TERMS");
      return env != nullptr && env[0] == '1';
    }();
    const bool useSparseTerms = SkipZeroR2HDMTensorTerms &&
                                Model == ModelID::ModelIDs::R2HDM &&
                                NHiggs == 8;
    static const bool UseUpperOnlyR2HDMDerivative = []
    {
      const char *env =
          std::getenv("BSMPT_SAFE_R2HDM_HIGGS_DERIVATIVE_UPPER_ONLY");
      return env != nullptr && env[0] == '1';
    }();
    const bool useUpperOnly = UseUpperOnlyR2HDMDerivative &&
                              Model == ModelID::ModelIDs::R2HDM &&
                              NHiggs == 8;
    std::size_t x0 = diff - 1;
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      const std::size_t jBegin = useUpperOnly ? i : 0;
      for (std::size_t j = jBegin; j < NHiggs; j++)
      {
        res(i, j) = Curvature_Higgs_L3[i][j][x0];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          if ((SkipZeroFields && v[k] == 0.0) ||
              (useSparseTerms &&
               Curvature_Higgs_L4[i][j][x0][k] == 0.0))
            continue;
          res(i, j) += Curvature_Higgs_L4[i][j][x0][k] * v[k];
        }
      }
    }
    if (useUpperOnly)
    {
      for (std::size_t i = 1; i < NHiggs; ++i)
      {
        for (std::size_t j = 0; j < i; ++j)
          res(i, j) = res(j, i);
      }
    }
  }
  else if (diff == -1)
  {
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        res(i, j) = 2 * DebyeHiggs[i][j] * Temp;
      }
    }
  }
  return res;
}

std::vector<double>
Class_Potential_Origin::HiggsMassesSquared(const std::vector<double> &v,
                                           const double &Temp,
                                           const int &diff) const
{
  ProfileExactModelRepeat(CalcGWProfiler::ExactRepeatMetric::HiggsMasses,
                          this, v, Temp, diff, 0);
  std::vector<double> res;
  res.reserve((diff > 0 || diff == -1) ? 2 * NHiggs : NHiggs);

  auto MassMatrix = HiggsMassMatrix(v, Temp);

  double ZeroMass = std::pow(10, -5);
  static const bool UseRealEigenDerivative = []
  {
    const char *env = std::getenv("BSMPT_USE_REAL_EIGEN_DERIVATIVE");
    return env != nullptr && env[0] == '1';
  }();

  if (diff == 0 and res.empty())
  {
    SelfAdjointEigenSolver<MatrixXd> es(MassMatrix, EigenvaluesOnly);
    const auto EV = es.eigenvalues();
    for (std::size_t i{0}; i < NHiggs; ++i)
    {
      if (std::abs(EV[i]) < ZeroMass)
      {
        res.push_back(0);
      }
      else
      {
        res.push_back(EV[i]);
      }
    }
  }
  else if (diff == 0 and res.size() == NHiggs)
  {
    SelfAdjointEigenSolver<MatrixXd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      double tmp = es.eigenvalues()[i];
      if (std::abs(tmp) < ZeroMass) tmp = 0;
      res[i] = tmp;
    }
  }
  else if (diff == 0 and res.size() != 0 and res.size() != NHiggs)
  {
    Logger::Write(LoggingLevel::Debug,
                  std::string("Something went wrong in ") + __func__ + ".\n" +
                      __func__ + "Is calculating the mass for " +
                      std::to_string(NHiggs) +
                      "fields but the resolution vector has a size of " +
                      std::to_string(res.size()) + ". This should be zero or " +
                      std::to_string(NHiggs));
  }
  else if (static_cast<std::size_t>(diff) <= NHiggs and diff > 0)
  {
    MatrixXd Diff(NHiggs, NHiggs);
    std::size_t x0 = diff - 1;
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        Diff(i, j) = Curvature_Higgs_L3[i][j][x0];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          Diff(i, j) += Curvature_Higgs_L4[i][j][x0][k] * v[k];
        }
      }
    }
    if (UseRealEigenDerivative)
    {
      res = FirstDerivativeOfRealEigenvalues(MassMatrix, Diff);
    }
    if (res.empty())
    {
      MatrixXcd MassCast(NHiggs, NHiggs);
      MassCast = MassMatrix;
      MatrixXcd DiffCast(NHiggs, NHiggs);
      DiffCast = Diff;
      res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
    }
  }
  else if (diff == -1)
  {
    MatrixXd Diff(NHiggs, NHiggs);
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        Diff(i, j) = 2 * DebyeHiggs[i][j] * Temp;
      }
    }

    MatrixXcd MassCast(NHiggs, NHiggs);
    MassCast = MassMatrix;
    MatrixXcd DiffCast(NHiggs, NHiggs);
    DiffCast = Diff;
    res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
  }

  return res;
}

std::vector<double>
Class_Potential_Origin::GaugeMassesSquared(const std::vector<double> &v,
                                           const double &Temp,
                                           const int &diff) const
{
  std::vector<double> res;
  res.reserve((diff > 0 || diff == -1) ? 2 * NGauge : NGauge);
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    res          = GaugeMassesSquared(Transformedv, Temp, diff);
    return res;
  }
  if (!SetCurvatureDone)
  {
    //        SetCurvatureArrays();
    std::string retmes = __func__;
    retmes += "was called while the model was not initialised correctly.\n";
    throw std::runtime_error(retmes);
  }
  MatrixXd MassMatrix(NGauge, NGauge);
  double ZeroMass = std::pow(10, -5);
  static const bool UseRealEigenDerivative = []
  {
    const char *env = std::getenv("BSMPT_USE_REAL_EIGEN_DERIVATIVE");
    return env != nullptr && env[0] == '1';
  }();
  static const bool UseUpperOnly = []
  {
    const char *env = std::getenv("BSMPT_SAFE_GAUGE_UPPER_ONLY");
    return env != nullptr && env[0] == '1';
  }();
  static const bool SkipZeroGaugeTensorTerms = []
  {
    const char *env = std::getenv("BSMPT_SKIP_ZERO_R2HDM_GAUGE_TERMS");
    return env != nullptr && env[0] == '1';
  }();
  const bool useSparseGaugeTerms = SkipZeroGaugeTensorTerms &&
                                   Model == ModelID::ModelIDs::R2HDM &&
                                   NHiggs == 8 && NGauge == 4;
  for (std::size_t a = 0; a < NGauge; a++)
  {
    const std::size_t bBegin = UseUpperOnly ? a : 0;
    for (std::size_t b = bBegin; b < NGauge; b++)
    {
      MassMatrix(a, b) = 0;
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          if (useSparseGaugeTerms &&
              Curvature_Gauge_G2H2[a][b][i][j] == 0.0)
            continue;
          MassMatrix(a, b) +=
              0.5 * Curvature_Gauge_G2H2[a][b][i][j] * v.at(i) * v.at(j);
        }
      }

      if (Temp != 0)
      {
        MassMatrix(a, b) += DebyeGauge[a][b] * std::pow(Temp, 2);
      }
    }
  }

  for (std::size_t a{1}; a < NGauge; ++a)
  {
    for (std::size_t b{0}; b < a; ++b)
    {
      MassMatrix(a, b) = MassMatrix(b, a);
    }
  }

  // R2HDM has an exact arrowhead gauge-mass matrix for arbitrary real CB/CP
  // field configurations.  Keep this path opt-in and validate the complete
  // matrix before using the closed form, so that a model change or a
  // convention mismatch falls back to the general Eigen implementation below.
  static const bool UseAnalyticR2HDMGauge = []
  {
    const char *env = std::getenv("BSMPT_USE_ANALYTIC_GAUGE_MASSES");
    return env != nullptr && env[0] == '1';
  }();
  if (UseAnalyticR2HDMGauge && diff == 0 &&
      Model == ModelID::ModelIDs::R2HDM && NHiggs == 8 && NGauge == 4)
  {
    const double g  = SMConstants.C_g;
    const double gp = SMConstants.C_gs;
    double fieldNormSquared = 0;
    for (std::size_t i = 0; i < 8; ++i)
      fieldNormSquared += v[i] * v[i];

    const double q = v[0] * v[4] + v[1] * v[5] + v[2] * v[6] + v[3] * v[7];
    const double p = v[0] * v[5] - v[1] * v[4] + v[2] * v[7] - v[3] * v[6];
    const double d = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]
                     - v[4] * v[4] - v[5] * v[5] - v[6] * v[6] - v[7] * v[7];

    const double halfGGp = 0.5 * g * gp;
    MatrixXd analyticMatrix = MatrixXd::Zero(4, 4);
    const double A = 0.25 * g * g * fieldNormSquared;
    const double B = 0.25 * gp * gp * fieldNormSquared;
    analyticMatrix(0, 0) = A;
    analyticMatrix(1, 1) = A;
    analyticMatrix(2, 2) = A;
    analyticMatrix(3, 3) = B;
    analyticMatrix(0, 3) = analyticMatrix(3, 0) = halfGGp * q;
    analyticMatrix(1, 3) = analyticMatrix(3, 1) = halfGGp * p;
    analyticMatrix(2, 3) = analyticMatrix(3, 2) = halfGGp * d;
    if (Temp != 0)
    {
      const double tempSquared = Temp * Temp;
      for (std::size_t a = 0; a < 4; ++a)
        analyticMatrix(a, a) += DebyeGauge[a][a] * tempSquared;
    }

    double matrixScale = 1;
    double matrixError = 0;
    bool finiteMatrices = true;
    for (std::size_t a = 0; a < 4; ++a)
    {
      for (std::size_t b = 0; b < 4; ++b)
      {
        if (!std::isfinite(MassMatrix(a, b)) ||
            !std::isfinite(analyticMatrix(a, b)))
        {
          finiteMatrices = false;
        }
        matrixScale = std::max(matrixScale, std::abs(MassMatrix(a, b)));
        matrixScale = std::max(matrixScale, std::abs(analyticMatrix(a, b)));
        matrixError = std::max(matrixError,
                               std::abs(MassMatrix(a, b) - analyticMatrix(a, b)));
      }
    }

    // The contraction loop and the invariant formula have different floating
    // point addition order.  This tolerance is only a structural validation;
    // it is many orders below the physical scale of the matrix.
    const double matrixTolerance = 1e-12 * matrixScale;
    if (finiteMatrices && std::isfinite(matrixError) &&
        matrixError <= matrixTolerance)
    {
      const double tNorm = std::hypot(halfGGp * q,
                                      std::hypot(halfGGp * p, halfGGp * d));
      const double AT = analyticMatrix(0, 0);
      const double BT = analyticMatrix(3, 3);
      const double delta = std::hypot(AT - BT, 2 * tNorm);
      const double lambdaPlusCorrect = 0.5 * (AT + BT + delta);
      double lambdaMinus = 0;
      const double formulaScale =
          std::max({1.0, std::abs(AT), std::abs(BT), tNorm});
      if (std::abs(lambdaPlusCorrect) > 1e-14 * formulaScale)
        lambdaMinus = (AT * BT - tNorm * tNorm) / lambdaPlusCorrect;
      else
        lambdaMinus = 0.5 * (AT + BT - delta);

      std::vector<double> analyticEigenvalues{
          lambdaMinus, AT, AT, lambdaPlusCorrect};
      std::sort(analyticEigenvalues.begin(), analyticEigenvalues.end());
      res.reserve(4);
      for (double value : analyticEigenvalues)
      {
        if (std::abs(value) < ZeroMass) value = 0;
        res.push_back(value);
      }
      return res;
    }
  }

  if (diff == 0)
  {

    SelfAdjointEigenSolver<MatrixXd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NGauge; i++)
    {
      double tmp = es.eigenvalues()[i];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t i = diff - 1;
    MatrixXd Diff(NGauge, NGauge);
    Diff = MatrixXd::Zero(NGauge, NGauge);
    for (std::size_t a = 0; a < NGauge; a++)
    {
      for (std::size_t b = 0; b < NGauge; b++)
      {
        for (std::size_t j = 0; j < NHiggs; j++)
          Diff(a, b) += Curvature_Gauge_G2H2[a][b][i][j] * v[j];
      }
    }
    if (UseRealEigenDerivative)
    {
      res = FirstDerivativeOfRealEigenvalues(MassMatrix, Diff);
    }
    if (res.empty())
    {
      MatrixXcd MassCast(NGauge, NGauge);
      MassCast = MassMatrix;
      MatrixXcd DiffCast(NGauge, NGauge);
      DiffCast = Diff;
      res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
    }
  }
  else if (diff == -1)
  {
    MatrixXd Diff(NGauge, NGauge);
    for (std::size_t i = 0; i < NGauge; i++)
    {
      for (std::size_t j = 0; j < NGauge; j++)
      {
        Diff(i, j) = 2 * DebyeGauge[i][j] * Temp;
      }
    }

    if (UseRealEigenDerivative)
    {
      res = FirstDerivativeOfRealEigenvalues(MassMatrix, Diff);
    }
    if (res.empty())
    {
      MatrixXcd MassCast(NGauge, NGauge);
      MassCast = MassMatrix;
      MatrixXcd DiffCast(NGauge, NGauge);
      DiffCast = Diff;
      res      = FirstDerivativeOfEigenvalues(MassCast, DiffCast);
    }
  }

  return res;
}

std::vector<double>
Class_Potential_Origin::QuarkMassesSquared(const std::vector<double> &v,
                                           const int &diff) const
{
  std::vector<double> res;
  QuarkMassesSquared(v, diff, res);
  return res;
}

void Class_Potential_Origin::QuarkMassesSquared(
    const std::vector<double> &v,
    const int &diff,
    std::vector<double> &res) const
{
  ProfileExactModelRepeat(CalcGWProfiler::ExactRepeatMetric::QuarkMasses,
                          this, v, 0.0, diff, 0);
  res.clear();
  res.reserve(diff > 0 ? 2 * NQuarks : NQuarks);
  double ZeroMass = std::pow(10, -10);
  static const bool UseR2HDMQuarkFixed12Diff0 = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_QUARK_FIXED12_DIFF0");
    return env != nullptr && env[0] == '1';
  }();
  static const bool UseR2HDMQuarkIndexCache = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_QUARK_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool useIndexCache =
      UseR2HDMQuarkIndexCache && Model == ModelID::ModelIDs::R2HDM &&
      NQuarks == 12 && QuarkTensorIndexCacheReady &&
      QuarkDerivativeIndexCache.size() == NHiggs * NQuarks * NQuarks;

  // The fixed-size path below used to call QuarkMassMatrix first and then
  // copy its dynamic MatrixXcd into a 12x12 object.  For the R2HDM input
  // shape (v.size()==NHiggs), assemble the same entries directly into the
  // fixed object.  This changes no arithmetic in the Yukawa sum and leaves
  // the original dynamic implementation as the fallback for all other
  // models/input layouts.  The option remains opt-in because Eigen's fixed
  // and dynamic kernels can differ in their last bit.
  const bool useDirectFixed12 =
      UseR2HDMQuarkFixed12Diff0 && diff <= 0 &&
      Model == ModelID::ModelIDs::R2HDM && NQuarks == 12 &&
      v.size() == NHiggs && SetCurvatureDone;
  using FixedQuarkMatrix = Eigen::Matrix<std::complex<double>, 12, 12>;
  FixedQuarkMatrix MIJFixed;
  if (useDirectFixed12)
  {
    MIJFixed.setZero();
    for (std::size_t i = 0; i < NQuarks; ++i)
    {
      for (std::size_t j = 0; j < NQuarks; ++j)
      {
        MIJFixed(i, j) = Curvature_Quark_F2[i][j];
        if (useIndexCache)
        {
          const auto &fields =
              QuarkMassFieldIndexCache[i * NQuarks + j];
          for (const std::size_t k : fields)
            MIJFixed(i, j) += Curvature_Quark_F2H1[i][j][k] * v[k];
        }
        else
        {
          for (std::size_t k = 0; k < NHiggs; ++k)
            MIJFixed(i, j) += Curvature_Quark_F2H1[i][j][k] * v[k];
        }
      }
    }
  }
  // Do not allocate dynamic fallback matrices on the fixed-size R2HDM path.
  // Eigen assignments below resize these objects when the fallback is used.
  MatrixXcd MassMatrix, MIJ;
  if (!useDirectFixed12)
    MIJ = QuarkMassMatrix(v);

  static const bool UseR2HDMBlock6 = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_QUARK_BLOCK6");
    return env != nullptr && env[0] == '1';
  }();
  if (UseR2HDMBlock6 && diff <= 0 &&
      Model == ModelID::ModelIDs::R2HDM && NQuarks == 12)
  {
    bool exact_chiral_structure = true;
    for (std::size_t i = 0; i < 6 && exact_chiral_structure; ++i)
    {
      for (std::size_t j = 0; j < 6; ++j)
      {
        if (MIJ(i, j) != std::complex<double>(0, 0) ||
            MIJ(i + 6, j + 6) != std::complex<double>(0, 0) ||
            MIJ(i, j + 6) != MIJ(j + 6, i))
        {
          exact_chiral_structure = false;
          break;
        }
      }
    }
    if (exact_chiral_structure)
    {
      const MatrixXcd K = MIJ.block(0, 6, 6, 6);
      const MatrixXcd reduced_mass_matrix = K.adjoint() * K;
      SelfAdjointEigenSolver<MatrixXcd> es(reduced_mass_matrix,
                                           EigenvaluesOnly);
      if (es.info() == Success)
      {
        for (std::size_t i = 0; i < 6; ++i)
        {
          double value = es.eigenvalues().real()[i];
          if (std::abs(value) < ZeroMass) value = 0;
          res.push_back(value);
          res.push_back(value);
        }
        return;
      }
    }
  }

  // For the standard R2HDM quark basis the zero-temperature mass matrix is
  // always 12x12.  Keep the same complete M^* M construction as the generic
  // path, but let Eigen see the compile-time dimensions.  This is opt-in
  // because fixed-size and dynamic Eigen expressions can differ in the last
  // bit, and the dynamic path remains the exact reference implementation.
  if (UseR2HDMQuarkFixed12Diff0 && diff <= 0 &&
      Model == ModelID::ModelIDs::R2HDM && NQuarks == 12)
  {
    if (!useDirectFixed12)
      MIJFixed = MIJ;
    const FixedQuarkMatrix MassMatrixFixed =
        MIJFixed.conjugate() * MIJFixed;
    SelfAdjointEigenSolver<FixedQuarkMatrix> es(MassMatrixFixed,
                                                EigenvaluesOnly);
    for (std::size_t i = 0; i < 12; ++i)
    {
      double tmp = es.eigenvalues().real()[i];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
    return;
  }

  MassMatrix = MIJ.conjugate() * MIJ;

  if (diff <= 0) // No temperature dependent part here
  {
    SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NQuarks; i++)
    {
      double tmp = es.eigenvalues().real()[i];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
  }
  else if (static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t m = diff - 1;
    MatrixXcd Diff(NQuarks, NQuarks);
    Diff = MatrixXcd::Zero(NQuarks, NQuarks);
    for (std::size_t a = 0; a < NQuarks; a++)
    {
      for (std::size_t b = 0; b < NQuarks; b++)
      {
        if (useIndexCache)
        {
          const auto &indices =
              QuarkDerivativeIndexCache[(m * NQuarks + a) * NQuarks + b];
          for (const auto &term : indices)
          {
            const std::size_t i = term.i;
            if (term.firstNonzero)
              Diff(a, b) +=
                  std::conj(Curvature_Quark_F2H1[a][i][m]) * MIJ(i, b);
            if (term.secondNonzero)
              Diff(a, b) +=
                  std::conj(MIJ(a, i)) * Curvature_Quark_F2H1[i][b][m];
          }
        }
        else
        {
          for (std::size_t i = 0; i < NQuarks; i++)
          {
            Diff(a, b) +=
                std::conj(Curvature_Quark_F2H1[a][i][m]) * MIJ(i, b);
            Diff(a, b) +=
                std::conj(MIJ(a, i)) * Curvature_Quark_F2H1[i][b][m];
          }
        }
      }
    }

    res = FirstDerivativeOfEigenvalues(MassMatrix, Diff);

    for (std::size_t j = 0; j < res.size(); j++)
    {
      if (std::isnan(res.at(j)))
      {
        std::stringstream ss;
        ss << "MassMatrix = \n"
           << MassMatrix << "\nDiff = \n"
           << Diff << std::endl;
        ss << "Fermion Masses : ";
        for (std::size_t i = 0; i < NQuarks; i++)
          ss << std::sqrt(std::abs(res.at(i))) << sep;
        ss << std::endl;
        ss << "VEV fields : ";
        for (std::size_t i = 0; i < v.size(); i++)
          ss << v.at(i) << sep;
        ss << std::endl;

        for (std::size_t l = 0; l < NHiggs; l++)
        {

          ss << "Curvature_Quark * v an Higgs  =  :" << l << "\n";
          for (std::size_t a = 0; a < NQuarks; a++)
          {
            for (std::size_t i = 0; i < NQuarks; i++)
            {
              ss << Curvature_Quark_F2H1[a][i][l] * v[l] << sep;
            }
            ss << std::endl;
          }
          ss << "conj Curvature_Quark an Higgs = :" << l << "\n";
          for (std::size_t a = 0; a < NQuarks; a++)
          {
            for (std::size_t i = 0; i < NQuarks; i++)
            {
              ss << std::conj(Curvature_Quark_F2H1[a][i][l]) * v[l] << sep;
            }
            ss << std::endl;
          }
        }

        Logger::Write(LoggingLevel::Debug, ss.str());

        std::string retmessage = "Nan found in ";
        retmessage += __func__;
        retmessage += " at deriv number ";
        retmessage += std::to_string(j);
        retmessage += " and m = ";
        retmessage += std::to_string(m);
        throw std::runtime_error(retmessage);
      }
    }
  }

  return;
}

std::vector<double>
Class_Potential_Origin::LeptonMassesSquared(const std::vector<double> &v,
                                            const int &diff) const
{
  std::vector<double> res;
  LeptonMassesSquared(v, diff, res);
  return res;
}

void Class_Potential_Origin::LeptonMassesSquared(
    const std::vector<double> &v,
    const int &diff,
    std::vector<double> &res) const
{
  res.clear();
  res.reserve(diff > 0 ? 2 * NLepton : NLepton);
  // The fixed-size R2HDM path returns before MassMatrix is needed.  Leaving
  // both dynamic matrices initially empty avoids dead heap allocations while
  // preserving the fallback assignments and all floating-point operations.
  MatrixXcd MassMatrix, MIJ;
  double ZeroMass = std::pow(10, -10);

  // The validated fixed-size R2HDM path does not need a dynamic intermediate.
  // Assemble the same entries directly in stack storage, preserving the field
  // index order used by LeptonMassMatrix and retaining the complete 9x9 Eigen
  // solve.  Keep the analytic experiment's historical precedence when both
  // environment switches are set.
  static const bool UseDirectR2HDMLeptonFixed9Diff0 = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_LEPTON_FIXED9_DIFF0");
    return env != nullptr && env[0] == '1';
  }();
  static const bool UseDirectR2HDMLeptonIndexCache = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_LEPTON_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  static const bool UseLeptonAnalyticBeforeFixed = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_LEPTON_ANALYTIC_DIFF0");
    return env != nullptr && env[0] == '1';
  }();
  const bool useDirectFixed9 =
      UseDirectR2HDMLeptonFixed9Diff0 &&
      !UseLeptonAnalyticBeforeFixed && diff <= 0 &&
      Model == ModelID::ModelIDs::R2HDM && NLepton == 9 &&
      v.size() == NHiggs && SetCurvatureDone;
  if (useDirectFixed9)
  {
    using Matrix9cd = Eigen::Matrix<std::complex<double>, 9, 9>;
    Matrix9cd MIJFixed;
    const bool useDirectIndexCache =
        UseDirectR2HDMLeptonIndexCache && LeptonTensorIndexCacheReady &&
        LeptonMassFieldIndexCache.size() == NLepton * NLepton;
    for (std::size_t i = 0; i < NLepton; ++i)
    {
      for (std::size_t j = 0; j < NLepton; ++j)
      {
        MIJFixed(i, j) = Curvature_Lepton_F2[i][j];
        if (useDirectIndexCache)
        {
          const auto &fields =
              LeptonMassFieldIndexCache[i * NLepton + j];
          for (const std::size_t k : fields)
            MIJFixed(i, j) += Curvature_Lepton_F2H1[i][j][k] * v[k];
        }
        else
        {
          for (std::size_t k = 0; k < NHiggs; ++k)
            MIJFixed(i, j) += Curvature_Lepton_F2H1[i][j][k] * v[k];
        }
      }
    }

    const Matrix9cd MassMatrixFixed = MIJFixed.conjugate() * MIJFixed;
    SelfAdjointEigenSolver<Matrix9cd> es(MassMatrixFixed, EigenvaluesOnly);
    for (std::size_t i = 0; i < 9; ++i)
    {
      double tmp = es.eigenvalues()[static_cast<Eigen::Index>(i)];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
    return;
  }
  MIJ             = LeptonMassMatrix(v);

  // In the R2HDM lepton basis, the mass matrix consists of three independent
  // (charged-left, charged-right, neutrino-left) 3x3 blocks.  For each
  // generation the block has the form
  //
  //       [ 0  b  0 ]
  //       [ b  0  a ] ,
  //       [ 0  a  0 ]
  //
  // and therefore contributes 0, |a|^2 and |a|^2 + |b|^2 to M^dagger M.
  // Keep this opt-in: unlike Eigen's full solver, the closed form does not
  // reproduce its floating-point roundoff bit-for-bit.
  static const bool UseR2HDMLeptonAnalyticDiff0 = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_LEPTON_ANALYTIC_DIFF0");
    return env != nullptr && env[0] == '1';
  }();
  const bool useAnalyticDiff0 =
      UseR2HDMLeptonAnalyticDiff0 && diff <= 0 &&
      Model == ModelID::ModelIDs::R2HDM && NLepton == 9;

  if (useAnalyticDiff0)
  {
    // Check the complete sparse/symmetric structure before using the block
    // formula.  Any unexpected entry (including a non-symmetric allowed
    // entry) falls through to the original 9x9 Eigen path below.
    bool structureOk = true;
    for (std::size_t i = 0; i < NLepton && structureOk; ++i)
    {
      for (std::size_t j = 0; j < NLepton; ++j)
      {
        bool allowed = false;
        if (i < 6)
        {
          const std::size_t generation = i / 2;
          const std::size_t chargedLeft  = 2 * generation;
          const std::size_t chargedRight = chargedLeft + 1;
          const std::size_t neutrino     = 6 + generation;
          allowed = (i == chargedLeft && j == chargedRight) ||
                    (i == chargedRight && j == chargedLeft) ||
                    (i == chargedRight && j == neutrino);
        }
        else
        {
          const std::size_t generation = i - 6;
          if (generation < 3)
            allowed = j == 2 * generation + 1;
        }
        if (!allowed && MIJ(i, j) != std::complex<double>(0, 0))
        {
          structureOk = false;
          break;
        }
      }
    }

    for (std::size_t generation = 0; generation < 3 && structureOk;
         ++generation)
    {
      const std::size_t chargedLeft  = 2 * generation;
      const std::size_t chargedRight = chargedLeft + 1;
      const std::size_t neutrino     = 6 + generation;
      if (MIJ(chargedLeft, chargedRight) !=
              MIJ(chargedRight, chargedLeft) ||
          MIJ(chargedRight, neutrino) != MIJ(neutrino, chargedRight))
        structureOk = false;
    }

    if (structureOk)
    {
      std::vector<double> analyticMasses;
      analyticMasses.reserve(NLepton);
      for (std::size_t generation = 0; generation < 3; ++generation)
      {
        const std::size_t chargedLeft  = 2 * generation;
        const std::size_t chargedRight = chargedLeft + 1;
        const std::size_t neutrino     = 6 + generation;
        const std::complex<double> b = MIJ(chargedLeft, chargedRight);
        const std::complex<double> a = MIJ(chargedRight, neutrino);
        const double aSquared        = std::norm(a);
        const double bSquared        = std::norm(b);
        analyticMasses.push_back(0);
        analyticMasses.push_back(aSquared);
        analyticMasses.push_back(aSquared + bSquared);
      }
      std::sort(analyticMasses.begin(), analyticMasses.end());
      for (const double massSquared : analyticMasses)
      {
        if (std::abs(massSquared) < ZeroMass)
          res.push_back(0);
        else
          res.push_back(massSquared);
      }
      return;
    }
  }

  static const bool UseR2HDMLeptonIndexCache = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_LEPTON_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool useIndexCache =
      UseR2HDMLeptonIndexCache && Model == ModelID::ModelIDs::R2HDM &&
      NLepton == 9 && LeptonTensorIndexCacheReady &&
      LeptonDerivativeIndexCache.size() == NHiggs * NLepton * NLepton;

  // For the R2HDM default lepton sector the matrix dimension is fixed at 9.
  // Keep the original dynamic-matrix path as the default, but allow the
  // fixed-size Eigen path to be selected for the very hot diff <= 0 case.
  // The fixed matrix is populated from the same MIJ entries and uses the
  // same Hermitian product and EigenvaluesOnly solver; this preserves the
  // nine sorted eigenvalues and the existing ZeroMass treatment while
  // avoiding dynamic-size Eigen dispatch in this path.
  static const bool UseR2HDMLeptonFixed9Diff0 = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_LEPTON_FIXED9_DIFF0");
    return env != nullptr && env[0] == '1';
  }();
  const bool useFixed9Diff0 =
      UseR2HDMLeptonFixed9Diff0 && diff <= 0 &&
      Model == ModelID::ModelIDs::R2HDM && NLepton == 9;

  if (useFixed9Diff0)
  {
    using Matrix9cd = Eigen::Matrix<std::complex<double>, 9, 9>;
    Matrix9cd MIJFixed;
    for (Eigen::Index i = 0; i < 9; ++i)
      for (Eigen::Index j = 0; j < 9; ++j)
        MIJFixed(i, j) = MIJ(i, j);

    const Matrix9cd MassMatrixFixed = MIJFixed.conjugate() * MIJFixed;
    SelfAdjointEigenSolver<Matrix9cd> es(MassMatrixFixed, EigenvaluesOnly);
    for (std::size_t i = 0; i < 9; ++i)
    {
      double tmp = es.eigenvalues()[static_cast<Eigen::Index>(i)];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
    return;
  }

  MassMatrix = MIJ.conjugate() * MIJ;

  if (diff <= 0) // no temperature part here
  {
    SelfAdjointEigenSolver<MatrixXcd> es(MassMatrix, EigenvaluesOnly);
    for (std::size_t i = 0; i < NLepton; i++)
    {
      double tmp = es.eigenvalues().real()[i];
      if (std::abs(tmp) < ZeroMass)
        res.push_back(0);
      else
        res.push_back(tmp);
    }
  }
  else if (static_cast<size_t>(diff) <= NHiggs)
  {

    auto k         = diff - 1;
    MatrixXcd Diff = MatrixXcd::Zero(NLepton, NLepton);
    for (std::size_t I{0}; I < NLepton; ++I)
    {
      for (std::size_t J{0}; J < NLepton; ++J)
      {
        if (useIndexCache)
        {
          const auto &indices =
              LeptonDerivativeIndexCache[(k * NLepton + I) * NLepton + J];
          for (const auto &term : indices)
          {
            const std::size_t L = term.i;
            if (term.firstNonzero)
              Diff(I, J) +=
                  std::conj(Curvature_Lepton_F2H1[I][L][k]) * MIJ(L, J);
            if (term.secondNonzero)
              Diff(I, J) +=
                  std::conj(MIJ(I, L)) * Curvature_Lepton_F2H1[L][J][k];
          }
        }
        else
        {
          for (std::size_t L{0}; L < NLepton; ++L)
          {
            Diff(I, J) +=
                std::conj(Curvature_Lepton_F2H1[I][L][k]) * MIJ(L, J);
            Diff(I, J) +=
                std::conj(MIJ(I, L)) * Curvature_Lepton_F2H1[L][J][k];
          }
        }
      }
    }

    res = FirstDerivativeOfEigenvalues(MassMatrix, Diff);

    for (std::size_t j = 0; j < res.size(); j++)
    {
      if (std::isnan(res.at(j)))
      {
        std::stringstream ss;
        ss << "MassMatrix = \n"
           << MassMatrix << "\nDiff = \n"
           << Diff << std::endl;
        ss << "Fermion Masses : ";
        for (std::size_t i = 0; i < NLepton; i++)
          ss << std::sqrt(std::abs(res.at(i))) << sep;
        ss << std::endl;
        ss << "VEV fields : ";
        for (std::size_t i = 0; i < v.size(); i++)
          ss << v.at(i) << sep;
        ss << std::endl;

        for (std::size_t l = 0; l < NHiggs; l++)
        {

          ss << "Curvature_Lepton * v an Higgs  =  :" << l << "\n";
          for (std::size_t a = 0; a < NLepton; a++)
          {
            for (std::size_t i = 0; i < NLepton; i++)
            {
              ss << Curvature_Lepton_F2H1[a][i][l] * v[l] << sep;
            }
            ss << std::endl;
          }
          ss << "conj Curvature_Lepton an Higgs = :" << l << "\n";
          for (std::size_t a = 0; a < NLepton; a++)
          {
            for (std::size_t i = 0; i < NLepton; i++)
            {
              ss << std::conj(Curvature_Lepton_F2H1[a][i][l]) * v[l] << sep;
            }
            ss << std::endl;
          }
        }

        Logger::Write(LoggingLevel::Debug, ss.str());

        std::string retmessage = "Nan found in ";
        retmessage += __func__;
        retmessage += " at deriv number ";
        retmessage += std::to_string(j);
        retmessage += " and m = ";
        retmessage += std::to_string(diff - 1);
        throw std::runtime_error(retmessage);
      }
    }
  }

  return;
}

double Class_Potential_Origin::VTree(const std::vector<double> &v,
                                     int diff,
                                     bool ForceExplicitCalculation) const
{
  double res = 0;
  CalcGWProfiler::vtree_call();
  const bool profile_vtree = CalcGWProfiler::enabled();

  // Do not dispatch to the simplified virtual potential when the model has
  // explicitly disabled it.  The previous result was discarded in that
  // case, so this removes work without changing the explicit calculation.
  if (not ForceExplicitCalculation and UseVTreeSimplified)
  {
    res = VTreeSimplified(v);
    if (diff == 0)
    {

      return res;
    }
  }
  res = 0;

  static const bool UseR2HDMVTreeIndexCache = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_VTREE_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool useVTreeIndexCache =
      UseR2HDMVTreeIndexCache && Model == ModelID::ModelIDs::R2HDM &&
      NHiggs == 8 && VTreeTensorIndexCacheReady &&
      VTreeTensorIndexCache.size() == NHiggs;

  if (profile_vtree)
  {
    std::uint64_t terms = 0;
    if (diff == 0)
    {
      if (useVTreeIndexCache)
      {
        for (std::size_t i = 0; i < NHiggs; ++i)
        {
          if (v[i] == 0) continue;
          const auto &iTerms = VTreeTensorIndexCache[i];
          if (iTerms.l1Nonzero) ++terms;
          for (const auto &jTerm : iTerms.jTerms)
          {
            if (v[jTerm.j] == 0) continue;
            if (jTerm.l2Nonzero) ++terms;
            for (const auto &kTerm : jTerm.kTerms)
            {
              if (kTerm.l3Nonzero) ++terms;
              terms += kTerm.l4Nonzero.size();
            }
          }
        }
      }
      else
      {
        std::uint64_t active = 0;
        for (std::size_t i = 0; i < NHiggs; ++i)
          if (v[i] != 0) ++active;
        const std::uint64_t activeSquared = active * active;
        terms = active + activeSquared + activeSquared * NHiggs +
                activeSquared * NHiggs * NHiggs;
      }
    }
    else if (diff > 0 && static_cast<std::size_t>(diff) <= NHiggs)
    {
      if (useVTreeIndexCache)
      {
        const auto &iTerms = VTreeTensorIndexCache[diff - 1];
        terms = 1;
        for (const auto &jTerm : iTerms.jTerms)
        {
          if (jTerm.l2Nonzero) ++terms;
          for (const auto &kTerm : jTerm.kTerms)
          {
            if (kTerm.l3Nonzero) ++terms;
            terms += kTerm.l4Nonzero.size();
          }
        }
      }
      else
      {
        const std::uint64_t n = NHiggs;
        terms = 1 + n + n * n + n * n * n;
      }
    }
    CalcGWProfiler::vtree_terms(terms);
  }

  if (diff == 0)
  {
    if (useVTreeIndexCache)
    {
      for (std::size_t i = 0; i < NHiggs; ++i)
      {
        if (v[i] == 0) continue;
        const auto &iTerms = VTreeTensorIndexCache[i];
        if (iTerms.l1Nonzero)
          res += Curvature_Higgs_L1[i] * v[i];
        for (const auto &jTerm : iTerms.jTerms)
        {
          const std::size_t j = jTerm.j;
          if (v[j] == 0) continue;
          if (jTerm.l2Nonzero)
            res += 0.5 * Curvature_Higgs_L2[i][j] * v[i] * v[j];
          for (const auto &kTerm : jTerm.kTerms)
          {
            const std::size_t k = kTerm.k;
            if (kTerm.l3Nonzero)
              res += 1.0 / 6.0 * Curvature_Higgs_L3[i][j][k] * v[i] * v[j] *
                     v[k];
            for (const std::size_t l : kTerm.l4Nonzero)
            {
              res += 1.0 / 24.0 * Curvature_Higgs_L4[i][j][k][l] * v[i] *
                     v[j] * v[k] * v[l];
            }
          }
        }
      }
    }
    else
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        if (v[i] != 0)
        {
          res += Curvature_Higgs_L1[i] * v[i];
          for (std::size_t j = 0; j < NHiggs; j++)
          {
            if (v[j] != 0)
            {
              res += 0.5 * Curvature_Higgs_L2[i][j] * v[i] * v[j];
              for (std::size_t k = 0; k < NHiggs; k++)
              {
                res += 1.0 / 6.0 * Curvature_Higgs_L3[i][j][k] * v[i] *
                       v[j] * v[k];
                for (std::size_t l = 0; l < NHiggs; l++)
                {
                  res += 1.0 / 24.0 * Curvature_Higgs_L4[i][j][k][l] * v[i] *
                         v[j] * v[k] * v[l];
                }
              }
            }
          }
        }
      }
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t i = diff - 1;
    res           = Curvature_Higgs_L1[i];
    if (useVTreeIndexCache)
    {
      const auto &iTerms = VTreeTensorIndexCache[i];
      for (const auto &jTerm : iTerms.jTerms)
      {
        const std::size_t j = jTerm.j;
        if (jTerm.l2Nonzero)
          res += Curvature_Higgs_L2[i][j] * v[j];
        for (const auto &kTerm : jTerm.kTerms)
        {
          const std::size_t k = kTerm.k;
          if (kTerm.l3Nonzero)
            res += 0.5 * Curvature_Higgs_L3[i][j][k] * v[j] * v[k];
          for (const std::size_t l : kTerm.l4Nonzero)
            res += 1.0 / 6.0 * Curvature_Higgs_L4[i][j][k][l] * v[j] * v[k] *
                   v[l];
        }
      }
    }
    else
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        res += Curvature_Higgs_L2[i][j] * v[j];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          res += 0.5 * Curvature_Higgs_L3[i][j][k] * v[j] * v[k];
          for (std::size_t l = 0; l < NHiggs; l++)
          {
            res += 1.0 / 6.0 * Curvature_Higgs_L4[i][j][k][l] * v[j] * v[k] *
                   v[l];
          }
        }
      }
    }
  }

  return res;
}

double Class_Potential_Origin::CounterTerm(const std::vector<double> &v,
                                           int diff,
                                           bool ForceExplicitCalculation) const
{

  double res = 0;
  CalcGWProfiler::counterterm_call();
  const bool profile_counterterm = CalcGWProfiler::enabled();
  if (not ForceExplicitCalculation and UseVCounterSimplified)
  {
    res = VCounterSimplified(v);
    if (UseVCounterSimplified and diff == 0) return res;
  }

  res = 0;
  static const bool UseR2HDMCounterTermIndexCache = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_COUNTERTERM_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  static const bool UseR2HDMCounterTermFlatIndexCache = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_COUNTERTERM_FLAT_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_r2hdm_counterterm_index_cache =
      UseR2HDMCounterTermIndexCache &&
      Model == ModelID::ModelIDs::R2HDM && NHiggs == 8 &&
      CounterTermIndexCacheReady &&
      CounterTermIndexCache.size() == NHiggs;
  const bool use_r2hdm_counterterm_flat_index_cache =
      UseR2HDMCounterTermFlatIndexCache &&
      use_r2hdm_counterterm_index_cache && CounterTermFlatIndexCacheReady &&
      CounterTermFlatIndexCache.size() == NHiggs;
  if (profile_counterterm)
  {
    std::uint64_t terms = 0;
    if (diff == 0)
    {
      if (use_r2hdm_counterterm_index_cache)
      {
        for (const auto &iTerms : CounterTermIndexCache)
        {
          if (iTerms.l1Nonzero) ++terms;
          for (const auto &jTerms : iTerms.jTerms)
          {
            if (jTerms.l2Nonzero) ++terms;
            for (const auto &kTerms : jTerms.kTerms)
            {
              if (kTerms.l3Nonzero) ++terms;
              terms += kTerms.l4Nonzero.size();
            }
          }
        }
      }
      else
      {
        const std::uint64_t n = NHiggs;
        terms = n + n * n + n * n * n + n * n * n * n;
      }
    }
    else if (diff > 0 && static_cast<std::size_t>(diff) <= NHiggs)
    {
      if (use_r2hdm_counterterm_index_cache)
      {
        const auto &iTerms = CounterTermIndexCache[diff - 1];
        terms = 1;
        for (const auto &jTerms : iTerms.jTerms)
        {
          if (jTerms.l2Nonzero) ++terms;
          for (const auto &kTerms : jTerms.kTerms)
          {
            if (kTerms.l3Nonzero) ++terms;
            terms += kTerms.l4Nonzero.size();
          }
        }
      }
      else
      {
        const std::uint64_t n = NHiggs;
        terms = 1 + n + n * n + n * n * n;
      }
    }
    CalcGWProfiler::counterterm_terms(terms);
  }
  if (diff == 0)
  {
    if (use_r2hdm_counterterm_flat_index_cache)
    {
      for (std::size_t i = 0; i < NHiggs; ++i)
      {
        const auto &iTerms = CounterTermIndexCache[i];
        if (iTerms.l1Nonzero)
          res += Curvature_Higgs_CT_L1[i] * v[i];
        for (const auto &term : CounterTermFlatIndexCache[i])
        {
          const std::size_t j = term.j;
          const std::size_t k = term.k;
          if (term.rank == 2)
          {
            res += 0.5 * Curvature_Higgs_CT_L2[i][j] * v[i] * v[j];
          }
          else if (term.rank == 3)
          {
            res += 1.0 / 6.0 * Curvature_Higgs_CT_L3[i][j][k] * v[i] *
                   v[j] * v[k];
          }
          else
          {
            res += 1.0 / 24.0 * Curvature_Higgs_CT_L4[i][j][k][term.l] *
                   v[i] * v[j] * v[k] * v[term.l];
          }
        }
      }
    }
    else if (use_r2hdm_counterterm_index_cache)
    {
      for (std::size_t i = 0; i < NHiggs; ++i)
      {
        const auto &iTerms = CounterTermIndexCache[i];
        if (iTerms.l1Nonzero)
          res += Curvature_Higgs_CT_L1[i] * v[i];
        for (const auto &jTerms : iTerms.jTerms)
        {
          const std::size_t j = jTerms.j;
          if (jTerms.l2Nonzero)
            res += 0.5 * Curvature_Higgs_CT_L2[i][j] * v[i] * v[j];
          for (const auto &kTerms : jTerms.kTerms)
          {
            const std::size_t k = kTerms.k;
            if (kTerms.l3Nonzero)
              res += 1.0 / 6.0 * Curvature_Higgs_CT_L3[i][j][k] * v[i] *
                     v[j] * v[k];
            for (const std::size_t l : kTerms.l4Nonzero)
            {
              res += 1.0 / 24.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v[i] *
                     v[j] * v[k] * v[l];
            }
          }
        }
      }
    }
    else
    {
      for (std::size_t i = 0; i < NHiggs; i++)
      {
        res += Curvature_Higgs_CT_L1[i] * v[i];
        for (std::size_t j = 0; j < NHiggs; j++)
        {
          res += 0.5 * Curvature_Higgs_CT_L2[i][j] * v[i] * v[j];
          for (std::size_t k = 0; k < NHiggs; k++)
          {
            res += 1.0 / 6.0 * Curvature_Higgs_CT_L3[i][j][k] * v[i] * v[j] *
                   v[k];
            for (std::size_t l = 0; l < NHiggs; l++)
            {
              res += 1.0 / 24.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v[i] *
                     v[j] * v[k] * v[l];
            }
          }
        }
      }
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    std::size_t i = diff - 1;
    res           = Curvature_Higgs_CT_L1[i];
    if (use_r2hdm_counterterm_flat_index_cache)
    {
      for (const auto &term : CounterTermFlatIndexCache[i])
      {
        const std::size_t j = term.j;
        const std::size_t k = term.k;
        if (term.rank == 2)
        {
          res += Curvature_Higgs_CT_L2[i][j] * v[j];
        }
        else if (term.rank == 3)
        {
          res += 0.5 * Curvature_Higgs_CT_L3[i][j][k] * v[j] * v[k];
        }
        else
        {
          res += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][term.l] * v[j] *
                 v[k] * v[term.l];
        }
      }
    }
    else if (use_r2hdm_counterterm_index_cache)
    {
      const auto &iTerms = CounterTermIndexCache[i];
      for (const auto &jTerms : iTerms.jTerms)
      {
        const std::size_t j = jTerms.j;
        if (jTerms.l2Nonzero) res += Curvature_Higgs_CT_L2[i][j] * v[j];
        for (const auto &kTerms : jTerms.kTerms)
        {
          const std::size_t k = kTerms.k;
          if (kTerms.l3Nonzero)
            res += 0.5 * Curvature_Higgs_CT_L3[i][j][k] * v[j] * v[k];
          for (const std::size_t l : kTerms.l4Nonzero)
          {
            res += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v[j] *
                   v[k] * v[l];
          }
        }
      }
    }
    else
    {
      for (std::size_t j = 0; j < NHiggs; j++)
      {
        res += Curvature_Higgs_CT_L2[i][j] * v[j];
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          res += 0.5 * Curvature_Higgs_CT_L3[i][j][k] * v[j] * v[k];
          for (std::size_t l = 0; l < NHiggs; l++)
          {
            res += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v[j] *
                   v[k] * v[l];
          }
        }
      }
    }
  }

  return res;
}

double Class_Potential_Origin::VEff(const std::vector<double> &v,
                                    double Temp,
                                    int diff,
                                    int Order) const
{
  return VEff(v, Temp, diff, (Order == 0) ? Order::TreeLevel : Order::OneLoop);
}

double Class_Potential_Origin::VEff(const std::vector<double> &v,
                                    double Temp,
                                    int diff,
                                    const Order &order) const
{
  ProfileExactModelRepeat(CalcGWProfiler::ExactRepeatMetric::VEff, this, v,
                          Temp, diff, static_cast<int>(order));
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    return VEff(Transformedv, Temp, diff);
  }

  double resOut = 0;
  {
    CalcGWProfiler::ScopedTimer timer(
        CalcGWProfiler::TimingMetric::VTree);
    resOut = VTree(v, diff);
  }
  if (order != Order::TreeLevel and not UseTreeLevel)
  {
    {
      CalcGWProfiler::ScopedTimer timer(
          CalcGWProfiler::TimingMetric::CounterTerm);
      resOut += CounterTerm(v, diff);
    }
    resOut += V1Loop(v, Temp, diff);
  }
  // for(std::size_t i=0;i<NHiggs;i++) resOut +=
  // DebyeHiggs[i][i]*0.5*std::pow(v.at(i),2)*std::pow(Temp,2);
  return resOut;
}

double Class_Potential_Origin::V1Loop(const std::vector<double> &v,
                                      double Temp,
                                      int diff) const
{
  double res = 0;
  static const bool UseCachedProfileGate = []
  {
    const char *env = std::getenv("BSMPT_USE_CACHED_PROFILE_GATE");
    return env != nullptr && env[0] == '1';
  }();
  static const bool CachedProfileGate = CalcGWProfiler::timing_enabled();
  const auto ProfileV1LoopTiming = []() noexcept
  {
    return UseCachedProfileGate ? CachedProfileGate
                                : CalcGWProfiler::timing_enabled();
  };

  // V1Loop evaluates the same thermal prefactors for every particle.  Keep
  // this context opt-in and inactive at T=0 so the default path, including its
  // exact edge-case behavior, remains byte-for-byte unchanged.
  const bool use_thermal_context = UseV1LoopThermalContext() && Temp != 0;
  const V1LoopThermalConstants thermal_constants{
      use_thermal_context ? std::pow(Temp, 2) : 0.0,
      use_thermal_context ? std::pow(Temp, 3) : 0.0,
      use_thermal_context ? std::pow(Temp, 4) : 0.0,
      use_thermal_context ? 2 * std::pow(M_PI, 2) : 0.0};
  static const bool UseV1LoopMassContributionCache = []
  {
    const char *env =
        std::getenv("BSMPT_USE_V1LOOP_MASS_CONTRIBUTION_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  static const bool UseV1LoopLinearMassContributionCache = []
  {
    const char *env = std::getenv(
        "BSMPT_USE_V1LOOP_LINEAR_MASS_CONTRIBUTION_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_mass_contribution_hash_cache =
      UseV1LoopMassContributionCache;
  const bool use_mass_contribution_linear_cache =
      UseV1LoopLinearMassContributionCache &&
      !use_mass_contribution_hash_cache;
  std::unordered_map<BosonContributionCacheKey,
                     double,
                     BosonContributionCacheKeyHash>
      bosonContributionCache;
  std::unordered_map<FermionContributionCacheKey,
      double,
      FermionContributionCacheKeyHash>
      fermionContributionCache;
  std::vector<std::pair<BosonContributionCacheKey, double>>
      bosonLinearContributionCache;
  std::vector<std::pair<FermionContributionCacheKey, double>>
      fermionLinearContributionCache;
  if (use_mass_contribution_hash_cache)
  {
    bosonContributionCache.reserve(32);
    fermionContributionCache.reserve(32);
  }
  if (use_mass_contribution_linear_cache)
  {
    bosonLinearContributionCache.reserve(32);
    fermionLinearContributionCache.reserve(32);
  }
  const auto boson = [&](double MassSquared,
                         double temperature,
                         double cb,
                         int derivative,
                         double DevMassSquared = 0.0)
  {
    CalcGWProfiler::ScopedTimer timer(
        CalcGWProfiler::TimingMetric::BosonThermal, ProfileV1LoopTiming());
    if (use_mass_contribution_hash_cache ||
        use_mass_contribution_linear_cache)
    {
      const BosonContributionCacheKey key{
          ExactDoubleBits(MassSquared),
          ExactDoubleBits(temperature),
          ExactDoubleBits(cb),
          ExactDoubleBits(DevMassSquared),
          derivative};
      if (use_mass_contribution_hash_cache)
      {
        const auto cached = bosonContributionCache.find(key);
        if (cached != bosonContributionCache.end()) return cached->second;
      }
      else
      {
        for (const auto &entry : bosonLinearContributionCache)
          if (entry.first == key) return entry.second;
      }
      const double value =
          use_thermal_context
              ? BosonWithV1LoopThermalContext(*this,
                                               MassSquared,
                                               temperature,
                                               cb,
                                               derivative,
                                               DevMassSquared,
                                               thermal_constants)
              : this->boson(MassSquared,
                            temperature,
                            cb,
                            derivative,
                            DevMassSquared);
      if (use_mass_contribution_hash_cache)
        bosonContributionCache.emplace(key, value);
      else
        bosonLinearContributionCache.emplace_back(key, value);
      return value;
    }
    return use_thermal_context
               ? BosonWithV1LoopThermalContext(*this,
                                                MassSquared,
                                                temperature,
                                                cb,
                                                derivative,
                                                DevMassSquared,
                                                thermal_constants)
               : this->boson(MassSquared,
                             temperature,
                             cb,
                             derivative,
                             DevMassSquared);
  };
  const auto fermion = [&](double MassSquared,
                           double temperature,
                           int derivative)
  {
    CalcGWProfiler::ScopedTimer timer(
        CalcGWProfiler::TimingMetric::FermionThermal, ProfileV1LoopTiming());
    if (use_mass_contribution_hash_cache ||
        use_mass_contribution_linear_cache)
    {
      const FermionContributionCacheKey key{ExactDoubleBits(MassSquared),
                                            ExactDoubleBits(temperature),
                                            derivative};
      if (use_mass_contribution_hash_cache)
      {
        const auto cached = fermionContributionCache.find(key);
        if (cached != fermionContributionCache.end()) return cached->second;
      }
      else
      {
        for (const auto &entry : fermionLinearContributionCache)
          if (entry.first == key) return entry.second;
      }
      const double value =
          use_thermal_context
              ? FermionWithV1LoopThermalContext(*this,
                                                 MassSquared,
                                                 temperature,
                                                 derivative,
                                                 thermal_constants)
              : this->fermion(MassSquared, temperature, derivative);
      if (use_mass_contribution_hash_cache)
        fermionContributionCache.emplace(key, value);
      else
        fermionLinearContributionCache.emplace_back(key, value);
      return value;
    }
    return use_thermal_context
               ? FermionWithV1LoopThermalContext(*this,
                                                  MassSquared,
                                                  temperature,
                                                  derivative,
                                                  thermal_constants)
               : this->fermion(MassSquared, temperature, derivative);
  };

  /**
   * In case of diff != 0 the mass vectors will directly have the derivatives of
   * the masses
   */
  std::vector<double> HiggsMassesVec, QuarkMassesVec, GaugeMassesVec,
      LeptonMassesVec, HiggsMassesZeroTempVec, GaugeMassesZeroTempVec;
  static const bool UseV1LoopMassVectorReserve = []
  {
    const char *env = std::getenv("BSMPT_USE_V1LOOP_MASS_VECTOR_RESERVE");
    return env != nullptr && env[0] == '1';
  }();
  static const bool UseR2HDMHiggsPair = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_HIGGS_PAIR");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_r2hdm_higgs_pair =
      UseR2HDMHiggsPair && !C_UseParwani && diff == 0 &&
      Model == ModelID::ModelIDs::R2HDM && NHiggs == 8 && nVEV == 4;
  static const bool UseR2HDMHiggsPairDiff = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_HIGGS_PAIR_DIFF");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_r2hdm_higgs_pair_diff =
      UseR2HDMHiggsPairDiff && !C_UseParwani && diff > 0 &&
      diff <= static_cast<int>(NHiggs) &&
      Model == ModelID::ModelIDs::R2HDM && NHiggs == 8 && nVEV == 4;
  static const bool UseR2HDMGaugePair = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_GAUGE_PAIR");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_r2hdm_gauge_pair =
      UseR2HDMGaugePair && !C_UseParwani && diff >= 0 &&
      diff <= static_cast<int>(NHiggs) &&
      Model == ModelID::ModelIDs::R2HDM && NHiggs == 8 && nVEV == 4 &&
      NGauge == 4;
  static const bool UseR2HDMGaugeIndexCache = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_GAUGE_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_r2hdm_gauge_index_cache =
      use_r2hdm_gauge_pair && UseR2HDMGaugeIndexCache &&
      GaugeTensorIndexCacheReady &&
      GaugeTensorIndexCache.size() == NGauge * NGauge;
  static const bool UseR2HDMFixedEigenSolver = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_FIXED_EIGENSOLVER");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_r2hdm_fixed_eigensolver =
      UseR2HDMFixedEigenSolver && diff == 0 && !C_UseParwani &&
      Model == ModelID::ModelIDs::R2HDM && NHiggs == 8 && NGauge == 4;
  static const bool UseR2HDMPairStackMatrices = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_PAIR_STACK_MATRICES");
    return env != nullptr && env[0] == '1';
  }();
  const bool use_r2hdm_pair_stack_matrices =
      UseR2HDMPairStackMatrices && use_r2hdm_higgs_pair &&
      use_r2hdm_gauge_pair;

  // In the Arnold--Espinosa, diff == 0 path the R2HDM Higgs spectrum needs
  // both M_H(v,T) and M_H(v,0).  Their field-dependent polynomial part is
  // identical; only the Debye term differs.  Build M_H(v,0) once, add the
  // Debye term to a copy, and keep two independent eigensolver calls.  The
  // opt-in guard deliberately falls back to the original path for every
  // other model, derivative, or resummation prescription.
  {
  CalcGWProfiler::ScopedTimer higgs_timer(
      CalcGWProfiler::TimingMetric::HiggsMasses, ProfileV1LoopTiming());
  if (use_r2hdm_higgs_pair)
  {
    if (use_r2hdm_pair_stack_matrices && v.size() == NHiggs)
    {
      // Assemble the two R2HDM Higgs matrices in stack storage, then reuse a
      // single dynamic workspace for the unchanged dynamic eigensolver.  The
      // contraction order is the same as HiggsMassMatrix; only allocation
      // and matrix-copy traffic are changed by this opt-in experiment.
      Matrix<double, 8, 8> HiggsMassZero;
      HiggsMassZero.setZero();
      for (std::size_t i = 0; i < NHiggs; ++i)
      {
        for (std::size_t j = i; j < NHiggs; ++j)
        {
          HiggsMassZero(i, j) = Curvature_Higgs_L2[i][j];
          for (std::size_t k = 0; k < NHiggs; ++k)
          {
            HiggsMassZero(i, j) += Curvature_Higgs_L3[i][j][k] * v[k];
            for (std::size_t l = 0; l < NHiggs; ++l)
            {
              HiggsMassZero(i, j) +=
                  0.5 * Curvature_Higgs_L4[i][j][k][l] * v[k] * v[l];
            }
          }
        }
      }
      for (std::size_t i = 1; i < NHiggs; ++i)
      {
        for (std::size_t j = 0; j < i; ++j)
          HiggsMassZero(i, j) = HiggsMassZero(j, i);
      }

      Matrix<double, 8, 8> HiggsMassThermal = HiggsMassZero;
      if (Temp != 0)
      {
        for (std::size_t i = 0; i < NHiggs; ++i)
        {
          for (std::size_t j = i; j < NHiggs; ++j)
          {
            HiggsMassThermal(i, j) +=
                DebyeHiggs[i][j] * std::pow(Temp, 2);
          }
        }
        for (std::size_t i = 1; i < NHiggs; ++i)
        {
          for (std::size_t j = 0; j < i; ++j)
            HiggsMassThermal(i, j) = HiggsMassThermal(j, i);
        }
      }

      MatrixXd eigensolverMatrix(NHiggs, NHiggs);
      const auto StackHiggsEigenvalues =
      [&](const Matrix<double, 8, 8> &massMatrix)
      {
        std::vector<double> result;
        if (UseV1LoopMassVectorReserve) result.reserve(NHiggs);
        const double ZeroMass = std::pow(10, -5);
        eigensolverMatrix = massMatrix;
        SelfAdjointEigenSolver<MatrixXd> es(eigensolverMatrix,
                                            EigenvaluesOnly);
        const auto EV = es.eigenvalues();
        for (std::size_t i = 0; i < NHiggs; ++i)
        {
          if (std::abs(EV[i]) < ZeroMass)
            result.push_back(0);
          else
            result.push_back(EV[i]);
        }
        return result;
      };
      HiggsMassesVec         = StackHiggsEigenvalues(HiggsMassThermal);
      HiggsMassesZeroTempVec = StackHiggsEigenvalues(HiggsMassZero);
    }
    else
    {
    MatrixXd HiggsMassZero = HiggsMassMatrix(v, 0);
    const auto HiggsEigenvalues = [&](const MatrixXd &massMatrix)
    {
      std::vector<double> result;
      if (UseV1LoopMassVectorReserve) result.reserve(NHiggs);
      const double ZeroMass = std::pow(10, -5);
      if (use_r2hdm_fixed_eigensolver)
      {
        Matrix<double, 8, 8> fixedMatrix = massMatrix;
        SelfAdjointEigenSolver<Matrix<double, 8, 8>> es(fixedMatrix,
                                                        EigenvaluesOnly);
        const auto EV = es.eigenvalues();
        for (std::size_t i = 0; i < NHiggs; ++i)
        {
          if (std::abs(EV[i]) < ZeroMass)
            result.push_back(0);
          else
            result.push_back(EV[i]);
        }
      }
      else
      {
        SelfAdjointEigenSolver<MatrixXd> es(massMatrix, EigenvaluesOnly);
        const auto EV = es.eigenvalues();
        for (std::size_t i = 0; i < NHiggs; ++i)
        {
          if (std::abs(EV[i]) < ZeroMass)
            result.push_back(0);
          else
            result.push_back(EV[i]);
        }
      }
      return result;
    };

    MatrixXd HiggsMassThermal = HiggsMassZero;
    if (Temp != 0)
    {
      for (std::size_t i = 0; i < NHiggs; ++i)
      {
        for (std::size_t j = i; j < NHiggs; ++j)
          HiggsMassThermal(i, j) +=
              DebyeHiggs[i][j] * std::pow(Temp, 2);
      }
      for (std::size_t i = 1; i < NHiggs; ++i)
        for (std::size_t j = 0; j < i; ++j)
          HiggsMassThermal(i, j) = HiggsMassThermal(j, i);
    }

    HiggsMassesVec         = HiggsEigenvalues(HiggsMassThermal);
    HiggsMassesZeroTempVec = HiggsEigenvalues(HiggsMassZero);
    }
  }
  else if (use_r2hdm_higgs_pair_diff)
  {
    const MatrixXd HiggsMassZero = HiggsMassMatrix(v, 0);
    MatrixXd HiggsMassThermal   = HiggsMassZero;
    if (Temp != 0)
    {
      for (std::size_t i = 0; i < NHiggs; ++i)
      {
        for (std::size_t j = i; j < NHiggs; ++j)
        {
          HiggsMassThermal(i, j) +=
              DebyeHiggs[i][j] * std::pow(Temp, 2);
        }
      }
      for (std::size_t i = 1; i < NHiggs; ++i)
      {
        for (std::size_t j = 0; j < i; ++j)
        {
          HiggsMassThermal(i, j) = HiggsMassThermal(j, i);
        }
      }
    }

    // The field derivative is independent of temperature.  Assemble it once,
    // but keep separate complex copies as inputs to the two independent
    // FirstDerivativeOfEigenvalues calls below.
    MatrixXd HiggsMassDerivative(NHiggs, NHiggs);
    const std::size_t x0 = static_cast<std::size_t>(diff - 1);
    for (std::size_t i = 0; i < NHiggs; ++i)
    {
      for (std::size_t j = 0; j < NHiggs; ++j)
      {
        HiggsMassDerivative(i, j) = Curvature_Higgs_L3[i][j][x0];
        for (std::size_t k = 0; k < NHiggs; ++k)
        {
          HiggsMassDerivative(i, j) +=
              Curvature_Higgs_L4[i][j][x0][k] * v[k];
        }
      }
    }

    MatrixXcd MassCastThermal(NHiggs, NHiggs);
    MassCastThermal = HiggsMassThermal;
    MatrixXcd DiffCastThermal(NHiggs, NHiggs);
    DiffCastThermal = HiggsMassDerivative;
    HiggsMassesVec =
        FirstDerivativeOfEigenvalues(MassCastThermal, DiffCastThermal);

    MatrixXcd MassCastZero(NHiggs, NHiggs);
    MassCastZero = HiggsMassZero;
    MatrixXcd DiffCastZero(NHiggs, NHiggs);
    DiffCastZero = HiggsMassDerivative;
    HiggsMassesZeroTempVec =
        FirstDerivativeOfEigenvalues(MassCastZero, DiffCastZero);
  }
  else
  {
    HiggsMassesVec = HiggsMassesSquared(v, Temp, diff);
  }
  }
  {
  CalcGWProfiler::ScopedTimer gauge_timer(
      CalcGWProfiler::TimingMetric::GaugeMasses, ProfileV1LoopTiming());
  if (use_r2hdm_gauge_pair)
  {
    if (UseV1LoopMassVectorReserve)
    {
      GaugeMassesVec.reserve(NGauge);
      GaugeMassesZeroTempVec.reserve(NGauge);
    }
    if (use_r2hdm_pair_stack_matrices && v.size() == NHiggs)
    {
      // Keep the gauge assembly in fixed stack storage and reuse one dynamic
      // matrix for both unchanged dynamic eigensolver calls.  This avoids a
      // second MatrixXd allocation while retaining the original solver and
      // the original upper/full-triangle policy.
      static const bool UseUpperOnly = []
      {
        const char *env = std::getenv("BSMPT_SAFE_GAUGE_UPPER_ONLY");
        return env != nullptr && env[0] == '1';
      }();
      Matrix<double, 4, 4> GaugeMassZero;
      GaugeMassZero.setZero();
      for (std::size_t a = 0; a < NGauge; ++a)
      {
        const std::size_t bBegin = UseUpperOnly ? a : 0;
        for (std::size_t b = bBegin; b < NGauge; ++b)
        {
          if (use_r2hdm_gauge_index_cache)
          {
            const auto &terms = GaugeTensorIndexCache[a * NGauge + b].mass;
            for (const auto &term : terms)
            {
              GaugeMassZero(a, b) +=
                  0.5 * Curvature_Gauge_G2H2[a][b][term.first][term.second] *
                  v.at(term.first) * v.at(term.second);
            }
          }
          else
          {
            for (std::size_t i = 0; i < NHiggs; ++i)
            {
              for (std::size_t j = 0; j < NHiggs; ++j)
              {
                GaugeMassZero(a, b) +=
                    0.5 * Curvature_Gauge_G2H2[a][b][i][j] * v.at(i) * v.at(j);
              }
            }
          }
        }
      }
      for (std::size_t a = 1; a < NGauge; ++a)
      {
        for (std::size_t b = 0; b < a; ++b)
          GaugeMassZero(a, b) = GaugeMassZero(b, a);
      }

      Matrix<double, 4, 4> GaugeMassThermal = GaugeMassZero;
      if (Temp != 0)
      {
        for (std::size_t a = 0; a < NGauge; ++a)
        {
          const std::size_t bBegin = UseUpperOnly ? a : 0;
          for (std::size_t b = bBegin; b < NGauge; ++b)
            GaugeMassThermal(a, b) +=
                DebyeGauge[a][b] * std::pow(Temp, 2);
        }
        for (std::size_t a = 1; a < NGauge; ++a)
        {
          for (std::size_t b = 0; b < a; ++b)
            GaugeMassThermal(a, b) = GaugeMassThermal(b, a);
        }
      }

      MatrixXd eigensolverMatrix(NGauge, NGauge);
      const auto StackGaugeEigenvalues =
      [&](const Matrix<double, 4, 4> &massMatrix)
      {
        std::vector<double> result;
        const double ZeroMass = std::pow(10, -5);
        eigensolverMatrix = massMatrix;
        SelfAdjointEigenSolver<MatrixXd> es(eigensolverMatrix,
                                            EigenvaluesOnly);
        for (std::size_t i = 0; i < NGauge; ++i)
        {
          double tmp = es.eigenvalues()[i];
          if (std::abs(tmp) < ZeroMass) tmp = 0;
          result.push_back(tmp);
        }
        return result;
      };
      GaugeMassesVec         = StackGaugeEigenvalues(GaugeMassThermal);
      GaugeMassesZeroTempVec = StackGaugeEigenvalues(GaugeMassZero);
    }
    else
    {
    // GaugeMassesSquared(v,T) and GaugeMassesSquared(v,0) have the same
    // field-dependent matrix.  Keep the original upper/full-triangle policy,
    // copy the zero-temperature matrix, add Debye terms to the thermal copy,
    // and retain separate eigensolver calls for the two matrices.
    static const bool UseUpperOnly = []
    {
      const char *env = std::getenv("BSMPT_SAFE_GAUGE_UPPER_ONLY");
      return env != nullptr && env[0] == '1';
    }();
    MatrixXd GaugeMassZero(NGauge, NGauge);
    for (std::size_t a = 0; a < NGauge; ++a)
    {
      const std::size_t bBegin = UseUpperOnly ? a : 0;
      for (std::size_t b = bBegin; b < NGauge; ++b)
      {
        GaugeMassZero(a, b) = 0;
        if (use_r2hdm_gauge_index_cache)
        {
          const auto &terms = GaugeTensorIndexCache[a * NGauge + b].mass;
          for (const auto &term : terms)
          {
            GaugeMassZero(a, b) +=
                0.5 * Curvature_Gauge_G2H2[a][b][term.first][term.second] *
                v.at(term.first) * v.at(term.second);
          }
        }
        else
        {
          for (std::size_t i = 0; i < NHiggs; ++i)
          {
            for (std::size_t j = 0; j < NHiggs; ++j)
            {
              GaugeMassZero(a, b) +=
                  0.5 * Curvature_Gauge_G2H2[a][b][i][j] * v.at(i) * v.at(j);
            }
          }
        }
      }
    }
    for (std::size_t a = 1; a < NGauge; ++a)
    {
      for (std::size_t b = 0; b < a; ++b)
      {
        GaugeMassZero(a, b) = GaugeMassZero(b, a);
      }
    }

    MatrixXd GaugeMassThermal = GaugeMassZero;
    if (Temp != 0)
    {
      for (std::size_t a = 0; a < NGauge; ++a)
      {
        const std::size_t bBegin = UseUpperOnly ? a : 0;
        for (std::size_t b = bBegin; b < NGauge; ++b)
        {
          GaugeMassThermal(a, b) +=
              DebyeGauge[a][b] * std::pow(Temp, 2);
        }
      }
      for (std::size_t a = 1; a < NGauge; ++a)
      {
        for (std::size_t b = 0; b < a; ++b)
        {
          GaugeMassThermal(a, b) = GaugeMassThermal(b, a);
        }
      }
    }

    const double ZeroMass = std::pow(10, -5);
    if (diff == 0)
    {
      if (use_r2hdm_fixed_eigensolver)
      {
        Matrix<double, 4, 4> fixedThermal = GaugeMassThermal;
        SelfAdjointEigenSolver<Matrix<double, 4, 4>> esThermal(
            fixedThermal, EigenvaluesOnly);
        for (std::size_t i = 0; i < NGauge; ++i)
        {
          double tmp = esThermal.eigenvalues()[i];
          if (std::abs(tmp) < ZeroMass) tmp = 0;
          GaugeMassesVec.push_back(tmp);
        }
        Matrix<double, 4, 4> fixedZero = GaugeMassZero;
        SelfAdjointEigenSolver<Matrix<double, 4, 4>> esZero(fixedZero,
                                                            EigenvaluesOnly);
        for (std::size_t i = 0; i < NGauge; ++i)
        {
          double tmp = esZero.eigenvalues()[i];
          if (std::abs(tmp) < ZeroMass) tmp = 0;
          GaugeMassesZeroTempVec.push_back(tmp);
        }
      }
      else
      {
        SelfAdjointEigenSolver<MatrixXd> esThermal(GaugeMassThermal,
                                                    EigenvaluesOnly);
        for (std::size_t i = 0; i < NGauge; ++i)
        {
          double tmp = esThermal.eigenvalues()[i];
          if (std::abs(tmp) < ZeroMass) tmp = 0;
          GaugeMassesVec.push_back(tmp);
        }
        SelfAdjointEigenSolver<MatrixXd> esZero(GaugeMassZero,
                                                 EigenvaluesOnly);
        for (std::size_t i = 0; i < NGauge; ++i)
        {
          double tmp = esZero.eigenvalues()[i];
          if (std::abs(tmp) < ZeroMass) tmp = 0;
          GaugeMassesZeroTempVec.push_back(tmp);
        }
      }
    }
    else // diff > 0 by the opt-in guard above
    {
      const std::size_t x0 = static_cast<std::size_t>(diff - 1);
      MatrixXd GaugeMassDerivative(NGauge, NGauge);
      GaugeMassDerivative = MatrixXd::Zero(NGauge, NGauge);
      for (std::size_t a = 0; a < NGauge; ++a)
      {
        for (std::size_t b = 0; b < NGauge; ++b)
        {
          if (use_r2hdm_gauge_index_cache)
          {
            const auto &terms =
                GaugeTensorIndexCache[a * NGauge + b].derivativeByField[x0];
            for (const std::size_t j : terms)
            {
              GaugeMassDerivative(a, b) +=
                  Curvature_Gauge_G2H2[a][b][x0][j] * v[j];
            }
          }
          else
          {
            for (std::size_t j = 0; j < NHiggs; ++j)
            {
              GaugeMassDerivative(a, b) +=
                  Curvature_Gauge_G2H2[a][b][x0][j] * v[j];
            }
          }
        }
      }

      MatrixXcd MassCastThermal(NGauge, NGauge);
      MassCastThermal = GaugeMassThermal;
      MatrixXcd DiffCastThermal(NGauge, NGauge);
      DiffCastThermal = GaugeMassDerivative;
      GaugeMassesVec =
          FirstDerivativeOfEigenvalues(MassCastThermal, DiffCastThermal);

      MatrixXcd MassCastZero(NGauge, NGauge);
      MassCastZero = GaugeMassZero;
      MatrixXcd DiffCastZero(NGauge, NGauge);
      DiffCastZero = GaugeMassDerivative;
      GaugeMassesZeroTempVec =
          FirstDerivativeOfEigenvalues(MassCastZero, DiffCastZero);
    }
    }
  }
  else
  {
    GaugeMassesVec         = GaugeMassesSquared(v, Temp, diff);
    GaugeMassesZeroTempVec = GaugeMassesSquared(v, 0, diff);
  }
  }
  {
    CalcGWProfiler::ScopedTimer timer(
        CalcGWProfiler::TimingMetric::QuarkMasses, ProfileV1LoopTiming());
    QuarkMassesVec = QuarkMassesSquared(v, diff);
  }
  {
    CalcGWProfiler::ScopedTimer timer(
        CalcGWProfiler::TimingMetric::LeptonMasses, ProfileV1LoopTiming());
    LeptonMassesVec = LeptonMassesSquared(v, diff);
  }

  // Length = N if diff = 0 else Length = 2N (mass_i, deriv_mass_i)
  const size_t NMultiplier         = (diff == 0) ? 1 : 2;
  const size_t NMultiplierFermions = (diff <= 0) ? 1 : 2;
  if (HiggsMassesVec.size() != NHiggs * NMultiplier)
    throw("Missmatch NHiggs [V1Loop]");
  if (GaugeMassesVec.size() != NGauge * NMultiplier)
    throw("Missmatch NGauge [V1Loop]");
  if (GaugeMassesZeroTempVec.size() != NGauge * NMultiplier)
    throw("Missmatch NGauge T = 0 [V1Loop]");
  if (QuarkMassesVec.size() != NQuarks * NMultiplierFermions)
    throw("Missmatch NQuarks [V1Loop]");
  if (LeptonMassesVec.size() != NLepton * NMultiplierFermions)
    throw("Missmatch NLepton [V1Loop]");

  if (diff == 0)
  {
    if (C_UseParwani)
    {
      // CW + Jbf with (field dependent + thermal) masses
      for (std::size_t k = 0; k < NHiggs; k++)
        res += boson(HiggsMassesVec[k], Temp, C_CWcbHiggs, 0);
      for (std::size_t k = 0; k < NGauge; k++)
        res += boson(GaugeMassesVec[k], Temp, C_CWcbGB, 0);
      for (std::size_t k = 0; k < NGauge; k++)
        res += 2 * boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, 0);
      for (std::size_t k = 0; k < NQuarks; k++)
        res += -6 * fermion(QuarkMassesVec[k], Temp, 0);
      for (std::size_t k = 0; k < NLepton; k++)
        res += -2 * fermion(LeptonMassesVec[k], Temp, 0);
    }
    else
    {
      // CW + Jbf with field dependent masses
      if (!use_r2hdm_higgs_pair && !use_r2hdm_higgs_pair_diff)
        HiggsMassesZeroTempVec = HiggsMassesSquared(v, 0);
      if (HiggsMassesVec.size() != NHiggs * NMultiplier)
        throw("Missmatch NHiggs T = 0 [V1Loop]");
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += boson(HiggsMassesZeroTempVec[k], Temp, C_CWcbHiggs, 0);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 3 * boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, 0);
      }
      double AddContQuark = 0;
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        AddContQuark += -2 * fermion(QuarkMassesVec[k], Temp, 0);
      }
      res += NColour * AddContQuark;
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * fermion(LeptonMassesVec[k], Temp, 0);
      }
      // Vdaisy
      double VDebye = 0;
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        if (HiggsMassesVec[k] > 0) VDebye += std::pow(HiggsMassesVec[k], 1.5);
        if (HiggsMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(HiggsMassesZeroTempVec[k], 1.5);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        if (GaugeMassesVec[k] > 0) VDebye += std::pow(GaugeMassesVec[k], 1.5);
        if (GaugeMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(GaugeMassesZeroTempVec[k], 1.5);
      }

      VDebye *= -Temp / (12 * M_PI);
      res += VDebye;
    }
  }
  else if (diff > 0 and static_cast<size_t>(diff) <= NHiggs)
  {
    if (C_UseParwani)
    {
      // dm^2/dphi * d(CW + Jbf)/dm^2
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += HiggsMassesVec[k + NHiggs] *
               boson(HiggsMassesVec[k], Temp, C_CWcbHiggs, diff);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += GaugeMassesVec[k + NGauge] *
               boson(GaugeMassesVec[k], Temp, C_CWcbGB, diff);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 2 * GaugeMassesZeroTempVec[k + NGauge] *
               boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, diff);
      }
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        res += -6 * QuarkMassesVec[k + NQuarks] *
               fermion(QuarkMassesVec[k], Temp, diff);
      }
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * LeptonMassesVec[k + NLepton] *
               fermion(LeptonMassesVec[k], Temp, diff);
      }
    }
    else
    {
      // dm^2/dphi * d(CW and Jbf)/dm^2
      if (!use_r2hdm_higgs_pair_diff)
        HiggsMassesZeroTempVec = HiggsMassesSquared(v, 0, diff);
      if (HiggsMassesVec.size() != NHiggs * NMultiplier)
        throw("Missmatch NHiggs T = 0 [V1Loop]");
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += HiggsMassesZeroTempVec[k + NHiggs] *
               boson(HiggsMassesZeroTempVec[k], Temp, C_CWcbHiggs, diff);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 3 * GaugeMassesZeroTempVec[k + NGauge] *
               boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, diff);
      }
      double AddContQuark = 0;
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        AddContQuark += -2 * QuarkMassesVec[k + NQuarks] *
                        fermion(QuarkMassesVec[k], Temp, diff);
      }
      for (std::size_t k = 0; k < NColour; k++)
        res += AddContQuark;
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * LeptonMassesVec[k + NLepton] *
               fermion(LeptonMassesVec[k], Temp, diff);
      }
      // dm^2/dphi * d(Vdaisy)/dm^2
      double VDebye = 0;
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        if (HiggsMassesVec[k] > 0)
        {
          VDebye += 1.5 * HiggsMassesVec[k + NHiggs] *
                    std::pow(HiggsMassesVec[k], 0.5);
        }
        if (HiggsMassesZeroTempVec[k] > 0)
        {
          VDebye += -1.5 * HiggsMassesZeroTempVec[k + NHiggs] *
                    std::pow(HiggsMassesZeroTempVec[k], 0.5);
        }
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        if (GaugeMassesVec[k] > 0)
        {
          VDebye += 1.5 * GaugeMassesVec[k + NGauge] *
                    std::pow(GaugeMassesVec[k], 0.5);
        }
        if (GaugeMassesZeroTempVec[k] > 0)
        {
          VDebye += -1.5 * GaugeMassesZeroTempVec[k + NGauge] *
                    std::pow(GaugeMassesZeroTempVec[k], 0.5);
        }
      }

      VDebye *= -Temp / (12 * M_PI);
      res += VDebye;
    }
  }
  else if (diff == -1)
  {
    if (C_UseParwani)
    {
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        // dm^2/dT * d(CW)/dm^2
        res += HiggsMassesVec[k + NHiggs] * boson(HiggsMassesVec[k],
                                                  Temp,
                                                  C_CWcbHiggs,
                                                  1,
                                                  HiggsMassesVec[k + NHiggs]);
        // d(Jbf)/dT
        res += boson(HiggsMassesVec[k],
                     Temp,
                     C_CWcbHiggs,
                     -1,
                     HiggsMassesVec[k + NHiggs]);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += GaugeMassesVec[k + NGauge] * boson(GaugeMassesVec[k],
                                                  Temp,
                                                  C_CWcbGB,
                                                  1,
                                                  GaugeMassesVec[k + NGauge]);
        res += boson(
            GaugeMassesVec[k], Temp, C_CWcbGB, -1, GaugeMassesVec[k + NGauge]);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 2 * boson(GaugeMassesZeroTempVec[k], Temp, C_CWcbGB, -1, 0.0);
      }
      for (std::size_t k = 0; k < NQuarks; k++)
      {
        res += -6 * fermion(QuarkMassesVec[k], Temp, -1);
      }
      for (std::size_t k = 0; k < NLepton; k++)
      {
        res += -2 * fermion(LeptonMassesVec[k], Temp, -1);
      }
    }
    else
    {
      // d(Jbf)/dT
      HiggsMassesZeroTempVec = HiggsMassesSquared(v, 0, diff);
      if (HiggsMassesVec.size() != NHiggs * NMultiplier)
        throw("Missmatch NHiggs T = 0 [V1Loop]");
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        res += boson(HiggsMassesZeroTempVec[k],
                     Temp,
                     C_CWcbHiggs,
                     -1,
                     0.0 /* 0 added for clarity */);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        res += 3 * boson(GaugeMassesZeroTempVec[k],
                         Temp,
                         C_CWcbGB,
                         -1,
                         0.0 /* 0 added for clarity */);
      }
      double AddContQuark = 0;
      for (std::size_t k = 0; k < NQuarks; k++)
        AddContQuark += -2 * fermion(QuarkMassesVec[k], Temp, -1);
      for (std::size_t k = 0; k < NColour; k++)
        res += AddContQuark;
      for (std::size_t k = 0; k < NLepton; k++)
        res += -2 * fermion(LeptonMassesVec[k], Temp, -1);

      // d(Vdaisy)/dT
      double VDebye = 0;
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        if (HiggsMassesVec[k] > 0)
        {
          // Explicit T part
          VDebye += std::pow(HiggsMassesVec[k], 1.5);
          // Thermal mass contribution
          VDebye += 1.5 * Temp * HiggsMassesVec[k + NHiggs] *
                    std::sqrt(HiggsMassesVec[k]);
        }
        // Explicit T part
        if (HiggsMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(HiggsMassesZeroTempVec[k], 1.5);
      }
      for (std::size_t k = 0; k < NGauge; k++)
      {
        if (GaugeMassesVec[k] > 0)
        {
          VDebye += std::pow(GaugeMassesVec[k], 1.5);
          VDebye += 1.5 * Temp * GaugeMassesVec[k + NGauge] *
                    std::sqrt(GaugeMassesVec[k]);
        }
        if (GaugeMassesZeroTempVec[k] > 0)
          VDebye += -std::pow(GaugeMassesZeroTempVec[k], 1.5);
      }

      VDebye *= -1.0 / (12 * M_PI);
      res += VDebye;
    }
  }

  return res;
}

void Class_Potential_Origin::CalculateDebye(bool forceCalculation)
{
  if (!SetCurvatureDone) SetCurvatureArrays();

  bool Calculate = forceCalculation or not CalculateDebyeSimplified();
  if (Calculate)
  {
    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = i; j < NHiggs; j++)
      {
        DebyeHiggs[i][j] = 0;
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          DebyeHiggs[i][j] += 0.5 * Curvature_Higgs_L4[i][j][k][k] / 12.0;
        }
        for (std::size_t k = 0; k < NGauge; k++)
        {
          DebyeHiggs[i][j] += 3 * 0.5 * Curvature_Gauge_G2H2[k][k][i][j] / 12.0;
        }

        for (std::size_t a = 0; a < NQuarks; a++)
        {
          for (std::size_t b = 0; b < NQuarks; b++)
          {
            double tmp = 0.5 * (std::conj(Curvature_Quark_F2H1[a][b][j]) *
                                    Curvature_Quark_F2H1[a][b][i] +
                                std::conj(Curvature_Quark_F2H1[a][b][i]) *
                                    Curvature_Quark_F2H1[a][b][j])
                                   .real();
            DebyeHiggs[i][j] += 6.0 / 24.0 * tmp;
          }
        }

        for (std::size_t a = 0; a < NLepton; a++)
        {
          for (std::size_t b = 0; b < NLepton; b++)
          {
            double tmp = 0.5 * (std::conj(Curvature_Lepton_F2H1[a][b][j]) *
                                    Curvature_Lepton_F2H1[a][b][i] +
                                std::conj(Curvature_Lepton_F2H1[a][b][i]) *
                                    Curvature_Lepton_F2H1[a][b][j])
                                   .real();
            DebyeHiggs[i][j] += 2.0 / 24.0 * tmp;
          }
        }

        //	            if(i==j) DebyeHiggs[i][j] *= 0.5;
      }
    }

    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = i; j < NHiggs; j++)
      {
        if (std::abs(DebyeHiggs[i][j]) <= 1e-5) DebyeHiggs[i][j] = 0;
      }
    }

    for (std::size_t i = 0; i < NHiggs; i++)
    {
      for (std::size_t j = 0; j < i; j++)
      {
        DebyeHiggs[i][j] = DebyeHiggs[j][i];
      }
    }
  }
}

void Class_Potential_Origin::CalculateDebyeGauge()
{
  for (std::size_t i = 0; i < NGauge; i++)
  {
    for (std::size_t j = 0; j < NGauge; j++)
      DebyeGauge[i][j] = 0;
  }

  bool Done = CalculateDebyeGaugeSimplified();
  if (Done) return;

  std::size_t nGaugeHiggs = 0;

  for (std::size_t i = 0; i < NHiggs; i++)
  {
    if (Curvature_Gauge_G2H2[0][0][i][i] != 0)
    {
      nGaugeHiggs++;
    }
  }
  for (std::size_t i = 0; i < NGauge; i++)
  {
    double GaugeFac = 0;
    for (std::size_t k = 0; k < NHiggs; k++)
    {
      GaugeFac += Curvature_Gauge_G2H2[i][i][k][k];
    }
    GaugeFac *= 1.0 / nGaugeHiggs;
    DebyeGauge[i][i] = 2.0 / 3.0 * (nGaugeHiggs / 8.0 + 5) * GaugeFac;
  }

  for (std::size_t i = 0; i < NGauge; i++)
  {
    for (std::size_t j = 0; j < NGauge; j++)
    {
      if (std::abs(DebyeGauge[i][j]) <= 1e-5) DebyeGauge[i][j] = 0;
    }
  }
}

void Class_Potential_Origin::initVectors()
{
  using vec2 = std::vector<std::vector<double>>;
  using vec3 = std::vector<std::vector<std::vector<double>>>;
  using vec4 = std::vector<std::vector<std::vector<std::vector<double>>>>;

  using vec1Complex = std::vector<std::complex<double>>;
  using vec2Complex = std::vector<std::vector<std::complex<double>>>;
  using vec3Complex =
      std::vector<std::vector<std::vector<std::complex<double>>>>;
  using vec4Complex =
      std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>>;

  Curvature_Higgs_L1 = std::vector<double>(NHiggs, 0);
  Curvature_Higgs_L2 = vec2{NHiggs, std::vector<double>(NHiggs, 0)};
  Curvature_Higgs_L3 =
      vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};
  Curvature_Higgs_L4 =
      vec4{NHiggs, vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}}};

  Curvature_Higgs_CT_L1 = std::vector<double>(NHiggs, 0);
  Curvature_Higgs_CT_L2 = vec2{NHiggs, std::vector<double>(NHiggs, 0)};
  Curvature_Higgs_CT_L3 =
      vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};
  Curvature_Higgs_CT_L4 =
      vec4{NHiggs, vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}}};

  VEVSymmetric = std::vector<double>(NHiggs, 0);

  DebyeHiggs = vec2{NHiggs, std::vector<double>(NHiggs, 0)};

  LambdaHiggs_3    = vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};
  LambdaHiggs_3_CT = vec3{NHiggs, vec2{NHiggs, std::vector<double>(NHiggs, 0)}};

  Curvature_Gauge_G2H2 =
      vec4{NGauge, vec3{NGauge, vec2{NHiggs, std::vector<double>(NHiggs, 0)}}};
  DebyeGauge    = vec2{NGauge, std::vector<double>(NGauge, 0)};
  LambdaGauge_3 = vec3{NGauge, vec2{NGauge, std::vector<double>(NHiggs, 0)}};

  Curvature_Lepton_F2 = vec2Complex{NLepton, vec1Complex(NLepton, 0)};
  Curvature_Lepton_F2H1 =
      vec3Complex{NLepton, vec2Complex{NLepton, vec1Complex(NHiggs, 0)}};
  LambdaLepton_3 =
      vec3Complex{NLepton, vec2Complex{NLepton, vec1Complex(NHiggs, 0)}};
  LambdaLepton_4 = vec4Complex{
      NLepton,
      vec3Complex{NLepton, vec2Complex{NHiggs, vec1Complex(NHiggs, 0)}}};

  Curvature_Quark_F2 = vec2Complex{NQuarks, vec1Complex(NQuarks, 0)};
  Curvature_Quark_F2H1 =
      vec3Complex{NQuarks, vec2Complex{NQuarks, vec1Complex(NHiggs, 0)}};
  LambdaQuark_3 =
      vec3Complex{NQuarks, vec2Complex{NQuarks, vec1Complex(NHiggs, 0)}};
  LambdaQuark_4 = vec4Complex{
      NQuarks,
      vec3Complex{NQuarks, vec2Complex{NHiggs, vec1Complex(NHiggs, 0)}}};

  HiggsVev = std::vector<double>(NHiggs, 0);

  HiggsRotationMatrixEnsuredConvention =
      std::vector<std::vector<double>>{NHiggs, std::vector<double>(NHiggs, 0)};

  HiggsTensorIndexCacheReady = false;
  HiggsTensorIndexCache.clear();
  GaugeTensorIndexCacheReady = false;
  GaugeTensorIndexCache.clear();
  QuarkTensorIndexCacheReady = false;
  QuarkMassFieldIndexCache.clear();
  QuarkDerivativeIndexCache.clear();
  LeptonTensorIndexCacheReady = false;
  LeptonMassFieldIndexCache.clear();
  LeptonDerivativeIndexCache.clear();
  VTreeTensorIndexCacheReady = false;
  VTreeTensorIndexCache.clear();
  CounterTermIndexCacheReady = false;
  CounterTermIndexCache.clear();
  CounterTermFlatIndexCacheReady = false;
  CounterTermFlatIndexCache.clear();
}

void Class_Potential_Origin::BuildHiggsTensorIndexCache()
{
  HiggsTensorIndexCache.clear();
  HiggsTensorIndexCache.resize(NHiggs * NHiggs);
  for (std::size_t i = 0; i < NHiggs; ++i)
  {
    for (std::size_t j = 0; j < NHiggs; ++j)
    {
      auto &terms = HiggsTensorIndexCache[i * NHiggs + j];
      terms.reserve(NHiggs);
      for (std::size_t k = 0; k < NHiggs; ++k)
      {
        HiggsTensorIndexTerms current;
        current.k          = k;
        current.l3Nonzero  = Curvature_Higgs_L3[i][j][k] != 0.0;
        current.l4Nonzero.reserve(NHiggs);
        for (std::size_t l = 0; l < NHiggs; ++l)
        {
          if (Curvature_Higgs_L4[i][j][k][l] != 0.0)
            current.l4Nonzero.push_back(l);
        }
        if (current.l3Nonzero || !current.l4Nonzero.empty())
          terms.push_back(std::move(current));
      }
    }
  }
  HiggsTensorIndexCacheReady = true;
}

void Class_Potential_Origin::BuildQuarkTensorIndexCache()
{
  QuarkMassFieldIndexCache.clear();
  QuarkMassFieldIndexCache.resize(NQuarks * NQuarks);
  for (std::size_t i = 0; i < NQuarks; ++i)
  {
    for (std::size_t j = 0; j < NQuarks; ++j)
    {
      auto &fields = QuarkMassFieldIndexCache[i * NQuarks + j];
      fields.reserve(NHiggs);
      for (std::size_t k = 0; k < NHiggs; ++k)
      {
        if (Curvature_Quark_F2H1[i][j][k] != std::complex<double>(0, 0))
          fields.push_back(k);
      }
    }
  }

  QuarkDerivativeIndexCache.clear();
  QuarkDerivativeIndexCache.resize(NHiggs * NQuarks * NQuarks);
  for (std::size_t m = 0; m < NHiggs; ++m)
  {
    for (std::size_t a = 0; a < NQuarks; ++a)
    {
      for (std::size_t b = 0; b < NQuarks; ++b)
      {
        auto &indices =
            QuarkDerivativeIndexCache[(m * NQuarks + a) * NQuarks + b];
        indices.reserve(NQuarks);
        for (std::size_t i = 0; i < NQuarks; ++i)
        {
          const bool firstNonzero =
              Curvature_Quark_F2H1[a][i][m] != std::complex<double>(0, 0);
          const bool secondNonzero =
              Curvature_Quark_F2H1[i][b][m] != std::complex<double>(0, 0);
          if (firstNonzero || secondNonzero)
            indices.push_back({i, firstNonzero, secondNonzero});
        }
      }
    }
  }
  QuarkTensorIndexCacheReady = true;
}

void Class_Potential_Origin::BuildLeptonTensorIndexCache()
{
  LeptonMassFieldIndexCache.clear();
  LeptonMassFieldIndexCache.resize(NLepton * NLepton);
  for (std::size_t i = 0; i < NLepton; ++i)
  {
    for (std::size_t j = 0; j < NLepton; ++j)
    {
      auto &fields = LeptonMassFieldIndexCache[i * NLepton + j];
      fields.reserve(NHiggs);
      for (std::size_t k = 0; k < NHiggs; ++k)
      {
        if (Curvature_Lepton_F2H1[i][j][k] != std::complex<double>(0, 0))
          fields.push_back(k);
      }
    }
  }

  LeptonDerivativeIndexCache.clear();
  LeptonDerivativeIndexCache.resize(NHiggs * NLepton * NLepton);
  for (std::size_t m = 0; m < NHiggs; ++m)
  {
    for (std::size_t a = 0; a < NLepton; ++a)
    {
      for (std::size_t b = 0; b < NLepton; ++b)
      {
        auto &indices =
            LeptonDerivativeIndexCache[(m * NLepton + a) * NLepton + b];
        indices.reserve(NLepton);
        for (std::size_t i = 0; i < NLepton; ++i)
        {
          const bool firstNonzero =
              Curvature_Lepton_F2H1[a][i][m] != std::complex<double>(0, 0);
          const bool secondNonzero =
              Curvature_Lepton_F2H1[i][b][m] != std::complex<double>(0, 0);
          if (firstNonzero || secondNonzero)
            indices.push_back({i, firstNonzero, secondNonzero});
        }
      }
    }
  }
  LeptonTensorIndexCacheReady = true;
}

void Class_Potential_Origin::BuildVTreeTensorIndexCache()
{
  VTreeTensorIndexCache.clear();
  VTreeTensorIndexCache.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; ++i)
  {
    auto &iTerms = VTreeTensorIndexCache[i];
    iTerms.l1Nonzero = Curvature_Higgs_L1[i] != 0.0;
    iTerms.jTerms.reserve(NHiggs);
    for (std::size_t j = 0; j < NHiggs; ++j)
    {
      VTreeJIndexTerms jTerm;
      jTerm.j          = j;
      jTerm.l2Nonzero  = Curvature_Higgs_L2[i][j] != 0.0;
      jTerm.kTerms.reserve(NHiggs);
      for (std::size_t k = 0; k < NHiggs; ++k)
      {
        VTreeKIndexTerms kTerm;
        kTerm.k          = k;
        kTerm.l3Nonzero  = Curvature_Higgs_L3[i][j][k] != 0.0;
        kTerm.l4Nonzero.reserve(NHiggs);
        for (std::size_t l = 0; l < NHiggs; ++l)
        {
          if (Curvature_Higgs_L4[i][j][k][l] != 0.0)
            kTerm.l4Nonzero.push_back(l);
        }
        if (kTerm.l3Nonzero || !kTerm.l4Nonzero.empty())
          jTerm.kTerms.push_back(std::move(kTerm));
      }
      if (jTerm.l2Nonzero || !jTerm.kTerms.empty())
        iTerms.jTerms.push_back(std::move(jTerm));
    }
  }
  VTreeTensorIndexCacheReady = true;
}

void Class_Potential_Origin::BuildGaugeTensorIndexCache()
{
  GaugeTensorIndexCache.clear();
  GaugeTensorIndexCache.resize(NGauge * NGauge);
  for (std::size_t a = 0; a < NGauge; ++a)
  {
    for (std::size_t b = 0; b < NGauge; ++b)
    {
      auto &terms = GaugeTensorIndexCache[a * NGauge + b];
      terms.mass.reserve(NHiggs * NHiggs);
      terms.derivativeByField.resize(NHiggs);
      for (std::size_t i = 0; i < NHiggs; ++i)
      {
        terms.derivativeByField[i].reserve(NHiggs);
        for (std::size_t j = 0; j < NHiggs; ++j)
        {
          if (Curvature_Gauge_G2H2[a][b][i][j] == 0.0) continue;
          terms.mass.emplace_back(i, j);
          terms.derivativeByField[i].push_back(j);
        }
      }
    }
  }
  GaugeTensorIndexCacheReady = true;
}

void Class_Potential_Origin::BuildCounterTermIndexCache()
{
  CounterTermIndexCache.clear();
  CounterTermIndexCache.resize(NHiggs);
  CounterTermFlatIndexCache.clear();
  CounterTermFlatIndexCache.resize(NHiggs);
  for (std::size_t i = 0; i < NHiggs; ++i)
  {
    auto &iTerms = CounterTermIndexCache[i];
    auto &flatTerms = CounterTermFlatIndexCache[i];
    iTerms.l1Nonzero = Curvature_Higgs_CT_L1[i] != 0.0;
    iTerms.jTerms.reserve(NHiggs);
    flatTerms.reserve(NHiggs * NHiggs * NHiggs);
    for (std::size_t j = 0; j < NHiggs; ++j)
    {
      CounterTermJIndexTerms jTerms;
      jTerms.j         = j;
      jTerms.l2Nonzero = Curvature_Higgs_CT_L2[i][j] != 0.0;
      if (jTerms.l2Nonzero)
        flatTerms.push_back({j, 0, 0, 2});
      jTerms.kTerms.reserve(NHiggs);
      for (std::size_t k = 0; k < NHiggs; ++k)
      {
        CounterTermKIndexTerms kTerms;
        kTerms.k         = k;
        kTerms.l3Nonzero = Curvature_Higgs_CT_L3[i][j][k] != 0.0;
        if (kTerms.l3Nonzero)
          flatTerms.push_back({j, k, 0, 3});
        kTerms.l4Nonzero.reserve(NHiggs);
        for (std::size_t l = 0; l < NHiggs; ++l)
        {
          if (Curvature_Higgs_CT_L4[i][j][k][l] != 0.0)
          {
            kTerms.l4Nonzero.push_back(l);
            flatTerms.push_back({j, k, l, 4});
          }
        }
        if (kTerms.l3Nonzero || !kTerms.l4Nonzero.empty())
          jTerms.kTerms.push_back(std::move(kTerms));
      }
      if (jTerms.l2Nonzero || !jTerms.kTerms.empty())
        iTerms.jTerms.push_back(std::move(jTerms));
    }
  }
  CounterTermIndexCacheReady = true;
  CounterTermFlatIndexCacheReady = true;
}

void Class_Potential_Origin::sym2Dim(
    std::vector<std::vector<double>> &Tensor2Dim,
    std::size_t Nk1,
    std::size_t Nk2)
{
  for (std::size_t k1 = 0; k1 < Nk1; k1++)
  {
    for (std::size_t k2 = k1; k2 < Nk2; k2++)
    {
      Tensor2Dim[k2][k1] = Tensor2Dim[k1][k2];
    }
  }
}

void Class_Potential_Origin::sym3Dim(
    std::vector<std::vector<std::vector<double>>> &Tensor3Dim,
    std::size_t Nk1,
    std::size_t Nk2,
    std::size_t Nk3)
{
  for (std::size_t k1 = 0; k1 < Nk1; k1++)
  {
    for (std::size_t k2 = k1; k2 < Nk2; k2++)
    {
      for (std::size_t k3 = k2; k3 < Nk3; k3++)
      {
        Tensor3Dim[k1][k3][k2] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k2][k1][k3] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k2][k3][k1] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k3][k1][k2] = Tensor3Dim[k1][k2][k3];
        Tensor3Dim[k3][k2][k1] = Tensor3Dim[k1][k2][k3];
      }
    }
  }
}

void Class_Potential_Origin::sym4Dim(
    std::vector<std::vector<std::vector<std::vector<double>>>> &Tensor4Dim,
    std::size_t Nk1,
    std::size_t Nk2,
    std::size_t Nk3,
    std::size_t Nk4)
{
  for (std::size_t k1 = 0; k1 < Nk1; k1++)
  {
    for (std::size_t k2 = k1; k2 < Nk2; k2++)
    {
      for (std::size_t k3 = k2; k3 < Nk3; k3++)
      {
        for (std::size_t k4 = k3; k4 < Nk4; k4++)
        {
          Tensor4Dim[k1][k2][k4][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k3][k2][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k3][k4][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k4][k2][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k1][k4][k3][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k1][k3][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k1][k4][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k3][k1][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k3][k4][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k4][k1][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k2][k4][k3][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k1][k2][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k1][k4][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k2][k1][k4] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k2][k4][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k4][k1][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k3][k4][k2][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k1][k2][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k1][k3][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k2][k1][k3] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k2][k3][k1] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k3][k1][k2] = Tensor4Dim[k1][k2][k3][k4];
          Tensor4Dim[k4][k3][k2][k1] = Tensor4Dim[k1][k2][k3][k4];
        }
      }
    }
  }
}

void Class_Potential_Origin::resetbools()
{
  SetCurvatureDone          = false;
  CalcCouplingsDone         = false;
  CalculatedTripleCopulings = false;
  parStored.clear();
  parCTStored.clear();
}

bool Class_Potential_Origin::CheckNLOVEV(const std::vector<double> &v) const
{
  // std::vector<double> vPotential;
  double MaxDiff           = 0;
  double AllowedDifference = 1;
  for (std::size_t i = 0; i < nVEV; i++)
  {
    double tmp = std::abs(std::abs(v[i]) - std::abs(vevTreeMin[i]));
    if (tmp > MaxDiff) MaxDiff = tmp;
  }

  return (MaxDiff < AllowedDifference);
}

double Class_Potential_Origin::EWSBVEV(const std::vector<double> &v) const
{
  double res = 0;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    double checkgauge = 0;
    for (std::size_t j = 0; j < NGauge; j++)
    {
      checkgauge += std::abs(Curvature_Gauge_G2H2[j][j][i][i]);
    }
    if (checkgauge != 0) res += std::pow(v.at(i), 2);
  }
  res = std::sqrt(res);

  if (res <= 0.5)
  {
    return 0;
  }
  return res;
}

void Class_Potential_Origin::SetEWVEVZero(std::vector<double> &sol) const
{
  int count = 0;
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    double checkgauge = 0;
    for (std::size_t j = 0; j < NGauge; j++)
    {
      checkgauge += std::abs(Curvature_Gauge_G2H2[j][j][i][i]);
    }
    if (checkgauge > 1e-10 and i == VevOrder[count]) // true for doublet vevs
    {
      sol.at(count) = 0;
      count += 1;
    }
  }
}

void Class_Potential_Origin::setUseIndexCol(std::string legend)
{
  UseIndexCol = legend.rfind(sep, 0) == 0;
}

bool Class_Potential_Origin::getUseIndexCol()
{
  return UseIndexCol;
}

void Class_Potential_Origin::FindSignSymmetries()
{
  SignSymmetries.clear();
  std::vector<double> testvev, testvevPotential;
  for (std::size_t i = 0; i < nVEV; i++)
    testvev.push_back(i + 1);
  testvevPotential      = MinimizeOrderVEV(testvev);
  double referenceValue = VEff(testvevPotential);

  std::vector<double> vevdummy, vevdummyPotential;
  double VEffDummy;

  // Fill a dummy vector with a certain amout of -1 and look at all possible
  // permutations of it
  for (std::size_t countNegative = 1; countNegative <= nVEV; countNegative++)
  {
    std::vector<double> tmpSymmetry;
    for (std::size_t i = 0; i < countNegative; i++)
      tmpSymmetry.push_back(-1);
    for (std::size_t i = countNegative; i < nVEV; i++)
      tmpSymmetry.push_back(1);

    do
    {
      vevdummy.clear();
      for (std::size_t i = 0; i < nVEV; i++)
        vevdummy.push_back(tmpSymmetry.at(i) * testvev.at(i));
      vevdummyPotential = MinimizeOrderVEV(vevdummy);
      VEffDummy         = VEff(vevdummyPotential);
      if (std::abs(VEffDummy - referenceValue) <=
          1e-3 * std::abs(referenceValue))
        SignSymmetries.push_back(tmpSymmetry);
    } while (std::next_permutation(tmpSymmetry.begin(), tmpSymmetry.end()));
  }
}

void Class_Potential_Origin::SetUseTreeLevel(bool val)
{
  UseTreeLevel = val;
}

std::pair<std::vector<double>, std::vector<double>>
Class_Potential_Origin::initModel(std::string linestr)
{
  std::vector<double> par(nPar), parCT(nParCT);
  resetbools();
  ReadAndSet(linestr, par);
  parCT = initModel(par);

  parStored   = par;
  parCTStored = parCT;

  std::pair<std::vector<double>, std::vector<double>> res;
  res.first  = par;
  res.second = parCT;

  return res;
}

std::vector<double>
Class_Potential_Origin::initModel(const std::vector<double> &par,
                                  const bool &adjust_rot_matrix)
{
  std::vector<double> parCT(nParCT);
  resetbools();
  set_gen(par);
  CalculatePhysicalCouplings();
  parCT = calc_CT();
  set_CT_Pot_Par(parCT);
  CalculateDebye();
  CalculateDebyeGauge();

  if (adjust_rot_matrix)
  {
    AdjustRotationMatrix();
  }

  parStored   = par;
  parCTStored = parCT;

  return parCT;
}

std::vector<double> Class_Potential_Origin::resetScale(const double &newScale)
{
  scale      = newScale;
  auto parCT = calc_CT();
  set_CT_Pot_Par(parCT);

  parCTStored = parCT;

  return parCT;
}

Eigen::MatrixXcd
Class_Potential_Origin::QuarkMassMatrix(const std::vector<double> &v) const
{
  MatrixXcd MIJ(NQuarks, NQuarks);
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    MIJ          = QuarkMassMatrix(Transformedv);
    return MIJ;
  }
  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " is called before SetCurvatureArrays() is called. \n";
    throw std::runtime_error(retmes);
  }

  MIJ = MatrixXcd::Zero(NQuarks, NQuarks);

  static const bool UseR2HDMQuarkIndexCache = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_QUARK_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool useIndexCache =
      UseR2HDMQuarkIndexCache && Model == ModelID::ModelIDs::R2HDM &&
      NQuarks == 12 && QuarkTensorIndexCacheReady &&
      QuarkMassFieldIndexCache.size() == NQuarks * NQuarks;

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    for (std::size_t j = 0; j < NQuarks; j++)
    {
      MIJ(i, j) = Curvature_Quark_F2[i][j];
      if (useIndexCache)
      {
        const auto &fields = QuarkMassFieldIndexCache[i * NQuarks + j];
        for (const std::size_t k : fields)
          MIJ(i, j) += Curvature_Quark_F2H1[i][j][k] * v[k];
      }
      else
      {
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          MIJ(i, j) += Curvature_Quark_F2H1[i][j][k] * v[k];
        }
      }
    }
  }

  return MIJ;
}

std::vector<std::complex<double>>
Class_Potential_Origin::QuarkMasses(const std::vector<double> &v) const
{
  std::vector<std::complex<double>> res;
  double ZeroMass = std::pow(10, -10);

  auto MIJ = QuarkMassMatrix(v);

  ComplexEigenSolver<MatrixXcd> es(MIJ, false);

  for (std::size_t i = 0; i < NQuarks; i++)
  {
    auto tmp = es.eigenvalues()[i];
    if (std::abs(tmp) < ZeroMass)
      res.push_back(0);
    else
      res.push_back(tmp);
  }

  return res;
}

MatrixXcd
Class_Potential_Origin::LeptonMassMatrix(const std::vector<double> &v) const
{
  MatrixXcd res = MatrixXcd::Zero(NLepton, NLepton);
  if (v.size() != nVEV and v.size() != NHiggs)
  {
    std::string ErrorString =
        std::string("You have called ") + std::string(__func__) +
        std::string(
            " with an invalid vev configuration. Your vev is of dimension ") +
        std::to_string(v.size()) + std::string(" and it should be ") +
        std::to_string(NHiggs) + std::string(".");
    throw std::runtime_error(ErrorString);
  }
  if (v.size() == nVEV and nVEV != NHiggs)
  {
    std::stringstream ss;
    ss << __func__
       << " is being called with a wrong sized vev configuration. It "
          "has the dimension of "
       << nVEV << " while it should have " << NHiggs
       << ". For now this is transformed but please fix this to reduce "
          "the runtime."
       << std::endl;
    Logger::Write(LoggingLevel::Default, ss.str());
    std::vector<double> Transformedv;
    Transformedv = MinimizeOrderVEV(v);
    res          = LeptonMassMatrix(Transformedv);
    return res;
  }
  if (!SetCurvatureDone)
  {
    std::string retmes = __func__;
    retmes += " is called before SetCurvatureArrays();\n";
    throw std::runtime_error(retmes);
  }

  static const bool UseR2HDMLeptonIndexCache = []
  {
    const char *env = std::getenv("BSMPT_USE_R2HDM_LEPTON_INDEX_CACHE");
    return env != nullptr && env[0] == '1';
  }();
  const bool useIndexCache =
      UseR2HDMLeptonIndexCache && Model == ModelID::ModelIDs::R2HDM &&
      NLepton == 9 && LeptonTensorIndexCacheReady &&
      LeptonMassFieldIndexCache.size() == NLepton * NLepton;

  for (std::size_t i = 0; i < NLepton; i++)
  {
    for (std::size_t j = 0; j < NLepton; j++)
    {
      res(i, j) = Curvature_Lepton_F2[i][j];
      if (useIndexCache)
      {
        const auto &fields = LeptonMassFieldIndexCache[i * NLepton + j];
        for (const std::size_t k : fields)
          res(i, j) += Curvature_Lepton_F2H1[i][j][k] * v[k];
      }
      else
      {
        for (std::size_t k = 0; k < NHiggs; k++)
        {
          res(i, j) += Curvature_Lepton_F2H1[i][j][k] * v[k];
        }
      }
    }
  }

  return res;
}

std::vector<std::complex<double>>
Class_Potential_Origin::LeptonMasses(const std::vector<double> &v) const
{
  std::vector<std::complex<double>> res;

  MatrixXcd MIJ = LeptonMassMatrix(v);

  ComplexEigenSolver<MatrixXcd> es(MIJ, false);
  double ZeroMass = 1e-10;
  for (std::size_t i = 0; i < NLepton; i++)
  {
    auto tmp = es.eigenvalues()[i];
    if (std::abs(tmp) < ZeroMass)
      res.push_back(0);
    else
      res.push_back(tmp);
  }
  return res;
}

double Class_Potential_Origin::CalculateRatioAlpha(
    const std::vector<double> &vev_symmetric,
    const std::vector<double> &vev_broken,
    const double &Temp) const
{
  (void)vev_symmetric;
  (void)vev_broken;
  (void)Temp;
  //  double res                          = 0;
  //  double PotentialSymmetricPhaseValue = VEff(vev_symmetric, Temp, 0);
  //  double PotentialSymmetricPhaseDeriv = VEff(vev_symmetric, Temp, -1);
  //  double PotentialBrokenPhaseValue    = VEff(vev_broken, Temp, 0);
  //  double PotentialBrokenPhaseDeriv    = VEff(vev_broken, Temp, -1);
  //  res = -(PotentialBrokenPhaseValue - PotentialSymmetricPhaseValue) +
  //        Temp * (PotentialBrokenPhaseDeriv - PotentialSymmetricPhaseDeriv);
  // TODO:: Unfinished!
  throw std::runtime_error("The CalculateRatioAlpha function is not finished "
                           "and you should not be using it!");
}

std::vector<double> Class_Potential_Origin::MinimizeOrderVEV(
    const std::vector<double> &vevMinimizer) const
{
  std::vector<double> vevFunction;

  std::size_t count = 0;
  for (std::size_t i = 0; i < NHiggs; ++i)
  {
    if (i == VevOrder[count])
    {
      vevFunction.push_back(vevMinimizer.at(count));
      count++;
    }
    else
      vevFunction.push_back(0);
  }
  return vevFunction;
}

Eigen::VectorXd
Class_Potential_Origin::NablaVCT(const std::vector<double> &v) const
{
  VectorXd result(NHiggs);
  for (std::size_t i{0}; i < NHiggs; ++i)
  {
    result(i) = Curvature_Higgs_CT_L1[i];
    for (std::size_t j{0}; j < NHiggs; ++j)
    {
      result(i) += Curvature_Higgs_CT_L2[i][j] * v.at(j);
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        result(i) += 0.5 * Curvature_Higgs_CT_L3[i][j][k] * v.at(j) * v.at(k);
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          result(i) += 1.0 / 6.0 * Curvature_Higgs_CT_L4[i][j][k][l] * v.at(j) *
                       v.at(k) * v.at(l);
        }
      }
    }
  }
  return result;
}

Eigen::MatrixXd
Class_Potential_Origin::HessianCT(const std::vector<double> &v) const
{
  Eigen::MatrixXd result(NHiggs, NHiggs);
  for (std::size_t i = 0; i < NHiggs; i++)
  {
    for (std::size_t j = 0; j < NHiggs; j++)
    {
      result(i, j) = Curvature_Higgs_CT_L2[i][j];
      for (std::size_t k = 0; k < NHiggs; k++)
      {
        result(i, j) += Curvature_Higgs_CT_L3[i][j][k] * v.at(k);
        for (std::size_t l = 0; l < NHiggs; l++)
        {
          result(i, j) +=
              0.5 * Curvature_Higgs_CT_L4[i][j][k][l] * v.at(k) * v.at(l);
        }
      }
    }
  }
  return result;
}

std::vector<double> Class_Potential_Origin::GetCTIdentities() const
{
  return std::vector<double>();
}

} // namespace BSMPT

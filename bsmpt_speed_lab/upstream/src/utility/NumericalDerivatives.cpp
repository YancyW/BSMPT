#include <BSMPT/utility/NumericalDerivatives.h>
#include <array>
#include <cstdlib>

namespace BSMPT
{

namespace
{
// The public result type intentionally remains vector-based because all
// existing callers consume std::vector/std::vector<std::vector>.  This helper
// only moves the two axial scratch arrays to the stack for the common
// four-field R2HDM case.  The callback still accepts std::vector<double>, so
// the sample coordinate itself cannot be changed to std::array without adding
// a conversion/allocation at every V call.
std::pair<std::vector<double>, std::vector<std::vector<double>>>
NablaAndHessianNumericalFixed4(
    const std::vector<double> &phi,
    const std::function<double(const std::vector<double> &)> &V,
    double eps,
    double gradient_offset,
    const std::vector<double> *active_hessian_direction)
{
  constexpr std::size_t dimension = 4;
  std::vector<double> gradient(dimension);
  std::vector<std::vector<double>> hessian(
      dimension, std::vector<double>(dimension));

  const double value_at_phi = V(phi);
  std::vector<double> xp(phi);
  std::array<double, dimension> fp2{};
  std::array<double, dimension> fm2{};

  for (std::size_t i = 0; i < dimension; ++i)
  {
    xp[i] = phi[i] + eps;
    const double fp1 = V(xp) - gradient_offset;
    xp[i] = phi[i] - eps;
    const double fm1 = V(xp) - gradient_offset;
    xp[i] = phi[i] + 2 * eps;
    fp2[i] = V(xp);
    const double fp2_gradient = fp2[i] - gradient_offset;
    xp[i] = phi[i] - 2 * eps;
    fm2[i] = V(xp);
    const double fm2_gradient = fm2[i] - gradient_offset;
    xp[i] = phi[i];

    gradient[i] =
        (-fp2_gradient + 8 * fp1 - 8 * fm1 + fm2_gradient) / (12 * eps);
  }

  for (std::size_t i = 0; i < dimension; ++i)
  {
    if (active_hessian_direction != nullptr &&
        active_hessian_direction->size() == dimension &&
        (*active_hessian_direction)[i] == 0.0)
      continue;

    double val = 0;
    val += fp2[i];
    val -= 2 * value_at_phi;
    val += fm2[i];
    hessian[i][i] = val / (4 * eps * eps);

    for (std::size_t j = i + 1; j < dimension; ++j)
    {
      if (active_hessian_direction != nullptr &&
          active_hessian_direction->size() == dimension &&
          (*active_hessian_direction)[j] == 0.0)
        continue;

      double r = 0;
      xp[i] = phi[i] + eps;
      xp[j] = phi[j] + eps;
      r += V(xp);

      xp[j] = phi[j] - eps;
      r -= V(xp);

      xp[i] = phi[i] - eps;
      xp[j] = phi[j] + eps;
      r -= V(xp);

      xp[j] = phi[j] - eps;
      r += V(xp);
      xp[i] = phi[i];
      xp[j] = phi[j];

      hessian[i][j] = r / (4 * eps * eps);
      hessian[j][i] = r / (4 * eps * eps);
    }
  }

  return {std::move(gradient), std::move(hessian)};
}
} // namespace

std::vector<double>
NablaNumerical(const std::vector<double> &phi,
               const std::function<double(std::vector<double>)> &f,
               const double &eps)
{
  std::vector<double> result(phi.size());
  std::vector<double> shifted(phi);
  static const bool use_central_two_point =
      std::getenv("BSMPT_USE_CENTRAL2_GRADIENT") != nullptr;

  for (size_t i = 0; i < phi.size(); i++)
  {
    shifted[i]       = phi[i] + eps;
    const double fp1 = f(shifted);
    shifted[i]       = phi[i] - eps;
    const double fm1 = f(shifted);
    shifted[i]       = phi[i];
    if (use_central_two_point)
    {
      result[i] = (fp1 - fm1) / (2 * eps);
    }
    else
    {
      shifted[i]       = phi[i] + 2 * eps;
      const double fp2 = f(shifted);
      shifted[i]       = phi[i] - 2 * eps;
      const double fm2 = f(shifted);
      shifted[i]       = phi[i];
      result[i] = (-fp2 + 8 * fp1 - 8 * fm1 + fm2) / (12 * eps);
    }
  }
  return result;
}

std::vector<std::vector<double>>
HessianNumerical(const std::vector<double> &phi,
                 const std::function<double(std::vector<double>)> &V,
                 double eps)
{
  std::vector<std::vector<double>> result(phi.size(),
                                          std::vector<double>(phi.size()));
  const double value_at_phi = V(phi);
  std::vector<double> xp(phi);
  for (size_t i = 0; i < phi.size(); i++)
  {
    double val = 0;
    xp[i] = phi[i] + 2 * eps;
    val += V(xp);

    val -= 2 * value_at_phi;

    xp[i] = phi[i] - 2 * eps;
    val += V(xp);
    xp[i] = phi[i];

    result[i][i] = val / (4 * eps * eps);

    // https://en.wikipedia.org/wiki/Finite_difference
    for (size_t j = i + 1; j < phi.size(); j++)
    {
      double r = 0;

      xp[i] = phi[i] + eps; // F(x+h, y+h)
      xp[j] = phi[j] + eps;
      r += V(xp);

      xp[j] = phi[j] - eps; //-F(x+h, y-h)
      r -= V(xp);

      xp[i] = phi[i] - eps; //-F(x-h, y+h)
      xp[j] = phi[j] + eps;
      r -= V(xp);

      xp[j] = phi[j] - eps; // F(x-h, y-h)
      r += V(xp);
      xp[i] = phi[i];
      xp[j] = phi[j];

      result[i][j] = r / (4 * eps * eps);
      result[j][i] = r / (4 * eps * eps);
    }
  }

  return result;
}

std::pair<std::vector<double>, std::vector<std::vector<double>>>
NablaAndHessianNumerical(
    const std::vector<double> &phi,
    const std::function<double(const std::vector<double> &)> &V,
    double eps,
    double gradient_offset,
    const std::vector<double> *active_hessian_direction)
{
  const std::size_t dimension = phi.size();
  static const bool use_fixed_r2hdm4 = []
  {
    const char *env =
        std::getenv("BSMPT_USE_R2HDM_FIXED_COMBINED_DERIVATIVES");
    return env != nullptr && env[0] == '1';
  }();
  if (use_fixed_r2hdm4 && dimension == 4)
    return NablaAndHessianNumericalFixed4(phi,
                                          V,
                                          eps,
                                          gradient_offset,
                                          active_hessian_direction);

  std::vector<double> gradient(dimension);
  std::vector<std::vector<double>> hessian(
      dimension, std::vector<double>(dimension));

  // This is intentionally the unshifted value: HessianNumerical uses V(phi)
  // directly, while gradient_offset reproduces a constant shift in the
  // gradient's potential without changing the Hessian stencil.
  const double value_at_phi = V(phi);
  std::vector<double> xp(phi);
  std::vector<double> fp2(dimension);
  std::vector<double> fm2(dimension);

  // First collect all axial values.  The +/- eps evaluations are used only by
  // the four-point gradient; +/-2 eps are shared with the Hessian diagonal.
  for (std::size_t i = 0; i < dimension; ++i)
  {
    xp[i] = phi[i] + eps;
    const double fp1 = V(xp) - gradient_offset;
    xp[i] = phi[i] - eps;
    const double fm1 = V(xp) - gradient_offset;
    xp[i] = phi[i] + 2 * eps;
    fp2[i] = V(xp);
    const double fp2_gradient = fp2[i] - gradient_offset;
    xp[i] = phi[i] - 2 * eps;
    fm2[i] = V(xp);
    const double fm2_gradient = fm2[i] - gradient_offset;
    xp[i] = phi[i];

    gradient[i] =
        (-fp2_gradient + 8 * fp1 - 8 * fm1 + fm2_gradient) / (12 * eps);
  }

  // Keep the diagonal and mixed-coordinate accumulation order of
  // HessianNumerical exactly; only the already-computed axial values are
  // substituted for their repeated function evaluations.
  for (std::size_t i = 0; i < dimension; ++i)
  {
    if (active_hessian_direction != nullptr &&
        active_hessian_direction->size() == dimension &&
        (*active_hessian_direction)[i] == 0.0)
      continue;

    double val = 0;
    val += fp2[i];
    val -= 2 * value_at_phi;
    val += fm2[i];
    hessian[i][i] = val / (4 * eps * eps);

    for (std::size_t j = i + 1; j < dimension; ++j)
    {
      if (active_hessian_direction != nullptr &&
          active_hessian_direction->size() == dimension &&
          (*active_hessian_direction)[j] == 0.0)
        continue;

      double r = 0;

      xp[i] = phi[i] + eps; // F(x+h, y+h)
      xp[j] = phi[j] + eps;
      r += V(xp);

      xp[j] = phi[j] - eps; //-F(x+h, y-h)
      r -= V(xp);

      xp[i] = phi[i] - eps; //-F(x-h, y+h)
      xp[j] = phi[j] + eps;
      r -= V(xp);

      xp[j] = phi[j] - eps; // F(x-h, y-h)
      r += V(xp);
      xp[i] = phi[i];
      xp[j] = phi[j];

      hessian[i][j] = r / (4 * eps * eps);
      hessian[j][i] = r / (4 * eps * eps);
    }
  }

  return {std::move(gradient), std::move(hessian)};
}
} // namespace BSMPT

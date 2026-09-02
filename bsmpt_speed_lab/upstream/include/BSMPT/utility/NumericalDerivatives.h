#pragma once
#include <functional>
#include <utility>
#include <vector>

namespace BSMPT
{
/**
 * @brief Numerical method to calculate the
 * gradient of a function f using finite differences method.
 *
 * This method is used while BSMPT is not able to
 * calculate the potential derivative analytically. We used the 4th order
 * method
 *
 * \f$\frac{\partial f}{\partial \phi_i} = \frac{1}{12
 * \epsilon}\left(-f(\dots ,\vec{\phi}_i + 2  \epsilon ) + 8 f(\dots
 * ,\vec{\phi}_i + \epsilon )- 8 f(\dots ,\vec{\phi}_i - \epsilon ) +
 * f(\dots ,\vec{\phi}_i - 2  \epsilon )\right)\f$
 *
 * where \f$ \epsilon \f$ is a small step.
 *
 * @param phi Where we want to calculate the gradient
 * @param f function
 * @param eps Size of finite differences step
 * @return std::vector<double> The \f$ dim \times 1 \f$ gradient of V taken at
 * phi
 */
std::vector<double>
NablaNumerical(const std::vector<double> &phi,
               const std::function<double(std::vector<double>)> &f,
               const double &eps);

/**
 * @brief Numerical method to calculate the potential's (or other functions's)
 * hessian matrix using finite differences method.
 *
 * \f$\frac{\partial^2 V}{\partial \phi_i \phi_j} = \frac{1}{4
 * \epsilon^2}\left(V(\dots, \vec{\phi}_i + \epsilon , \vec{\phi}_j +
 * \epsilon) - V(\dots, \vec{\phi}_i - \epsilon , \vec{\phi}_j +
 * \epsilon) - V(\dots, \vec{\phi}_i + \epsilon , \vec{\phi}_j -
 * \epsilon) + V(\dots, \vec{\phi}_i - \epsilon , \vec{\phi}_j -
 * \epsilon) \right)\f$
 *
 * where \f$ \epsilon \f$ is a small step.
 *
 * @param phi Where we want to calculate the Hessian matrix
 * @param V Potential (or other function)
 * @param eps Size of finite differences step
 * @return std::vector<std::vector<double>> The \f$ dim \times \dim \f$
 *  hessian matrix of V taken at phi
 */
std::vector<std::vector<double>>
HessianNumerical(const std::vector<double> &phi,
                 const std::function<double(std::vector<double>)> &V,
                 double eps);

/**
 * @brief Calculate a numerical gradient and Hessian while sharing the center
 * and axial +/-2 eps evaluations.
 *
 * The gradient uses the same four-point stencil as NablaNumerical and the
 * Hessian uses the same diagonal and mixed-coordinate stencils as
 * HessianNumerical.  If gradient_offset is nonzero, it is subtracted from
 * each value used by the gradient only; this preserves the common BSMPT case
 * where the gradient uses a potential shifted by a constant while the
 * Hessian uses the unshifted potential.
 * If active_hessian_direction is provided with matching dimension, Hessian
 * rows/columns whose corresponding direction entry is exactly zero are left
 * zero; the gradient remains full-dimensional.
 */
std::pair<std::vector<double>, std::vector<std::vector<double>>>
NablaAndHessianNumerical(
    const std::vector<double> &phi,
    const std::function<double(const std::vector<double> &)> &V,
    double eps,
    double gradient_offset = 0.0,
    const std::vector<double> *active_hessian_direction = nullptr);
} // namespace BSMPT

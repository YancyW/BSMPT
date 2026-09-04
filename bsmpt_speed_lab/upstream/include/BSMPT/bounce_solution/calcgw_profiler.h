// SPDX-FileCopyrightText: 2026
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <cstdint>
#include <chrono>

namespace BSMPT
{

/**
 * @brief Very small, process-wide counters for CalcGW's bounce solver.
 *
 * The counters are disabled unless BSMPT_CALCGW_PROFILE is set to a non-zero
 * value.  They are deliberately kept separate from the model code so that a
 * profiling run cannot change any numerical path through the calculation.
 */
class CalcGWProfiler final
{
public:
  enum class TimingMetric : std::uint8_t
  {
    VTree,
    CounterTerm,
    HiggsMasses,
    GaugeMasses,
    QuarkMasses,
    LeptonMasses,
    BosonThermal,
    FermionThermal,
    Count
  };

  static bool enabled() noexcept;

  /**
   * @brief Whether low-overhead VEff/V1Loop timing is enabled.
   *
   * This is intentionally separate from the existing bounce counters.  It is
   * enabled only when BSMPT_CALCGW_PROFILE_VEFF is set to a non-zero value.
   */
  static bool timing_enabled() noexcept;
  static bool thermal_repeat_enabled() noexcept;
  static void thermal_repeat_call(bool fermion, int diff, bool hit) noexcept;


  static void record_timing(TimingMetric metric,
                            std::uint64_t nanoseconds) noexcept;

  /**
   * @brief Aggregate timer for one VEff/V1Loop hot-path call.
   *
   * The constructor performs no clock read when timing is disabled.  This
   * keeps the normal CalcGW path free of timing-system calls.
   */
  class ScopedTimer final
  {
  public:
    explicit ScopedTimer(TimingMetric metric) noexcept
        : metric_{metric}, active_{CalcGWProfiler::timing_enabled()}
    {
      if (active_) start_ = clock::now();
    }

    ScopedTimer(TimingMetric metric, bool active) noexcept
        : metric_{metric}, active_{active}
    {
      if (active_) start_ = clock::now();
    }

    ~ScopedTimer() noexcept
    {
      if (active_)
      {
        const auto elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(clock::now() -
                                                                  start_)
                .count();
        CalcGWProfiler::record_timing(
            metric_, static_cast<std::uint64_t>(elapsed));
      }
    }

    ScopedTimer(const ScopedTimer &)            = delete;
    ScopedTimer &operator=(const ScopedTimer &) = delete;

  private:
    using clock = std::chrono::steady_clock;
    TimingMetric metric_;
    bool active_ = false;
    clock::time_point start_{};
  };

  static void rasterized_call(std::uint64_t samples) noexcept;
  static void exact_threshold_call() noexcept;
  static void exact_threshold_iteration() noexcept;
  static void integrate_call() noexcept;
  static void rk5_step() noexcept;
  static void solve_1d_call() noexcept;
  static void kinetic_action_samples(std::uint64_t samples) noexcept;
  static void potential_action_samples(std::uint64_t samples) noexcept;
  static void vtree_call() noexcept;
  static void vtree_terms(std::uint64_t terms) noexcept;
  static void counterterm_call() noexcept;
  static void counterterm_terms(std::uint64_t terms) noexcept;
  static void active_hessian_dimensions(std::uint64_t dimensions) noexcept;

  /// Print one aggregate line at process exit when profiling is enabled.
  static void report() noexcept;
};

} // namespace BSMPT

// SPDX-License-Identifier: GPL-3.0-or-later

#include <BSMPT/bounce_solution/calcgw_profiler.h>

#include <atomic>
#include <array>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace BSMPT
{
namespace
{
struct Counters
{
  std::atomic<std::uint64_t> rasterized_calls{0};
  std::atomic<std::uint64_t> rasterized_samples{0};
  std::atomic<std::uint64_t> exact_threshold_calls{0};
  std::atomic<std::uint64_t> exact_threshold_iterations{0};
  std::atomic<std::uint64_t> integrate_calls{0};
  std::atomic<std::uint64_t> rk5_steps{0};
  std::atomic<std::uint64_t> solve_1d_calls{0};
  std::atomic<std::uint64_t> kinetic_action_samples{0};
  std::atomic<std::uint64_t> potential_action_samples{0};
  std::atomic<std::uint64_t> vtree_calls{0};
  std::atomic<std::uint64_t> vtree_terms{0};
  std::atomic<std::uint64_t> counterterm_calls{0};
  std::atomic<std::uint64_t> counterterm_terms{0};
  std::atomic<std::uint64_t> active_hessian_dimensions{0};
  std::array<std::atomic<std::uint64_t>, static_cast<std::size_t>(
                                             CalcGWProfiler::TimingMetric::Count)>
      timing_calls{};
  std::array<std::atomic<std::uint64_t>, static_cast<std::size_t>(
                                             CalcGWProfiler::TimingMetric::Count)>
      timing_nanoseconds{};
  std::array<std::atomic<std::uint64_t>, 8> thermal_repeat_calls{};
  std::array<std::atomic<std::uint64_t>, 8> thermal_repeat_hits{};
  std::array<std::atomic<std::uint64_t>, 3> exact_repeat_calls{};
  std::array<std::atomic<std::uint64_t>, 3> exact_repeat_hits{};

  Counters() noexcept
  {
    for (auto &value : timing_calls) value.store(0, std::memory_order_relaxed);
    for (auto &value : timing_nanoseconds)
      value.store(0, std::memory_order_relaxed);
    for (auto &value : thermal_repeat_calls)
      value.store(0, std::memory_order_relaxed);
    for (auto &value : thermal_repeat_hits)
      value.store(0, std::memory_order_relaxed);
    for (auto &value : exact_repeat_calls)
      value.store(0, std::memory_order_relaxed);
    for (auto &value : exact_repeat_hits)
      value.store(0, std::memory_order_relaxed);
  }
};

Counters &counters() noexcept
{
  static Counters value;
  return value;
}

bool read_enabled() noexcept
{
  const char *value = std::getenv("BSMPT_CALCGW_PROFILE");
  return value != nullptr && *value != '\0' && std::strcmp(value, "0") != 0;
}

bool read_timing_enabled() noexcept
{
  const char *value = std::getenv("BSMPT_CALCGW_PROFILE_VEFF");
  return value != nullptr && *value != '\0' && std::strcmp(value, "0") != 0;
}

bool read_thermal_repeat_enabled() noexcept
{
  const char *value = std::getenv("BSMPT_PROFILE_THERMAL_REPEATS");
  return value != nullptr && *value != '\0' && std::strcmp(value, "0") != 0;
}

bool read_exact_repeat_enabled() noexcept
{
  const char *value = std::getenv("BSMPT_PROFILE_VEFF_REPEATS");
  return value != nullptr && *value != '\0' && std::strcmp(value, "0") != 0;
}

void ensure_report_registration() noexcept
{
  static const bool registered = [] {
    if (read_enabled() || read_timing_enabled() ||
        read_thermal_repeat_enabled() || read_exact_repeat_enabled())
      std::atexit(&CalcGWProfiler::report);
    return true;
  }();
  (void)registered;
}
} // namespace

bool CalcGWProfiler::enabled() noexcept
{
  static const bool value = [] {
    const bool on = read_enabled();
    if (on) ensure_report_registration();
    return on;
  }();
  return value;
}

bool CalcGWProfiler::timing_enabled() noexcept
{
  static const bool value = [] {
    const bool on = read_timing_enabled();
    if (on) ensure_report_registration();
    return on;
  }();
  return value;
}

bool CalcGWProfiler::thermal_repeat_enabled() noexcept
{
  static const bool value = [] {
    const bool on = read_thermal_repeat_enabled();
    if (on) ensure_report_registration();
    return on;
  }();
  return value;
}

bool CalcGWProfiler::exact_repeat_enabled() noexcept
{
  static const bool value = [] {
    const bool on = read_exact_repeat_enabled();
    if (on) ensure_report_registration();
    return on;
  }();
  return value;
}

void CalcGWProfiler::exact_repeat_call(ExactRepeatMetric metric,
                                       bool hit) noexcept
{
  if (!exact_repeat_enabled()) return;
  const auto index = static_cast<std::size_t>(metric);
  if (index >= static_cast<std::size_t>(ExactRepeatMetric::Count)) return;
  counters().exact_repeat_calls[index].fetch_add(1,
                                                 std::memory_order_relaxed);
  if (hit)
    counters().exact_repeat_hits[index].fetch_add(1,
                                                  std::memory_order_relaxed);
}

void CalcGWProfiler::thermal_repeat_call(bool fermion,
                                         int diff,
                                         bool hit) noexcept
{
  if (!thermal_repeat_enabled()) return;
  const std::size_t diff_index =
      diff == 0 ? 0 : (diff == 1 ? 1 : (diff == -1 ? 2 : 3));
  const std::size_t index = (fermion ? 4 : 0) + diff_index;
  counters().thermal_repeat_calls[index].fetch_add(1,
                                                    std::memory_order_relaxed);
  if (hit)
    counters().thermal_repeat_hits[index].fetch_add(1,
                                                    std::memory_order_relaxed);
}

void CalcGWProfiler::record_timing(TimingMetric metric,
                                   std::uint64_t nanoseconds) noexcept
{
  if (!timing_enabled()) return;
  const auto index = static_cast<std::size_t>(metric);
  if (index >= static_cast<std::size_t>(TimingMetric::Count)) return;
  counters().timing_calls[index].fetch_add(1, std::memory_order_relaxed);
  counters().timing_nanoseconds[index].fetch_add(nanoseconds,
                                                 std::memory_order_relaxed);
}

void CalcGWProfiler::rasterized_call(std::uint64_t samples) noexcept
{
  if (!enabled()) return;
  counters().rasterized_calls.fetch_add(1, std::memory_order_relaxed);
  counters().rasterized_samples.fetch_add(samples, std::memory_order_relaxed);
}

void CalcGWProfiler::exact_threshold_call() noexcept
{
  if (!enabled()) return;
  counters().exact_threshold_calls.fetch_add(1, std::memory_order_relaxed);
}

void CalcGWProfiler::exact_threshold_iteration() noexcept
{
  if (!enabled()) return;
  counters().exact_threshold_iterations.fetch_add(1,
                                                  std::memory_order_relaxed);
}

void CalcGWProfiler::integrate_call() noexcept
{
  if (!enabled()) return;
  counters().integrate_calls.fetch_add(1, std::memory_order_relaxed);
}

void CalcGWProfiler::rk5_step() noexcept
{
  if (!enabled()) return;
  counters().rk5_steps.fetch_add(1, std::memory_order_relaxed);
}

void CalcGWProfiler::solve_1d_call() noexcept
{
  if (!enabled()) return;
  counters().solve_1d_calls.fetch_add(1, std::memory_order_relaxed);
}

void CalcGWProfiler::potential_action_samples(std::uint64_t samples) noexcept
{
  if (!enabled()) return;
  counters().potential_action_samples.fetch_add(samples,
                                                std::memory_order_relaxed);
}

void CalcGWProfiler::kinetic_action_samples(std::uint64_t samples) noexcept
{
  if (!enabled()) return;
  counters().kinetic_action_samples.fetch_add(samples,
                                              std::memory_order_relaxed);
}

void CalcGWProfiler::vtree_call() noexcept
{
  if (!enabled()) return;
  counters().vtree_calls.fetch_add(1, std::memory_order_relaxed);
}

void CalcGWProfiler::vtree_terms(std::uint64_t terms) noexcept
{
  if (!enabled()) return;
  counters().vtree_terms.fetch_add(terms, std::memory_order_relaxed);
}

void CalcGWProfiler::counterterm_call() noexcept
{
  if (!enabled()) return;
  counters().counterterm_calls.fetch_add(1, std::memory_order_relaxed);
}

void CalcGWProfiler::counterterm_terms(std::uint64_t terms) noexcept
{
  if (!enabled()) return;
  counters().counterterm_terms.fetch_add(terms, std::memory_order_relaxed);
}

void CalcGWProfiler::active_hessian_dimensions(
    std::uint64_t dimensions) noexcept
{
  if (!enabled()) return;
  counters().active_hessian_dimensions.fetch_add(dimensions,
                                                 std::memory_order_relaxed);
}

void CalcGWProfiler::report() noexcept
{
  if (!enabled() && !timing_enabled() && !thermal_repeat_enabled() &&
      !exact_repeat_enabled())
    return;
  const auto &c = counters();
  std::cerr << "BSMPT_CALCGW_PROFILE"
            << " rasterized_calls=" << c.rasterized_calls.load()
            << " rasterized_samples=" << c.rasterized_samples.load()
            << " threshold_calls=" << c.exact_threshold_calls.load()
            << " threshold_iterations=" << c.exact_threshold_iterations.load()
            << " integrate_calls=" << c.integrate_calls.load()
            << " rk5_steps=" << c.rk5_steps.load()
            << " solve1d_calls=" << c.solve_1d_calls.load()
            << " kinetic_action_samples="
            << c.kinetic_action_samples.load()
            << " potential_action_samples="
            << c.potential_action_samples.load()
            << " vtree_calls=" << c.vtree_calls.load()
            << " vtree_terms=" << c.vtree_terms.load()
            << " counterterm_calls=" << c.counterterm_calls.load()
            << " counterterm_terms=" << c.counterterm_terms.load()
            << " active_hessian_dimensions="
            << c.active_hessian_dimensions.load();
  if (timing_enabled())
  {
    static constexpr const char *labels[] = {
        "VTree",       "CounterTerm", "HiggsMasses", "GaugeMasses",
        "QuarkMasses", "LeptonMasses", "BosonThermal", "FermionThermal"};
    std::cerr << " veff_timing_ns=";
    for (std::size_t i = 0;
         i < static_cast<std::size_t>(TimingMetric::Count); ++i)
    {
      if (i != 0) std::cerr << ',';
      std::cerr << labels[i] << ':' << c.timing_nanoseconds[i].load();
    }
    std::cerr << " veff_calls=";
    for (std::size_t i = 0;
         i < static_cast<std::size_t>(TimingMetric::Count); ++i)
    {
      if (i != 0) std::cerr << ',';
      std::cerr << labels[i] << ':' << c.timing_calls[i].load();
    }
  }
  if (thermal_repeat_enabled())
  {
    static constexpr const char *labels[] = {
        "boson_d0", "boson_d1", "boson_dm1", "boson_other",
        "fermion_d0", "fermion_d1", "fermion_dm1", "fermion_other"};
    std::cerr << " thermal_repeat=";
    for (std::size_t i = 0; i < 8; ++i)
    {
      if (i != 0) std::cerr << ',';
      std::cerr << labels[i] << ':' << c.thermal_repeat_hits[i].load() << '/'
                << c.thermal_repeat_calls[i].load();
    }
  }
  if (exact_repeat_enabled())
  {
    static constexpr const char *labels[] = {"veff", "higgs", "quark"};
    std::cerr << " exact_repeat=";
    for (std::size_t i = 0; i < 3; ++i)
    {
      if (i != 0) std::cerr << ',';
      std::cerr << labels[i] << ':' << c.exact_repeat_hits[i].load() << '/'
                << c.exact_repeat_calls[i].load();
    }
  }
  std::cerr << '\n';
}

} // namespace BSMPT

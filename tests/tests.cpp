#include "fmi.h"
#include "signal.h"
#include "results.h"
#include "examples.h"

#include <gtest/gtest.h>

#include <cassert>
#include <cmath>
#include <vector>
#include <memory>

using namespace fmi;

static bool are_approx_equal(double expected, double actual, double tol) {
  return fabs(expected - actual) < tol;
}

TEST(Tests, test_shift_to_beginning_0) {
  const Interval interval{ 3, 5 };
  std::vector<fmi2Real> values({ 7 });
  const Derivatives derivatives{ values, interval.end };
  const Sample sample{ interval, derivatives };
  const Sample shifted = fmi::signal::shift_to_beginning(sample);
  assert(are_approx_equal(7, shifted.derivatives.values[0], 1e-6));
  assert(are_approx_equal(3, shifted.derivatives.t, 1e-6));
  assert(shifted.derivatives.values.size() == 1);
}

TEST(Tests, test_shift_to_beginning_1) {
  const Interval interval{ 3, 5 };
  std::vector<fmi2Real> values({ 7, -2 });
  const Derivatives derivatives{ values, interval.end };
  const Sample sample{ interval, derivatives };
  const Sample shifted = fmi::signal::shift_to_beginning(sample);
  assert(are_approx_equal(11, shifted.derivatives.values[0], 1e-6));
  assert(are_approx_equal(-2, shifted.derivatives.values[1], 1e-6));
  assert(are_approx_equal(3, shifted.derivatives.t, 1e-6));
  assert(shifted.derivatives.values.size() == 2);
}

TEST(Tests, test_shift_to_beginning_2) {
  const Interval interval{ 3, 5 };
  std::vector<fmi2Real> values({ 7, -2, 5 });
  const Derivatives derivatives{ values, interval.end };
  const Sample sample{ interval, derivatives };
  const Sample shifted = fmi::signal::shift_to_beginning(sample);
  assert(are_approx_equal(5, shifted.derivatives.values[2], 1e-6));
  assert(are_approx_equal(-12, shifted.derivatives.values[1], 1e-6));
  assert(are_approx_equal(21, shifted.derivatives.values[0], 1e-6));
  assert(are_approx_equal(3, shifted.derivatives.t, 1e-6));
  assert(shifted.derivatives.values.size() == 3);
}

TEST(Tests, test_tau2omega_hermite) {
  double h = 3.;

  CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2);
  std::shared_ptr<FMU> tau2omega = network->instances["tau2omega"];

  fmi::InputDerivatives u = {
    std::vector<fmi2Real>({7., 11., 13.}),
    0.
  };
  tau2omega->set_input(0, u);
  tau2omega->step(h);
  fmi::OutputDerivatives y = tau2omega->get_output(OMEGA_T2O_VR);
  fmi::OutputDerivatives yh = tau2omega->get_hermite(OMEGA_T2O_VR, OMEGA_T2O_H_VR, h);

  assert(y.values.size() == 3);
  assert(yh.values.size() == 4);

  for (size_t i = 0; i < 3; i++) {
    assert(are_approx_equal(y.values[i], yh.values[i], 1e-6));
  }

  assert(are_approx_equal(fmi::signal::eval(y, h), fmi::signal::eval(yh, h), 1e-6));
}

TEST(Tests, test_eval) {
  double t0 = 0., tEnd = 0.91;
  fmi::Interval interval = { t0, tEnd };
  fmi::Derivatives derivatives = {
    std::vector<fmi2Real>({5., 3., 4.}),
    t0
  };
  fmi::Sample sample = { interval, derivatives };

  // polyval([2, 3, 5], 0.7)
  assert(are_approx_equal(fmi::signal::eval(sample.derivatives, 0.7), 8.08, 1e-6));

  sample.derivatives.t = tEnd;
  // polyval([2, 3, 5], -0.21)
  assert(are_approx_equal(fmi::signal::eval(sample.derivatives, 0.7), 4.4582, 1e-6));
}

TEST(Tests, test_variable_step_twomass) {
  auto network = fmi::examples::instantiate_twomass(2, 2);
  auto instance = network->instances["twomass"];
  double h = 0.13;
  instance->step(h);
  Results results;
  record_outputs(network, h, results);
  record_output_defects(results, network->monitoredOutputs);
  fmi::Signal signal = { fmi::Defect, {"twomass", OMEGA_2M_VR} };
  auto sample = results[signal].back();
  assert(fmi::signal::rms(sample) > 0.);
}

template <typename StepGenerator>
static void test_helper(StepGenerator stepGenerator, size_t extrapolationOrder) {
  double tEnd = 3;
  fmi::Results resultsGenerator = fmi::examples::analytic_experiment(stepGenerator, fmi::Both, extrapolationOrder);
  fmi::Signal signal0 = { fmi::SignalType::Values, { "twomass", OMEGA_2M_VR } };
  fmi::Signal signal1 = { fmi::SignalType::Values, { "twomass", TAU_2M_VR } };
  fmi::StepSequence stepSequence(resultsGenerator[signal0]);
  fmi::Results resultsSequence = fmi::examples::analytic_experiment(stepSequence, fmi::Both, extrapolationOrder);
  assert(are_approx_equal(0., fmi::results::calculate_error(resultsGenerator[signal0], resultsSequence[signal0]), 1e-6));
  assert(are_approx_equal(0., fmi::results::calculate_error(resultsGenerator[signal1], resultsSequence[signal1]), 1e-6));
}

TEST(Tests, test_step_sequence) {
  for (size_t extrapolationOrder = 0; extrapolationOrder <= 2; extrapolationOrder++)
  {
    double h0 = 0.7 * (extrapolationOrder + 1);
    test_helper(ConstantStepSize(h0, 3), extrapolationOrder);
    test_helper(VariableStepSize(h0, 3, 0.1), extrapolationOrder);
  }
}

TEST(Tests, test_square_1) {
  Derivatives sample = { std::vector<fmi2Real>({5.}), 3. };
  Derivatives squared = fmi::signal::square(sample);
  assert(are_approx_equal(3., squared.t, 1e-6));
  assert(squared.values.size() == 1);
  assert(are_approx_equal(25., squared.values[0], 1e-6));
}

TEST(Tests, test_square_2) {
  Derivatives sample = { std::vector<fmi2Real>({7., 11.}), 5. };
  Derivatives squared = fmi::signal::square(sample);
  assert(are_approx_equal(5., squared.t, 1e-6));
  assert(squared.values.size() == 3);
  assert(are_approx_equal(49., squared.values[0], 1e-6));
  assert(are_approx_equal(154., squared.values[1], 1e-6));
  assert(are_approx_equal(121. * 2., squared.values[2], 1e-6));
}

TEST(Tests, test_square_3) {
  Derivatives sample = { std::vector<fmi2Real>({2., 3., 5.}), 7. };
  Derivatives squared = fmi::signal::square(sample);
  assert(squared.values.size() == 5);
  assert(are_approx_equal(4., squared.values[0], 1e-6));
  assert(are_approx_equal(12., squared.values[1], 1e-6));
  assert(are_approx_equal(19. * 2., squared.values[2], 1e-6));
  assert(are_approx_equal(15. * 6., squared.values[3], 1e-6));
  assert(are_approx_equal(6.25 * 24., squared.values[4], 1e-6));
}

TEST(Tests, test_derivative_difference) {
  Derivatives sample1 = { std::vector<fmi2Real>({2., 3.}), 0. };
  Derivatives sample2 = { std::vector<fmi2Real>({5.}), 0. };

  Derivatives diff1 = fmi::signal::calculate_difference(sample1, sample2);
  assert(diff1.values.size() == 2);
  assert(are_approx_equal(-3., diff1.values[0], 1e-6));
  assert(are_approx_equal(3., diff1.values[1], 1e-6));

  Derivatives diff2 = fmi::signal::calculate_difference(sample2, sample1);
  assert(diff2.values.size() == 2);
  assert(are_approx_equal(3., diff2.values[0], 1e-6));
  assert(are_approx_equal(-3., diff2.values[1], 1e-6));
}

TEST(Tests, test_sample_rms_1) {
  Interval interval = { 0., 1. };
  Derivatives der1 = { std::vector<fmi2Real>({2., 3.}), 0. };
  Derivatives der2 = { std::vector<fmi2Real>({5.}), 0. };

  Sample sample1 = {interval, der1 }, sample2 = {interval, der2 };

  Samples samples1({sample1}), samples2({ sample2 });

  Samples diffSamples = fmi::signal::calculate_difference(samples1, samples2);
  Sample& diffSample = diffSamples[0];

  assert(diffSample.derivatives.values.size() == 2);
  assert(are_approx_equal(-3., diffSample.derivatives.values[0], 1e-6));
  assert(are_approx_equal(3., diffSample.derivatives.values[1], 1e-6));

  double rmsdiffSample = fmi::signal::rms(diffSample);
  double rmsdiffSamples = fmi::signal::rms(diffSamples);

  assert(are_approx_equal(sqrt(3.), rmsdiffSample, 1e-6));
  assert(are_approx_equal(sqrt(3.), rmsdiffSamples, 1e-6));
}

TEST(Tests, test_sample_rms_2) {
    Interval interval01 = { 0., 1. }, interval12 = { 1., 2. };
    Derivatives commonDer = { std::vector<fmi2Real>({2., 3.}), 0. };

    Derivatives der1 = { std::vector<fmi2Real>({5.}), 0. };
    Derivatives der2 = { std::vector<fmi2Real>({2., 3.}), 0. };

    Sample commonSample = { interval01, commonDer };
    Sample sample1 = { interval12, der1 }, sample2 = { interval12, der2 };

    Samples samples1({ commonSample, sample1 }), samples2({ commonSample, sample2 });

    Samples diffSamples = fmi::signal::calculate_difference(samples1, samples2);
    Sample& diffSample = diffSamples[1];

    double rmsdiffSample = fmi::signal::rms(diffSample);
    double rmsdiffSamples = fmi::signal::rms(diffSamples);

    assert(are_approx_equal(sqrt(3.), rmsdiffSample, 1e-6));
    assert(are_approx_equal(sqrt(3. / 2.), rmsdiffSamples, 1e-6));
}

TEST(Tests, test_t2o_0_euler) {
 CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  
  CosimulationNetworkPtr analyticNetwork = examples::instantiate_twomass(2);

  std::shared_ptr<FMU> t2o = network->instances["tau2omega"];

  double J_Omega2Tau = 10;
  double c_Omega2Tau = 1;
  double d_Omega2Tau = 1;
  double phi0_Omega2Tau = 0.1;
  double omega0_Omega2Tau = 0.1;
  double J_Tau2Omega = 10;
  double c_Tau2Omega = 1;
  double d_Tau2Omega = 2;
  double phi0_Tau2Omega = 0.2;
  double omega0_Tau2Omega = 0.1;
  double ck = 1;
  double dk = 2;

  OutputDerivatives t2oOmega = t2o->get_output(OMEGA_T2O_VR);

  double expectedTau = ck * phi0_Omega2Tau + dk * omega0_Omega2Tau - ck * phi0_Tau2Omega - dk * omega0_Tau2Omega;
  double expectedT2OAlpha = (-c_Tau2Omega * phi0_Tau2Omega - d_Tau2Omega * omega0_Tau2Omega + expectedTau) / J_Tau2Omega;
  double expectedT2OOmega[3] = {
    omega0_Tau2Omega,
    0.,
    0.
  };
  assert(are_approx_equal(t2oOmega.t, 0, 1e-6));
  for (size_t i = 0; i < 3; i++)
  {
    assert(are_approx_equal(t2oOmega.values[i], expectedT2OOmega[i], 1e-6));
  }
}

TEST(Tests, test_t2o_h_euler) {
 CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  
  CosimulationNetworkPtr analyticNetwork = examples::instantiate_twomass(2);

  std::shared_ptr<FMU> t2o = network->instances["tau2omega"];
  std::shared_ptr<FMU> twomass = analyticNetwork->instances["twomass"];

  double J_Omega2Tau = 10;
  double c_Omega2Tau = 1;
  double d_Omega2Tau = 1;
  double phi0_Omega2Tau = 0.1;
  double omega0_Omega2Tau = 0.1;
  double J_Tau2Omega = 10;
  double c_Tau2Omega = 1;
  double d_Tau2Omega = 2;
  double phi0_Tau2Omega = 0.2;
  double omega0_Tau2Omega = 0.1;
  double ck = 1;
  double dk = 2;

  OutputDerivatives twomassTau = twomass->get_output(TAU_2M_VR);

  double h = 0.5;
  t2o->step(h);
  OutputDerivatives t2oOmega = t2o->get_output(OMEGA_T2O_VR);

  double expectedDOmega = (-c_Tau2Omega * phi0_Tau2Omega - d_Tau2Omega * omega0_Tau2Omega + twomassTau.values[0]) / J_Tau2Omega;
  double expectedOmega = omega0_Omega2Tau + h * expectedDOmega;

  double expectedT2OOmega[3] = {
    expectedOmega,
    expectedDOmega,
    0.
  };
  assert(are_approx_equal(t2oOmega.t, h, 1e-6));
  for (size_t i = 0; i <= 2; i++)
  {
    assert(are_approx_equal(t2oOmega.values[i], expectedT2OOmega[i], 1e-6));
  }
}

TEST(Tests, test_o2t_0_euler) {
 CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  
  CosimulationNetworkPtr analyticNetwork = examples::instantiate_twomass(2);

  std::shared_ptr<FMU> o2t = network->instances["omega2tau"];
  std::shared_ptr<FMU> twomass = analyticNetwork->instances["twomass"];

  double J_Omega2Tau = 10;
  double c_Omega2Tau = 1;
  double d_Omega2Tau = 1;
  double phi0_Omega2Tau = 0.1;
  double omega0_Omega2Tau = 0.1;
  double J_Tau2Omega = 10;
  double c_Tau2Omega = 1;
  double d_Tau2Omega = 2;
  double phi0_Tau2Omega = 0.2;
  double omega0_Tau2Omega = 0.1;
  double ck = 1;
  double dk = 2;

  OutputDerivatives twomassOmega = twomass->get_output(OMEGA_2M_VR);

  OutputDerivatives o2tTau = o2t->get_output(TAU_O2T_VR);

  double expectedTau = ck * phi0_Omega2Tau + dk * omega0_Omega2Tau - ck * phi0_Tau2Omega - dk * twomassOmega.values[0];
  double expectedO2TTau[3] = {
    expectedTau,
    - dk * twomassOmega.values[1],
    - dk * twomassOmega.values[2]
  };
  assert(are_approx_equal(o2tTau.t, 0, 1e-6));
  for (size_t i = 0; i < 3; i++)
  {
    assert(are_approx_equal(o2tTau.values[i], expectedO2TTau[i], 1e-6));
  }
}

TEST(Tests, test_o2t_h_euler) {
 CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  
  CosimulationNetworkPtr analyticNetwork = examples::instantiate_twomass(2);

  std::shared_ptr<FMU> o2t = network->instances["omega2tau"];
  std::shared_ptr<FMU> twomass = analyticNetwork->instances["twomass"];

  double J_Omega2Tau = 10;
  double c_Omega2Tau = 1;
  double d_Omega2Tau = 1;
  double phi0_Omega2Tau = 0.1;
  double omega0_Omega2Tau = 0.1;
  double J_Tau2Omega = 10;
  double c_Tau2Omega = 1;
  double d_Tau2Omega = 2;
  double phi0_Tau2Omega = 0.2;
  double omega0_Tau2Omega = 0.1;
  double ck = 1;
  double dk = 2;

  OutputDerivatives twomassOmega = twomass->get_output(OMEGA_2M_VR);

  double h = 0.5;
  o2t->step(h);

  OutputDerivatives o2tTau = o2t->get_output(TAU_O2T_VR);

  double expectedDPhiO2T = omega0_Omega2Tau;
  double expectedPhiO2T[3] = {
    phi0_Omega2Tau + h * expectedDPhiO2T,
    expectedDPhiO2T,
    0.
  };

  double expectedDOmegaO2T =
    -(c_Omega2Tau + ck) / J_Omega2Tau * phi0_Omega2Tau
    - (d_Omega2Tau + dk) / J_Omega2Tau * omega0_Omega2Tau
    + ck / J_Omega2Tau * phi0_Tau2Omega
    + dk / J_Omega2Tau * twomassOmega.values[0];
  double expectedOmegaO2T[3] = {
    omega0_Omega2Tau + h * expectedDOmegaO2T,
    expectedDOmegaO2T,
    0.
  };

  double expectedDPhiT2O = twomassOmega.values[0];
  double expectedPhiT2O[3] = {
    phi0_Tau2Omega + h * expectedDPhiT2O,
    expectedDPhiT2O,
    0.
  };

  std::vector<fmi2Real> shiftedOmegaT2O = signal::shift_derivatives(twomassOmega.values, h);

  assert(are_approx_equal(o2tTau.t, h, 1e-6));
  for (size_t i = 0; i < 3; i++)
  {
    double expectedO2tTau = ck * expectedPhiO2T[i] + dk * expectedOmegaO2T[i] - ck * expectedPhiT2O[i] - dk * shiftedOmegaT2O[i];
    assert(are_approx_equal(o2tTau.values[i], expectedO2tTau, 1e-6));
  }
}

TEST(Tests, test_t2o_0_cvode) {
 CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Cvode);
  
  CosimulationNetworkPtr analyticNetwork = examples::instantiate_twomass(2);

  std::shared_ptr<FMU> t2o = network->instances["tau2omega"];

  double J_Omega2Tau = 10;
  double c_Omega2Tau = 1;
  double d_Omega2Tau = 1;
  double phi0_Omega2Tau = 0.1;
  double omega0_Omega2Tau = 0.1;
  double J_Tau2Omega = 10;
  double c_Tau2Omega = 1;
  double d_Tau2Omega = 2;
  double phi0_Tau2Omega = 0.2;
  double omega0_Tau2Omega = 0.1;
  double ck = 1;
  double dk = 2;

  OutputDerivatives t2oOmega = t2o->get_output(OMEGA_T2O_VR);

  double expectedTau = ck * phi0_Omega2Tau + dk * omega0_Omega2Tau - ck * phi0_Tau2Omega - dk * omega0_Tau2Omega;
  double expectedT2OAlpha = (-c_Tau2Omega * phi0_Tau2Omega - d_Tau2Omega * omega0_Tau2Omega + expectedTau) / J_Tau2Omega;
  double expectedT2OOmega[3] = {
    omega0_Tau2Omega,
    0.,
    0.
  };
  assert(are_approx_equal(t2oOmega.t, 0, 1e-6));
  for (size_t i = 0; i < 3; i++)
  {
    assert(are_approx_equal(t2oOmega.values[i], expectedT2OOmega[i], 1e-6));
  }
}

TEST(Tests, test_o2t_0_cvode) {
 CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Cvode);
  
  CosimulationNetworkPtr analyticNetwork = examples::instantiate_twomass(2);

  std::shared_ptr<FMU> o2t = network->instances["omega2tau"];
  std::shared_ptr<FMU> twomass = analyticNetwork->instances["twomass"];

  double J_Omega2Tau = 10;
  double c_Omega2Tau = 1;
  double d_Omega2Tau = 1;
  double phi0_Omega2Tau = 0.1;
  double omega0_Omega2Tau = 0.1;
  double J_Tau2Omega = 10;
  double c_Tau2Omega = 1;
  double d_Tau2Omega = 2;
  double phi0_Tau2Omega = 0.2;
  double omega0_Tau2Omega = 0.1;
  double ck = 1;
  double dk = 2;

  OutputDerivatives twomassOmega = twomass->get_output(OMEGA_2M_VR);

  OutputDerivatives o2tTau = o2t->get_output(TAU_O2T_VR);

  double expectedTau = ck * phi0_Omega2Tau + dk * omega0_Omega2Tau - ck * phi0_Tau2Omega - dk * twomassOmega.values[0];
  double expectedO2TTau[3] = {
    expectedTau,
    - dk * twomassOmega.values[1],
    - dk * twomassOmega.values[2]
  };
  assert(are_approx_equal(o2tTau.t, 0, 1e-6));
  for (size_t i = 0; i < 3; i++)
  {
    assert(are_approx_equal(o2tTau.values[i], expectedO2TTau[i], 1e-6));
  }
}

TEST(Tests, test_oscillator_3) {
  CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2);
  CosimulationNetworkPtr analyticNetwork = examples::instantiate_twomass(2);

  std::shared_ptr<FMU> t2o = network->instances["tau2omega"];
  std::shared_ptr<FMU> o2t = network->instances["omega2tau"];
  std::shared_ptr<FMU> twomass = analyticNetwork->instances["twomass"];

  double h = 1e-3;

  OutputDerivatives twomassOmega = twomass->get_output(OMEGA_2M_VR);
  OutputDerivatives twomassTau = twomass->get_output(TAU_2M_VR);

  t2o->set_input(0, twomassTau);
  o2t->set_input(0, twomassOmega);

  OutputDerivatives t2oOmega = t2o->get_output(OMEGA_T2O_VR);
  OutputDerivatives o2tTau = o2t->get_output(TAU_O2T_VR);
  twomassOmega = twomass->get_output(OMEGA_2M_VR);
  twomassTau = twomass->get_output(TAU_2M_VR);

  assert(are_approx_equal(t2oOmega.t, twomassOmega.t, 1e-6));
  assert(are_approx_equal(o2tTau.t, twomassTau.t, 1e-6));

  t2o->step(h);
  o2t->step(h);
  twomass->step(h);

  t2oOmega = t2o->get_output(OMEGA_T2O_VR);
  o2tTau = o2t->get_output(TAU_O2T_VR);
  twomassOmega = twomass->get_output(OMEGA_2M_VR);
  twomassTau = twomass->get_output(TAU_2M_VR);

  assert(are_approx_equal(t2oOmega.t, twomassOmega.t, 1e-6));
  assert(are_approx_equal(o2tTau.t, twomassTau.t, 1e-6));
  for (size_t i = 0; i <= 2; i++)
  {
      assert(are_approx_equal(t2oOmega.values[i], twomassOmega.values[i], 1e-3));
      assert(are_approx_equal(o2tTau.values[i], twomassTau.values[i], 1e-3));
  }

  t2o->step(h);
  o2t->step(h);
  twomass->step(h);

  t2oOmega = t2o->get_output(OMEGA_T2O_VR);
  o2tTau = o2t->get_output(TAU_O2T_VR);
  twomassOmega = twomass->get_output(OMEGA_2M_VR);
  twomassTau = twomass->get_output(TAU_2M_VR);

  assert(are_approx_equal(t2oOmega.t, twomassOmega.t, 1e-6));
  assert(are_approx_equal(o2tTau.t, twomassTau.t, 1e-6));
  for (size_t i = 0; i <= 2; i++)
  {
      assert(are_approx_equal(t2oOmega.values[i], twomassOmega.values[i], 1e-3));
      assert(are_approx_equal(o2tTau.values[i], twomassTau.values[i], 1e-3));
  }
}

TEST(Tests, test_t2o_euler_intermediate_out_samples) {
  double h = 1.;
  InputDerivatives in = { std::vector<fmi2Real>({ 2., 3., 5. }), 0. };

  CosimulationNetworkPtr networkHalf = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  std::shared_ptr<FMU> t2oHalf = networkHalf->instances["tau2omega"];
  t2oHalf->set_input(TAU_T2O_VR, in);
  t2oHalf->step(h / 2);
  OutputDerivatives t2oOmegaHalf = t2oHalf->get_output(OMEGA_T2O_VR);

  CosimulationNetworkPtr networkFull = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  std::shared_ptr<FMU> t2oFull = networkFull->instances["tau2omega"];
  t2oFull->set_input(TAU_T2O_VR, in);
  t2oFull->step(h);
  OutputDerivatives t2oOmegaFull = t2oFull->get_output(OMEGA_T2O_VR);
  OutputDerivatives t2oOmegaHermite = t2oFull->get_hermite(OMEGA_T2O_VR, OMEGA_T2O_H_VR, h);
  fmi2Real hermiteEval = fmi::signal::eval(t2oOmegaHermite, h / 2);

  assert(are_approx_equal(t2oOmegaHalf.values[0], hermiteEval, 1e-6));
  assert(are_approx_equal(t2oOmegaFull.values[0], t2oOmegaHermite.values[0], 1e-6));
  assert(are_approx_equal(t2oOmegaFull.values[1], t2oOmegaHermite.values[1], 1e-6));
  assert(are_approx_equal(t2oOmegaFull.values[2], t2oOmegaHermite.values[2], 1e-6));
}

TEST(Tests, test_t2o_cvode_intermediate_out_samples) {
  double h = 1.;
  InputDerivatives in = { std::vector<fmi2Real>({ 2., 3., 5. }), 0. };

  CosimulationNetworkPtr networkHalf = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Cvode);
  std::shared_ptr<FMU> t2oHalf = networkHalf->instances["tau2omega"];
  t2oHalf->set_input(TAU_T2O_VR, in);
  t2oHalf->step(h / 2);
  OutputDerivatives t2oOmegaHalf = t2oHalf->get_output(OMEGA_T2O_VR);

  CosimulationNetworkPtr networkFull = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Cvode);
  std::shared_ptr<FMU> t2oFull = networkFull->instances["tau2omega"];
  t2oFull->set_input(TAU_T2O_VR, in);
  t2oFull->step(h);
  OutputDerivatives t2oOmegaFull = t2oFull->get_output(OMEGA_T2O_VR);
  OutputDerivatives t2oOmegaHermite = t2oFull->get_hermite(OMEGA_T2O_VR, OMEGA_T2O_H_VR, h);
  fmi2Real hermiteEval = fmi::signal::eval(t2oOmegaHermite, h / 2);

  assert(are_approx_equal(t2oOmegaHalf.values[0], hermiteEval, 1e-6));
  assert(are_approx_equal(t2oOmegaFull.values[0], t2oOmegaHermite.values[0], 1e-6));
  assert(are_approx_equal(t2oOmegaFull.values[1], t2oOmegaHermite.values[1], 1e-6));
  assert(are_approx_equal(t2oOmegaFull.values[2], t2oOmegaHermite.values[2], 1e-6));
}

TEST(Tests, test_o2t_euler_intermediate_out_samples) {
  double h = 1.;
  InputDerivatives in = { std::vector<fmi2Real>({ 2., 3., 5. }), 0. };

  CosimulationNetworkPtr networkHalf = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  std::shared_ptr<FMU> o2tHalf = networkHalf->instances["omega2tau"];
  o2tHalf->set_input(OMEGA_T2O_VR, in);
  o2tHalf->step(h / 2);
  OutputDerivatives o2tTauHalf = o2tHalf->get_output(TAU_O2T_VR);

  CosimulationNetworkPtr networkFull = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);
  std::shared_ptr<FMU> o2tFull = networkFull->instances["omega2tau"];
  o2tFull->set_input(OMEGA_T2O_VR, in);
  o2tFull->step(h);
  OutputDerivatives o2tTauFull = o2tFull->get_output(TAU_O2T_VR);
  OutputDerivatives o2tTauHermite = o2tFull->get_hermite(TAU_O2T_VR, TAU_O2T_H_VR, h);
  fmi2Real hermiteEval = fmi::signal::eval(o2tTauHermite, h / 2);

  assert(are_approx_equal(o2tTauHalf.values[0], hermiteEval, 1e-6));
  assert(are_approx_equal(o2tTauFull.values[0], o2tTauHermite.values[0], 1e-6));
  assert(are_approx_equal(o2tTauFull.values[1], o2tTauHermite.values[1], 1e-6));
  assert(are_approx_equal(o2tTauFull.values[2], o2tTauHermite.values[2], 1e-6));
}

TEST(Tests, test_o2t_cvode_intermediate_out_samples) {
  double h = 1.;
  InputDerivatives in = { std::vector<fmi2Real>({ 2., 3., 5. }), 0. };

  CosimulationNetworkPtr networkHalf = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Cvode);
  std::shared_ptr<FMU> o2tHalf = networkHalf->instances["omega2tau"];
  o2tHalf->set_input(OMEGA_T2O_VR, in);
  o2tHalf->step(h / 2);
  OutputDerivatives o2tTauHalf = o2tHalf->get_output(TAU_O2T_VR);

  CosimulationNetworkPtr networkFull = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Cvode);
  std::shared_ptr<FMU> o2tFull = networkFull->instances["omega2tau"];
  o2tFull->set_input(OMEGA_T2O_VR, in);
  o2tFull->step(h);
  OutputDerivatives o2tTauFull = o2tFull->get_output(TAU_O2T_VR);
  OutputDerivatives o2tTauHermite = o2tFull->get_hermite(TAU_O2T_VR, TAU_O2T_H_VR, h);
  fmi2Real hermiteEval = fmi::signal::eval(o2tTauHermite, h / 2);

  assert(are_approx_equal(o2tTauHalf.values[0], hermiteEval, 1e-6));
  assert(are_approx_equal(o2tTauFull.values[0], o2tTauHermite.values[0], 1e-6));
  assert(are_approx_equal(o2tTauFull.values[1], o2tTauHermite.values[1], 1e-6));
  assert(are_approx_equal(o2tTauFull.values[2], o2tTauHermite.values[2], 1e-6));
}

TEST(Tests, test_euler_defects) {
  CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Euler);

  Results results = jacobi_co_simulation(ConstantStepSize(0.25, 10.), Both, network);

  Signal out{ Values, { "tau2omega", OMEGA_T2O_VR} };
  Signal outH{ Values, { "tau2omega", OMEGA_T2O_H_VR} };
  Signal outDefect{ Defect, { "tau2omega", OMEGA_T2O_VR} };

  for (size_t i = 0; i < results[out].size(); i++) {
    Sample outSample = results[out][i];
    Sample outHSample = results[outH][i];
    Sample outDefectSample = results[outDefect][i];
    assert(are_approx_equal(outSample.derivatives.values[0], outHSample.derivatives.values[0], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[1], outHSample.derivatives.values[1], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[2], outHSample.derivatives.values[2], 1e-6));
  }

  out = { Values, { "omega2tau", TAU_O2T_VR} };
  outH = { Values, { "omega2tau", TAU_O2T_H_VR} };
  outDefect = { Defect, { "omega2tau", TAU_O2T_VR} };

  for (size_t i = 0; i < results[out].size(); i++) {
    Sample outSample = results[out][i];
    Sample outHSample = results[outH][i];
    Sample outDefectSample = results[outDefect][i];
    assert(are_approx_equal(outSample.derivatives.values[0], outHSample.derivatives.values[0], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[1], outHSample.derivatives.values[1], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[2], outHSample.derivatives.values[2], 1e-6));
  }
}

TEST(Tests, test_cvode_defects) {
  CosimulationNetworkPtr network = examples::instantiate_oscillator_network(2, 2, examples::SimpleSolver::Cvode);

  Results results = jacobi_co_simulation(ConstantStepSize(0.25, 10.), Both, network);

  Signal out{ Values, { "tau2omega", OMEGA_T2O_VR} };
  Signal outH{ Values, { "tau2omega", OMEGA_T2O_H_VR} };
  Signal outDefect{ Defect, { "tau2omega", OMEGA_T2O_VR} };

  for (size_t i = 0; i < results[out].size(); i++) {
    Sample outSample = results[out][i];
    Sample outHSample = results[outH][i];
    Sample outDefectSample = results[outDefect][i];
    assert(are_approx_equal(outSample.derivatives.values[0], outHSample.derivatives.values[0], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[1], outHSample.derivatives.values[1], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[2], outHSample.derivatives.values[2], 1e-6));
  }

  out = { Values, { "omega2tau", TAU_O2T_VR} };
  outH = { Values, { "omega2tau", TAU_O2T_H_VR} };
  outDefect = { Defect, { "omega2tau", TAU_O2T_VR} };

  for (size_t i = 0; i < results[out].size(); i++) {
    Sample outSample = results[out][i];
    Sample outHSample = results[outH][i];
    Sample outDefectSample = results[outDefect][i];
    assert(are_approx_equal(outSample.derivatives.values[0], outHSample.derivatives.values[0], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[1], outHSample.derivatives.values[1], 1e-6));
    assert(are_approx_equal(outSample.derivatives.values[2], outHSample.derivatives.values[2], 1e-6));
  }
}

TEST(Tests, test_variable_step) {
  double h0 = 1e-3;
  double tol = 0.1;
  double tEnd = 3;
  VariableStepSize stepSizeGenerator(h0, tEnd, tol);
  double h = stepSizeGenerator.initial();
  while (!stepSizeGenerator.empty())
  {
    double eps = tol;
    h = stepSizeGenerator.next(eps);
    assert(are_approx_equal(h, h0, 1e-6));
  }
}

TEST(Tests, test_monolithic) {
  CosimulationNetworkPtr monolithic = fmi::examples::instantiate_twomass(2, 2);
  std::shared_ptr<FMU> tm = monolithic->instances["twomass"];
  tm->step(100);
  OutputDerivatives omega2M = tm->get_output(OMEGA_2M_VR);
  OutputDerivatives tau2M = tm->get_output(TAU_2M_VR);

  assert(are_approx_equal(omega2M.values[0], 0., 1e-3));
  assert(are_approx_equal(tau2M.values[0], 0., 1e-3));

}

TEST(Tests, test_values) {
  CosimulationNetworkPtr monolithic = fmi::examples::instantiate_twomass(2, 2);
  CosimulationNetworkPtr partitioned = examples::instantiate_oscillator_network(2);
  std::shared_ptr<FMU> tm = monolithic->instances["twomass"];
  std::shared_ptr<FMU> t2o = partitioned->instances["tau2omega"];
  std::shared_ptr<FMU> o2t = partitioned->instances["omega2tau"];

  double h = 1e-3;

  tm->step(h);
  t2o->step(h);
  o2t->step(h);

  OutputDerivatives omega2M = tm->get_output(OMEGA_2M_VR);
  OutputDerivatives tau2M = tm->get_output(TAU_2M_VR);

  OutputDerivatives omegaT2O = t2o->get_output(OMEGA_T2O_VR);
  OutputDerivatives tauO2T = o2t->get_output(TAU_O2T_VR);

  
  for (size_t i = 0; i < 3; i++) {
    assert(are_approx_equal(omega2M.values[i], omegaT2O.values[i], 1e-3));
    assert(are_approx_equal(tau2M.values[i], tauO2T.values[i], 1e-3));
  }
}

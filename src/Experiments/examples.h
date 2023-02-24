#pragma once

#include "fmi.h"

#define OMEGA_2M_VR 3
#define TAU_2M_VR 16

#define OMEGA_O2T_VR 0
#define TAU_O2T_VR 1
#define TAU_O2T_H_VR 11

#define TAU_T2O_VR 0
#define OMEGA_T2O_VR 1
#define OMEGA_T2O_H_VR 8


namespace fmi {
  namespace examples {

    CosimulationNetworkPtr instantiate_twomass(size_t extrapolationOrder);

    CosimulationNetworkPtr instantiate_twomass(size_t inputExtrapolationOrder, size_t outputExtrapolationOrder);

    enum SimpleSolver { Cvode = 0, Euler = 1 };

    CosimulationNetworkPtr instantiate_oscillator_network(size_t extrapolationOrder);

    CosimulationNetworkPtr instantiate_oscillator_network(size_t inputExtrapolationOrder, size_t outputExtrapolationOrder);

    CosimulationNetworkPtr instantiate_oscillator_network(size_t inputExtrapolationOrder, size_t outputExtrapolationOrder, SimpleSolver solver);

    template <typename StepGenerator>
    Results analytic_experiment(StepGenerator stepGenerator, DefectType defectType, size_t inputExtrapolationOrder, size_t outputExtrapolationOrder) {
      CosimulationNetworkPtr network = instantiate_twomass(inputExtrapolationOrder, outputExtrapolationOrder);
      return fmi::jacobi_co_simulation(stepGenerator, defectType, network);
    }

    template <typename StepGenerator>
    Results analytic_experiment(StepGenerator stepGenerator, DefectType defectType, size_t extrapolationOrder) {
      return analytic_experiment(stepGenerator, defectType, extrapolationOrder, extrapolationOrder);
    }

    template <typename StepGenerator>
    Results simple_experiment(StepGenerator stepGenerator, double h0, double tEnd, DefectType defectType, size_t extrapolationOrder) {
      CosimulationNetworkPtr network = instantiate_oscillator_network(extrapolationOrder);
      return fmi::jacobi_co_simulation(stepGenerator, h0, tEnd, defectType, network);
    }

    template <typename StepGenerator>
    void simple_example_run(StepGenerator stepGenerator, double h, double tEnd, size_t extrapolationOrder, const char* csvFile) {
      Results resultsData = simple_experiment(stepGenerator, h, tEnd, Both, extrapolationOrder);
      fmi::results::store_results(resultsData, csvFile);
    }

    template <typename StepGenerator>
    void analytic_example_run(StepGenerator stepGenerator, double h, double tEnd, size_t extrapolationOrder, const char* csvFile) {
      Results resultsData = examples::analytic_experiment(stepGenerator, h, tEnd, Both, extrapolationOrder);
      fmi::results::store_results(resultsData, csvFile);
    }
  }
}

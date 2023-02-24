#include "article.h"
#include "examples.h"
#include "results.h"
#include "fmi.h"
#include "signal.h"

#include <map>
#include <string>
#include <fstream>
#include <sstream>

#include <iostream>
#include <algorithm>
#include <cmath>

using namespace fmi;
using namespace fmi::examples;

template <typename StepGenerator>
static void produce_and_store_results(StepGenerator stepGenerator, CosimulationNetworkPtr network, const char *networkAnnotation, size_t extrapolationOrder)
{

  std::string csvFile = std::string(networkAnnotation) + "_" + std::to_string(extrapolationOrder) + ".csv";

  Results results = jacobi_co_simulation(stepGenerator, DefectType::Both, network);
  results::store_results(results, csvFile.c_str());
}

static std::map<std::string, CosimulationNetworkPtr> generate_networks(size_t inputExtrapolationOrder, size_t outputExtrapolationOrder)
{
  return {
      {"euler_subsystem", instantiate_oscillator_network(inputExtrapolationOrder, outputExtrapolationOrder, SimpleSolver::Euler)},
      {"cvode_subsystem", instantiate_oscillator_network(inputExtrapolationOrder, outputExtrapolationOrder, SimpleSolver::Cvode)},
  };
}

static std::map<std::string, CosimulationNetworkPtr> generate_networks(size_t extrapolationOrder)
{
  return generate_networks(extrapolationOrder, extrapolationOrder);
}

template <typename StepGenerator, typename ReferenceStepGenerator>
static void simulate_three_extrapolation_orders(StepGenerator stepGenerator, ReferenceStepGenerator referenceStepGenerator)
{
  for (size_t extrapolationOrder = 0; extrapolationOrder <= 2; extrapolationOrder++)
  {
    std::map<std::string, CosimulationNetworkPtr> networks = generate_networks(extrapolationOrder);

    for (auto it = networks.begin(); it != networks.end(); it++)
    {
      const char *networkAnnotation = it->first.c_str();
      CosimulationNetworkPtr network = it->second;
      produce_and_store_results(stepGenerator, network, networkAnnotation, extrapolationOrder);
    }
  }
  produce_and_store_results(referenceStepGenerator, instantiate_twomass(2), "cvode_monolithic", 2);
}

typedef void MyDirectoryAction();

static void execute_in_new_directory(const char *directory, MyDirectoryAction action)
{
#ifdef __linux__
#elif _WIN32
  char currentDirectory[MAX_PATH];
  GetCurrentDirectoryA(MAX_PATH, currentDirectory);
  CreateDirectoryA(directory, NULL);
  SetCurrentDirectoryA(directory);
  action();
  SetCurrentDirectoryA(currentDirectory);
#endif
}

static void experiment_1_simulation()
{
  simulate_three_extrapolation_orders(ConstantStepSize(1., 50.), ConstantStepSize(0.5, 50.));
}

void fmi::examples::article::experiment_1()
{
  execute_in_new_directory("experiment_1", experiment_1_simulation);
}

static void experiment_2_simulation()
{
  simulate_three_extrapolation_orders(VariableStepSize(1., 50., 0.005), ConstantStepSize(0.5, 50.));
}

void fmi::examples::article::experiment_2()
{
  execute_in_new_directory("experiment_2", experiment_2_simulation);
}

struct PortDescription
{
  std::string name;
  size_t inputExtrapolationOrder;
  size_t outputExtrapolationOrder;
  Port port;
  bool operator<(const PortDescription& rhs) const {
    if (name < rhs.name)
    {
      return true;
    }
    else if (name > rhs.name)
    {
      return false;
    }
    else if (inputExtrapolationOrder < rhs.inputExtrapolationOrder)
    {
      return true;
    }
    else if (inputExtrapolationOrder > rhs.inputExtrapolationOrder)
    {
      return false;
    }
    else if (outputExtrapolationOrder < rhs.outputExtrapolationOrder)
    {
      return true;
    }
    else if (outputExtrapolationOrder > rhs.outputExtrapolationOrder)
    {
      return false;
    }
    else if (port < rhs.port)
    {
      return true;
    }
    else
    {
      return false;
    }
  }
};

struct StepSizeStats {
  double min, max, avg;
};

struct StepSizeInfo
{
  std::vector<double> xAxis;
  std::vector<StepSizeStats> values;
};

struct ErrorInfo
{
  std::vector<double> xAxis;
  std::vector<double> values;
};

static void find_errors_and_defects(
  std::map<std::string, CosimulationNetworkPtr> networks,
  double stepSize, double tEnd,
  size_t inputExtrapolationOrder, size_t outputExtrapolationOrder,
  std::map<PortDescription, ErrorInfo>& errors,
  std::map<PortDescription, ErrorInfo>& defects
)
{
  for (auto networkIt = networks.begin(); networkIt != networks.end(); networkIt++)
  {

    CosimulationNetworkPtr network = networkIt->second;

    CosimulationNetworkPtr monolithic = instantiate_twomass(2);
    Results networkResults = jacobi_co_simulation(ConstantStepSize(stepSize, tEnd), DefectType::Both, network);
    Results monolithicResults = jacobi_co_simulation(ConstantStepSize(stepSize, tEnd), DefectType::Both, monolithic);

    Port t2oOmega = {"tau2omega", OMEGA_T2O_VR};
    Port twomassOmega = {"twomass", OMEGA_2M_VR};
    PortDescription t2oOmegaErrId = {networkIt->first, inputExtrapolationOrder, outputExtrapolationOrder, t2oOmega };
    errors[t2oOmegaErrId].xAxis.push_back(stepSize);
    errors[t2oOmegaErrId].values.push_back(
      fmi::results::calculate_error(
        networkResults[{fmi::SignalType::Values, t2oOmega}],
        monolithicResults[{fmi::SignalType::Values, twomassOmega}])
    );

    defects[t2oOmegaErrId].xAxis.push_back(stepSize);
    defects[t2oOmegaErrId].values.push_back(fmi::signal::rms(networkResults[{fmi::SignalType::Defect, t2oOmega}]));

    Port t2oTau = { "tau2omega", TAU_T2O_VR };
    PortDescription t2oTauErrId = { networkIt->first, inputExtrapolationOrder, outputExtrapolationOrder, t2oTau };
    defects[t2oTauErrId].xAxis.push_back(stepSize);
    defects[t2oTauErrId].values.push_back(fmi::signal::rms(networkResults[{fmi::SignalType::Defect, t2oTau}]));

    Port o2tTau = {"omega2tau", TAU_O2T_VR};
    Port twomassTau = {"twomass", TAU_2M_VR};
    PortDescription o2tTauErrId = {networkIt->first, inputExtrapolationOrder, outputExtrapolationOrder, o2tTau };
    errors[o2tTauErrId].xAxis.push_back(stepSize);
    errors[o2tTauErrId].values.push_back(
      fmi::results::calculate_error(
        networkResults[{fmi::SignalType::Values, o2tTau}],
        monolithicResults[{fmi::SignalType::Values, twomassTau}])
      );

    defects[o2tTauErrId].xAxis.push_back(stepSize);
    defects[o2tTauErrId].values.push_back(fmi::signal::rms(networkResults[{fmi::SignalType::Defect, o2tTau}]));

    Port o2tOmega = { "omega2tau", OMEGA_O2T_VR };
    PortDescription o2tOmegaErrId = { networkIt->first, inputExtrapolationOrder, outputExtrapolationOrder, o2tOmega };
    defects[o2tOmegaErrId].xAxis.push_back(stepSize);
    defects[o2tOmegaErrId].values.push_back(fmi::signal::rms(networkResults[{fmi::SignalType::Defect, o2tOmega}]));
  }
}

static StepSizeStats analyze_step_size(StepSequence& stepSequence)
{
  double h = stepSequence.initial();
  StepSizeStats stats = { h, h, h };
  int n = 1;

  while (!stepSequence.empty())
  {
    h = stepSequence.next(1.);
    if (h < stats.min)
    {
      stats.min = h;
    }
    if (h > stats.max)
    {
      stats.max = h;
    }
    stats.avg += h;
    n ++;
  }
  stats.avg /= n;
  return stats;
}

static void analyze_variable_step(
  std::map<std::string, CosimulationNetworkPtr> networks,
  double h0, double tol, double tEnd, size_t extrapolationOrder,
  std::map<PortDescription, StepSizeInfo>& stepSizeStats,
  std::map<PortDescription, ErrorInfo>& errors,
  std::map<PortDescription, ErrorInfo>& defects
)
{
  for (auto networkIt = networks.begin(); networkIt != networks.end(); networkIt++)
  {

    CosimulationNetworkPtr network = networkIt->second;

    CosimulationNetworkPtr analytic = instantiate_twomass(2);
    Results results = jacobi_co_simulation(VariableStepSize(h0, tEnd, tol), DefectType::Both, network);

    Port t2oOmega = {"tau2omega", OMEGA_T2O_VR};
    StepSequence stepSequence(results[{fmi::SignalType::Values, t2oOmega}]);
    StepSequence stepSequenceForAnalysis(results[{fmi::SignalType::Values, t2oOmega}]);

    Results analyticResults = jacobi_co_simulation(stepSequence, DefectType::Both, analytic);


    PortDescription t2oOmegaErrId = {networkIt->first, extrapolationOrder, extrapolationOrder, t2oOmega };
    errors[t2oOmegaErrId].xAxis.push_back(tol);
    errors[t2oOmegaErrId].values.push_back(
      fmi::results::calculate_error(
        results[{fmi::SignalType::Values, t2oOmega}],
        analyticResults[{fmi::SignalType::Values, {"twomass", OMEGA_2M_VR}}])
    );

    stepSizeStats[t2oOmegaErrId].xAxis.push_back(tol);
    stepSizeStats[t2oOmegaErrId].values.push_back(analyze_step_size(stepSequenceForAnalysis));

    defects[t2oOmegaErrId].xAxis.push_back(tol);
    defects[t2oOmegaErrId].values.push_back(fmi::signal::rms(results[{fmi::SignalType::Defect, t2oOmega}]));

    Port t2oTau = { "tau2omega", 0 };
    PortDescription t2oTauErrId = { networkIt->first, extrapolationOrder, extrapolationOrder, t2oTau };
    defects[t2oTauErrId].xAxis.push_back(tol);
    defects[t2oTauErrId].values.push_back(fmi::signal::rms(results[{fmi::SignalType::Defect, t2oTau}]));

    Port o2tTau = {"omega2tau", TAU_O2T_VR};
    PortDescription o2tTauErrId = {networkIt->first, extrapolationOrder, extrapolationOrder, o2tTau };
    errors[o2tTauErrId].xAxis.push_back(tol);
    errors[o2tTauErrId].values.push_back(
      fmi::results::calculate_error(
        results[{fmi::SignalType::Values, o2tTau}],
        analyticResults[{fmi::SignalType::Values, {"twomass", TAU_2M_VR}}])
      );

    defects[o2tTauErrId].xAxis.push_back(tol);
    defects[o2tTauErrId].values.push_back(fmi::signal::rms(results[{fmi::SignalType::Defect, o2tTau}]));

    Port o2tOmega = { "omega2tau", OMEGA_O2T_VR };
    PortDescription o2tOmegaErrId = { networkIt->first, extrapolationOrder, extrapolationOrder, o2tOmega };
    defects[o2tOmegaErrId].xAxis.push_back(tol);
    defects[o2tOmegaErrId].values.push_back(fmi::signal::rms(results[{fmi::SignalType::Defect, o2tOmega}]));
  }
}

static void store_errors_and_defects(
  const std::map<PortDescription, ErrorInfo>& errors,
  const std::map<PortDescription, ErrorInfo>& defects
)
{
  const char* delimiter = ",";
  for (auto errs: errors)
  {
    std::string csvFileName = "errors_" + errs.first.name + "_" + errs.first.port.instanceName + "_" + std::to_string(errs.first.port.portReference)
      + "_" + std::to_string(errs.first.inputExtrapolationOrder) + "_" + std::to_string(errs.first.outputExtrapolationOrder) + ".csv";
    std::ofstream errFile(csvFileName);
    std::stringstream ssStepSizes, ssValues;
    std::copy(errs.second.xAxis.begin(), errs.second.xAxis.end(), std::ostream_iterator<double>(ssStepSizes, delimiter));
    std::copy(errs.second.values.begin(), errs.second.values.end(), std::ostream_iterator<double>(ssValues, delimiter));
    errFile << ssStepSizes.str() << std::endl << ssValues.str();
  }

  for (auto errs : defects)
  {
      std::string csvFileName = "defects_" + errs.first.name + "_" + errs.first.port.instanceName + "_" + std::to_string(errs.first.port.portReference)
          + "_" + std::to_string(errs.first.inputExtrapolationOrder) + "_" + std::to_string(errs.first.outputExtrapolationOrder) + ".csv";
      std::ofstream errFile(csvFileName);
      std::stringstream ssStepSizes, ssValues;
      std::copy(errs.second.xAxis.begin(), errs.second.xAxis.end(), std::ostream_iterator<double>(ssStepSizes, delimiter));
      std::copy(errs.second.values.begin(), errs.second.values.end(), std::ostream_iterator<double>(ssValues, delimiter));
      errFile << ssStepSizes.str() << std::endl << ssValues.str();
  }

}

static void experiment_3_simulation()
{
  std::vector<double> stepSizes;
  for (double expStepSize = -3.; expStepSize < 0.1; expStepSize += 0.2)
  {
    stepSizes.push_back(pow(10, expStepSize));
  }

  std::map<PortDescription, ErrorInfo> errors;
  std::map<PortDescription, ErrorInfo> defects;
  double tEnd = 1.;

  for (size_t extrapolationOrder = 0; extrapolationOrder <= 2; extrapolationOrder++)
  {
    size_t inputExtrapolationOrder = extrapolationOrder;
    size_t outputExtrapolationOrder = extrapolationOrder;
    for (auto stepSize = stepSizes.begin(); stepSize != stepSizes.end(); stepSize++)
    {
      std::map<std::string, CosimulationNetworkPtr> networks = generate_networks(extrapolationOrder);
      find_errors_and_defects(
        networks, *stepSize, tEnd,
        inputExtrapolationOrder, outputExtrapolationOrder,
        errors, defects
      );
    }
  }
  store_errors_and_defects(errors, defects);
}

static void experiment_4_simulation()
{
  std::vector<double> stepSizes;
  for (double expStepSize = -3.; expStepSize < 0.1; expStepSize += 0.2)
  {
    stepSizes.push_back(pow(10, expStepSize));
  }

  std::map<PortDescription, ErrorInfo> errors;
  std::map<PortDescription, ErrorInfo> defects;
  size_t inputExtrapolationOrder = 0;
  double tEnd = 1.;

  for (size_t outputExtrapolationOrder = 0; outputExtrapolationOrder <= 2; outputExtrapolationOrder++)
  {
    for (auto stepSize = stepSizes.begin(); stepSize != stepSizes.end(); stepSize++)
    {
      std::map<std::string, CosimulationNetworkPtr> networks = generate_networks(inputExtrapolationOrder, outputExtrapolationOrder);
      find_errors_and_defects(
        networks, *stepSize, tEnd,
        inputExtrapolationOrder, outputExtrapolationOrder,
        errors, defects
      );
    }
  }
  store_errors_and_defects(errors, defects);
}

static void experiment_5_simulation()
{
  std::vector<double> stepSizes;
  for (double expStepSize = -3.; expStepSize < 0.1; expStepSize += 0.2)
  {
    stepSizes.push_back(pow(10, expStepSize));
  }

  std::map<PortDescription, ErrorInfo> errors;
  std::map<PortDescription, ErrorInfo> defects;
  size_t outputExtrapolationOrder = 0;
  double tEnd = 1.;

  for (size_t inputExtrapolationOrder = 0; inputExtrapolationOrder <= 2; inputExtrapolationOrder++)
  {
    for (auto stepSize = stepSizes.begin(); stepSize != stepSizes.end(); stepSize++)
    {
      std::map<std::string, CosimulationNetworkPtr> networks = generate_networks(inputExtrapolationOrder, outputExtrapolationOrder);
      find_errors_and_defects(
        networks, *stepSize, tEnd,
        inputExtrapolationOrder, outputExtrapolationOrder,
        errors, defects
      );
    }
  }
  store_errors_and_defects(errors, defects);
}


static void store_step_size_stats(std::map<PortDescription, StepSizeInfo> stepSizeStats)
{
  const char* delimiter = ",";
  for (auto stats: stepSizeStats)
  {
    std::string csvFileName = "step_size_" + stats.first.name + "_" + std::to_string(stats.first.inputExtrapolationOrder) + "_" + std::to_string(stats.first.outputExtrapolationOrder) + ".csv";
    std::ofstream errFile(csvFileName);
    std::stringstream ssStepSizes, ssMins, ssMaxs, ssAvgs;
    std::copy(stats.second.xAxis.begin(), stats.second.xAxis.end(), std::ostream_iterator<double>(ssStepSizes, delimiter));
    std::vector<double> havgs;
    std::transform(stats.second.values.begin(), stats.second.values.end(), std::back_inserter(havgs), [](StepSizeStats& ssStats) { return ssStats.avg; });
    std::copy(havgs.begin(), havgs.end(), std::ostream_iterator<double>(ssAvgs, delimiter));
    std::vector<double> hmaxs;
    std::transform(stats.second.values.begin(), stats.second.values.end(), std::back_inserter(hmaxs), [](StepSizeStats& ssStats) { return ssStats.max; });
    std::copy(hmaxs.begin(), hmaxs.end(), std::ostream_iterator<double>(ssMaxs, delimiter));
    std::vector<double> hmins;
    std::transform(stats.second.values.begin(), stats.second.values.end(), std::back_inserter(hmins), [](StepSizeStats& ssStats) { return ssStats.min; });
    std::copy(hmins.begin(), hmins.end(), std::ostream_iterator<double>(ssMins, delimiter));
    errFile << ssStepSizes.str() << std::endl << ssMins.str() << std::endl << ssAvgs.str() << std::endl << ssMaxs.str();
  }
}


static void experiment_6_simulation()
{
  std::vector<double> tolerances;
  for (double expTol = -3.; expTol < 0.1; expTol += 0.2)
  {
    tolerances.push_back(pow(10, expTol));
  }

  std::map<PortDescription, ErrorInfo> errors;
  std::map<PortDescription, ErrorInfo> defects;
  std::map<PortDescription, StepSizeInfo> stepSizeStats;
  double tEnd = 1.;

  for (size_t extrapolationOrder = 0; extrapolationOrder <= 2; extrapolationOrder++)
  {
    for (auto tol = tolerances.begin(); tol != tolerances.end(); tol++)
    {
      std::map<std::string, CosimulationNetworkPtr> networks = generate_networks(extrapolationOrder);
      analyze_variable_step(
        networks, 1e-4, *tol, tEnd,
        extrapolationOrder,
        stepSizeStats, errors, defects
      );
    }
  }
  store_errors_and_defects(errors, defects);
  store_step_size_stats(stepSizeStats);
}

void fmi::examples::article::experiment_3()
{
  execute_in_new_directory("experiment_3", experiment_3_simulation);
}

void fmi::examples::article::experiment_4()
{
  execute_in_new_directory("experiment_4", experiment_4_simulation);
}

void fmi::examples::article::experiment_5()
{
  execute_in_new_directory("experiment_5", experiment_5_simulation);
}

void fmi::examples::article::experiment_6()
{
  execute_in_new_directory("experiment_6", experiment_6_simulation);
}

void fmi::examples::article::results_for_all_figures()
{
  fmi::examples::article::experiment_1();
  fmi::examples::article::experiment_2();
  fmi::examples::article::experiment_3();
  fmi::examples::article::experiment_4();
  fmi::examples::article::experiment_5();
  fmi::examples::article::experiment_6();
}

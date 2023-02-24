#include "examples.h"
#include "results.h"

#include <iostream>

#ifdef __linux__
#include <unistd.h>
#include <limits.h>
#include <libgen.h>

static std::string exe_directory() {
  char buffer[PATH_MAX];
  ssize_t len = readlink("/proc/self/exe", buffer, sizeof(buffer)-1);
  buffer[len] = '\0';
  return std::string(dirname(buffer));
}

static std::string omega2tau_library() {
  return exe_directory() + "/libOscillatorOmega2Tau.so";
}

static std::string tau2omega_library() {
  return exe_directory() + "/libOscillatorTau2Omega.so";
}

static std::string two_mass_rotational_oscillator_library() {
  return exe_directory() + "/libTwoMassRotationalOscillator.so";
}

#elif _WIN32

static std::string omega2tau_library() {
  return "OscillatorOmega2Tau.dll";
}

static std::string tau2omega_library() {
  return "OscillatorTau2Omega.dll";
}

static std::string two_mass_rotational_oscillator_library() {
  return "TwoMassRotationalOscillator.dll";
}

#endif

using namespace fmi;

CosimulationNetworkPtr fmi::examples::instantiate_oscillator_network(size_t extrapolationOrder) {
  return instantiate_oscillator_network(extrapolationOrder, extrapolationOrder, Cvode);
}

CosimulationNetworkPtr fmi::examples::instantiate_oscillator_network(size_t inputExtrapolationOrder, size_t outputExtrapolationOrder)
{
  return instantiate_oscillator_network(inputExtrapolationOrder, outputExtrapolationOrder, Cvode);
}

static void initialize_simple_network(CosimulationNetworkPtr simpleNetwork, size_t inputExtrapolationOrder, size_t outputExtrapolationOrder)
{
  CosimulationNetworkPtr analytic = fmi::examples::instantiate_twomass(inputExtrapolationOrder, outputExtrapolationOrder);
  simpleNetwork->instances["tau2omega"]->set_input(TAU_T2O_VR,
    analytic->instances["twomass"]->get_output(TAU_2M_VR)
  );
  simpleNetwork->instances["omega2tau"]->set_input(OMEGA_O2T_VR,
    analytic->instances["twomass"]->get_output(OMEGA_2M_VR)
  );
}

CosimulationNetworkPtr fmi::examples::instantiate_oscillator_network(size_t inputExtrapolationOrder, size_t outputExtrapolationOrder, SimpleSolver solver)
{
  CosimulationNetworkPtr network = std::make_shared<CosimulationNetwork>();

  InstanceInfo omega2tau = {
    "omega2tau", omega2tau_library(), "resource",
    "{3063ef2f-10cf-48fd-898e-833e2d41174b}", inputExtrapolationOrder, outputExtrapolationOrder,
    { 
      {2, 10}, 
      {3, 1},
      {4, 1},
      {5, 1},
      {6, 2},
      {7, 0.1},
      {8, 0.1},
      {9, 0.2},
      {10, (fmi2Integer)solver}
    }
  };
  network->instances[omega2tau.instanceName] = std::make_shared<FMU>(omega2tau);

  InstanceInfo tau2omega = {
    "tau2omega", tau2omega_library(), "resource",
    "{3063ef2f-10cf-48fd-898e-833e2d41174a}", inputExtrapolationOrder, outputExtrapolationOrder,
    {
      {2, 10},
      {3, 1},
      {4, 2},
      {5, 0.2},
      {6, 0.1},
      {7, (fmi2Integer)solver}
    }
  };
  network->instances[tau2omega.instanceName] = std::make_shared<FMU>(tau2omega);


  Port input1{ "tau2omega", TAU_T2O_VR}, input2{ "omega2tau", OMEGA_O2T_VR };
  Port output1{ "tau2omega", OMEGA_T2O_VR }, output2{ "omega2tau", TAU_O2T_VR };
  OptionalOutput output1H{ { "tau2omega", OMEGA_T2O_H_VR}, true }, output2H{ { "omega2tau", TAU_O2T_H_VR }, true };

  network->connections.emplace(input1, output2);
  network->monitoredOutputs.emplace(output2, output2H);

  network->connections.emplace(input2, output1);
  network->monitoredOutputs.emplace(output1, output1H);

  initialize_simple_network(network, inputExtrapolationOrder, outputExtrapolationOrder);

  return network;
}

CosimulationNetworkPtr fmi::examples::instantiate_twomass(size_t extrapolationOrder) {
  return instantiate_twomass(extrapolationOrder, extrapolationOrder);
}

CosimulationNetworkPtr fmi::examples::instantiate_twomass(size_t inputExtrapolationOrder, size_t outputExtrapolationOrder) {
  CosimulationNetworkPtr network = std::make_shared<CosimulationNetwork>();
  RealParameters parameters = {
    {4, 10},
    {5, 1},
    {6, 1},
    {7, 0.1},
    {8, 0.1},
    {9, 10},
    {10, 1},
    {11, 2},
    {12, 0.2},
    {13, 0.1},
    {14, 1},
    {15, 2},
  };
  std::map<fmi2ValueReference, size_t> inputExtrapolationOrders;
  std::map<fmi2ValueReference, size_t> outputExtrapolationOrders = {
    {3, outputExtrapolationOrder}, {16, outputExtrapolationOrder},
    {17, outputExtrapolationOrder + 1}, {18, outputExtrapolationOrder + 1}
  };
  IntegerParameters intParameters;
  InstanceInfo info = { "twomass", two_mass_rotational_oscillator_library().c_str(), "resource", "{799f52bd-d05b-4e92-8325-81d2bf6f7b09}", inputExtrapolationOrder, outputExtrapolationOrder, parameters, intParameters };
  FMU simple(info);
  network->instances[info.instanceName] = std::make_shared<FMU>(info);
  Port output1{ "twomass", 3 }, output2{ "twomass", 16 };
  OptionalOutput output1H{ { "twomass", 17 }, true }, output2H{ {"twomass", 18 }, true };
  network->monitoredOutputs.emplace(output1, output1H);
  network->monitoredOutputs.emplace(output2, output2H);
  return network;
}

#ifndef FMI_H
#define FMI_H

#include <fmi2TypesPlatform.h>
#include <fmi2FunctionTypes.h>

#ifdef __linux__
#include <dlfcn.h>
#elif _WIN32
#define NOMINMAX
#include <Windows.h>
#endif


#include <map>
#include <vector>
#include <queue>
#include <memory>
#include <numeric>
#include <algorithm>
#include <string>
#include <iterator>

namespace fmi {
  typedef std::map<fmi2ValueReference, fmi2Real> RealParameters;
  typedef std::map<fmi2ValueReference, fmi2Integer> IntegerParameters;
  struct Derivatives {
    std::vector<fmi2Real> values;
    double t;
  };
  typedef Derivatives InputDerivatives;
  typedef Derivatives OutputDerivatives;

  struct InstanceInfo {
    std::string instanceName;
    std::string dllPath;
    std::string resourceLocation;
    std::string guid;
    size_t inputExtrapolationOrder;
    size_t outputExtrapolationOrder;
    RealParameters realParameters;
    IntegerParameters intParameters;
  };

  struct Port {
    std::string instanceName;
    fmi2ValueReference portReference;
    bool operator<(const Port &other)  const;
  };

  typedef Port Output;
  typedef Port Input;

  typedef double Timestamp;
  struct Interval {
    Timestamp begin, end;
  };
  struct Sample {
    Interval interval;
    Derivatives derivatives;
  };
  enum SignalType { Values = 0, Defect = 1 };
  struct Signal {
    SignalType type;
    Port port;
    bool operator<(const Signal& rhs) const;
  };
  typedef std::vector<Sample> Samples;
  typedef std::map<Signal, Samples> Results;

  class FMU {
#ifdef __linux__
    void* _handle;
#elif _WIN32
    HINSTANCE _handle;
#endif
    fmi2Component _component;
    fmi2ComponentEnvironment _componentEnvironment;
    fmi2ExitInitializationModeTYPE* _exitInitializationMode;
    fmi2DoStepTYPE* _doStep;
    fmi2SetRealTYPE* _setReal;
    fmi2SetRealInputDerivativesTYPE* _setRealInputDerivatives;
    fmi2SetIntegerTYPE* _setInteger;
    fmi2GetRealTYPE* _getReal;
    fmi2GetIntegerTYPE* _getInteger;
    fmi2GetRealOutputDerivativesTYPE* _getRealOutputDerivatives;
    fmi2TerminateTYPE* _terminate;
    fmi2FreeInstanceTYPE* _freeInstance;
    double _time;
    size_t _inputExtrapolationOrder;
    size_t _outputExtrapolationOrder;
    bool _initialized;
    FMU(const FMU& other) = delete; // diabling copies
    FMU& operator=(const FMU& other) = delete;
  public:
    FMU(const InstanceInfo& info);
    InputDerivatives set_input(fmi2ValueReference reference, const InputDerivatives& values);
    OutputDerivatives get_output(fmi2ValueReference reference);
    OutputDerivatives get_hermite(fmi2ValueReference reference, fmi2ValueReference midReference, double h);
    void step(double h);
    double get_time() const;
    ~FMU();
  };

  struct OptionalOutput
  {
    Output output;
    bool exists;
  };

  typedef std::map<std::string, std::shared_ptr<FMU>> Instances;
  typedef std::map<Input, Output> Connections;
  typedef std::map<Output, OptionalOutput> MonitoredOutputs;

  struct CosimulationNetwork {
    Instances instances;
    Connections connections;
    MonitoredOutputs monitoredOutputs;
    double t;
  };

  typedef std::shared_ptr<CosimulationNetwork> CosimulationNetworkPtr;

  class ConstantStepSize {
    double _tEnd, _h0;
    double _h, _t;
  public:
    ConstantStepSize(double h0, double tEnd);
    bool empty();
    double initial();
    double next(double eps);
  };

  class VariableStepSize {
    double _KP, _KI, _tol, _phimax, _I;
    double _h, _h0, _tEnd, _t;
  public:
    VariableStepSize(double h0, double tEnd, double tol);
    bool empty();
    double initial();
    double next(double eps);
  };

  class StepSequence {
    std::queue<double> _stepSequence;
  public:
    StepSequence(const Samples& samples);
    bool empty();
    double initial();
    double next(double eps);
  };

  void exchange_values(CosimulationNetworkPtr network, double h, bool recordInputs, Results& results);

  void record_outputs(CosimulationNetworkPtr network, double h, Results& results);

  void record_output_defects(Results& results, const MonitoredOutputs& monitoredOutputs);

  void record_connection_defects(Results& results, const Connections& connections);

  enum DefectType { OutputDefect, ConnectionDefect, Both };

  double calculate_output_defects(const Results& results, const MonitoredOutputs& defectOutputs);

  double calculate_connection_defects(const Results& results, const Connections& connections);

  template <typename StepGenerator>
  Results jacobi_co_simulation(StepGenerator stepGenerator, DefectType defectType, CosimulationNetworkPtr network) {
    Results results;
    std::vector<std::shared_ptr<FMU>> instances;
    std::transform(network->instances.begin(), network->instances.end(), std::back_inserter(instances),
      [](const std::pair<std::string, std::shared_ptr<FMU>>& kv) { return kv.second; }
    );

    double h = stepGenerator.initial();

    while (true)
    {
      exchange_values(network, h, true, results);
      for (auto it = instances.begin(); it != instances.end(); it++) {
        std::shared_ptr<FMU> instance = *it;
        instance->step(h);
      }
      record_outputs(network, h, results);
      network->t += h;

      record_connection_defects(results, network->connections);
      record_output_defects(results, network->monitoredOutputs);

      double eps = 0.;
      if (defectType == Both || defectType == OutputDefect) {
        eps = calculate_output_defects(results, network->monitoredOutputs);
      }
      if (defectType == Both || defectType == ConnectionDefect) {
      double connectionEps = calculate_connection_defects(results, network->connections);
        eps = std::max(eps, connectionEps);
      }
      if (stepGenerator.empty())
      {
          break;
      }
      h = stepGenerator.next(eps);
    }
    return results;
  }
}

#include <results.h>

#endif // FMI_H

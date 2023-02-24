#include "fmi.h"
#include "signal.h"

#include <cstring>
#include <set>
#include <cmath>
#include <cstdarg>

#ifdef __linux__
#include <dlfcn.h>
#endif

using namespace fmi;

#define MSG_BUFFER_MAX_SIZE 2048
#define STATUS_MAX_SIZE 30

void logger(fmi2ComponentEnvironment componentEnvironment, fmi2String instanceName, fmi2Status status, fmi2String category, fmi2String message, ...)
{
  char buffer[MSG_BUFFER_MAX_SIZE];
  char fmi2StatusString[STATUS_MAX_SIZE];

  switch (status)
  {
  case fmi2OK:
    strcpy(fmi2StatusString, "fmi2OK");
    break;
  case fmi2Warning:
    strcpy(fmi2StatusString, "fmi2Warning");
    break;
  case fmi2Discard:
    strcpy(fmi2StatusString, "fmi2Discard");
    break;
  case fmi2Error:
    strcpy(fmi2StatusString, "fmi2Error");
    break;
  case fmi2Fatal:
    strcpy(fmi2StatusString, "fmi2Fatal");
    break;
  case fmi2Pending:
    strcpy(fmi2StatusString, "fmi2Pending");
    break;
  default:
    break;
  }

  va_list args;
  va_start(args, message);

  vsnprintf(buffer, MSG_BUFFER_MAX_SIZE, message, args);
}

void* allocateMemory(size_t nobj, size_t size)
{
  return malloc(nobj * size);
}

void freeMemory(void* obj)
{
  free(obj);
}

fmi2CallbackFunctions callbackFunctions = { logger, allocateMemory, freeMemory, NULL, NULL };

FMU::FMU(const InstanceInfo& info)
  : _time(0), _initialized(false), _inputExtrapolationOrder(info.inputExtrapolationOrder), _outputExtrapolationOrder(info.outputExtrapolationOrder)
{
#ifdef __linux__
  _handle = dlopen(info.dllPath.c_str(), RTLD_NOW);

  fmi2InstantiateTYPE* instantiate = (fmi2InstantiateTYPE*)dlsym(_handle, "fmi2Instantiate");
  fmi2EnterInitializationModeTYPE* enterInitializationMode = (fmi2EnterInitializationModeTYPE*)dlsym(_handle, "fmi2EnterInitializationMode");
  _exitInitializationMode = (fmi2ExitInitializationModeTYPE*)dlsym(_handle, "fmi2ExitInitializationMode");
  _doStep = (fmi2DoStepTYPE*)dlsym(_handle, "fmi2DoStep");
  _terminate = (fmi2TerminateTYPE*)dlsym(_handle, "fmi2Terminate");
  _freeInstance = (fmi2FreeInstanceTYPE*)dlsym(_handle, "fmi2FreeInstance");

  _setReal = (fmi2SetRealTYPE*)dlsym(_handle, "fmi2SetReal");
  _setInteger = (fmi2SetIntegerTYPE*)dlsym(_handle, "fmi2SetInteger");
  _setRealInputDerivatives = (fmi2SetRealInputDerivativesTYPE*)dlsym(_handle, "fmi2SetRealInputDerivatives");
  _getReal = (fmi2GetRealTYPE*)dlsym(_handle, "fmi2GetReal");
  _getInteger = (fmi2GetIntegerTYPE*)dlsym(_handle, "fmi2GetInteger");
  _getRealOutputDerivatives = (fmi2GetRealOutputDerivativesTYPE*)dlsym(_handle, "fmi2GetRealOutputDerivatives");
#elif _WIN32
  _handle = LoadLibraryA(info.dllPath.c_str());

  fmi2InstantiateTYPE* instantiate = (fmi2InstantiateTYPE*)GetProcAddress(_handle, "fmi2Instantiate");
  fmi2EnterInitializationModeTYPE* enterInitializationMode = (fmi2EnterInitializationModeTYPE*)GetProcAddress(_handle, "fmi2EnterInitializationMode");
  _exitInitializationMode = (fmi2ExitInitializationModeTYPE*)GetProcAddress(_handle, "fmi2ExitInitializationMode");
  _doStep = (fmi2DoStepTYPE*)GetProcAddress(_handle, "fmi2DoStep");
  _terminate = (fmi2TerminateTYPE*)GetProcAddress(_handle, "fmi2Terminate");
  _freeInstance = (fmi2FreeInstanceTYPE*)GetProcAddress(_handle, "fmi2FreeInstance");

  _setReal = (fmi2SetRealTYPE*)GetProcAddress(_handle, "fmi2SetReal");
  _setInteger = (fmi2SetIntegerTYPE*)GetProcAddress(_handle, "fmi2SetInteger");
  _setRealInputDerivatives = (fmi2SetRealInputDerivativesTYPE*)GetProcAddress(_handle, "fmi2SetRealInputDerivatives");
  _getReal = (fmi2GetRealTYPE*)GetProcAddress(_handle, "fmi2GetReal");
  _getInteger = (fmi2GetIntegerTYPE*)GetProcAddress(_handle, "fmi2GetInteger");
  _getRealOutputDerivatives = (fmi2GetRealOutputDerivativesTYPE*)GetProcAddress(_handle, "fmi2GetRealOutputDerivatives");
#endif

  _component = instantiate(info.instanceName.c_str(), fmi2CoSimulation, info.guid.c_str(), info.resourceLocation.c_str(), &callbackFunctions, fmi2False, fmi2True);

  enterInitializationMode(_component);
  for (RealParameters::const_iterator it = info.realParameters.begin(); it != info.realParameters.end(); it++) {
    _setReal(_component, &it->first, 1, &it->second);
  }
  for (IntegerParameters::const_iterator it = info.intParameters.begin(); it != info.intParameters.end(); it++) {
    _setInteger(_component, &it->first, 1, &it->second);
  }
}

InputDerivatives FMU::set_input(fmi2ValueReference reference, const InputDerivatives& values) {
  InputDerivatives actuallySet = { std::vector<fmi2Real>(), values.t };
  size_t tmpOrder = _inputExtrapolationOrder <= values.values.size() - 1 ? _inputExtrapolationOrder : values.values.size() - 1;
  _setReal(_component, &reference, 1, &values.values[0]);
  actuallySet.values.push_back(values.values[0]);
  for (fmi2Integer i = 1; i <= tmpOrder; i++) {
    _setRealInputDerivatives(_component, &reference, 1, &i, &values.values[i]);
    actuallySet.values.push_back(values.values[i]);
  }
  fmi2Real zero = 0.;
  for (fmi2Integer i = tmpOrder + 1; i <= _inputExtrapolationOrder; i++) {
    _setRealInputDerivatives(_component, &reference, 1, &i, &zero);
  }
  return actuallySet;
}

OutputDerivatives FMU::get_output(fmi2ValueReference reference) {
  OutputDerivatives values = {
    std::vector<fmi2Real>(_outputExtrapolationOrder + 1),
    _time
  };
  _getReal(_component, &reference, 1, &values.values[0]);
  for (fmi2Integer i = 1; i <= _outputExtrapolationOrder; i++) {
    _getRealOutputDerivatives(_component, &reference, 1, &i, &values.values[i]);
  }
  return values;
}

OutputDerivatives FMU::get_hermite(fmi2ValueReference reference, fmi2ValueReference midReference, double h) {
  OutputDerivatives hermite = {
    std::vector<fmi2Real>(_outputExtrapolationOrder + 2),
    _time
  };
  fmi2Real yhMid;
  _getReal(_component, &reference, 1, &hermite.values[0]);
  _getReal(_component, &midReference, 1, &yhMid);
  for (fmi2Integer i = 1; i <= _outputExtrapolationOrder; i++) {
    _getRealOutputDerivatives(_component, &reference, 1, &i, &hermite.values[i]);
  }
  fmi2Real atMid = signal::eval(hermite, _time - h / 2);
  fmi2Real dyh = yhMid - atMid;
  for (size_t i = 1; i <= _outputExtrapolationOrder + 1; i++)
  {
    dyh *= i * 2 / h;
  }
  hermite.values[_outputExtrapolationOrder + 1] = -dyh;
  return hermite;
}

void FMU::step(double h) {
  if (!_initialized) {
    _initialized = true;
    _exitInitializationMode(_component);
  }
  _doStep(_component, _time, h, fmi2True);
  _time += h;
}

double FMU::get_time() const {
  return _time;
}

FMU::~FMU() {
#ifdef __linux__
  dlclose(_handle);
#elif _WIN32
  FreeLibrary(_handle);
#endif
}

bool Signal::operator< (const Signal& rhs) const {
  if (type < rhs.type) {
    return true;
  }
  else if (type > rhs.type) {
    return false;
  }
  else if (port.instanceName < rhs.port.instanceName) {
    return true;
  }
  else if (port.instanceName > rhs.port.instanceName) {
    return false;
  }
  else if (port.portReference < rhs.port.portReference) {
    return true;
  }
  else {
    return false;
  }
}


void fmi::exchange_values(CosimulationNetworkPtr network, double h, bool recordInputs, Results& results) {
  Interval interval = { network->t, network->t + h };

  for (auto kv = network->connections.begin(); kv != network->connections.end(); kv++) {
    auto dst = kv->first;
    std::string& dstInstanceName = dst.instanceName;
    fmi2ValueReference& dstValRef = dst.portReference;
    auto src = kv->second;
    std::string& srcInstanceName = src.instanceName;
    fmi2ValueReference& srcValRef = src.portReference;
    OutputDerivatives values = network->instances[srcInstanceName]->get_output(srcValRef);
    auto dstInstance = network->instances[dstInstanceName];
    Sample sample = { interval, dstInstance->set_input(dstValRef, values) };
    if (recordInputs) {
      Signal signal{ SignalType::Values, dst };
      results[signal].emplace_back(sample);
    }
  }
}

void fmi::record_outputs(CosimulationNetworkPtr network, double h, Results& results) {
  Interval interval = { network->t, network->t + h };
  for (auto outIt = network->monitoredOutputs.begin(); outIt != network->monitoredOutputs.end(); outIt++) {
    Output output = outIt->first;
    std::string& srcInstanceName = output.instanceName;
    fmi2ValueReference& srcValRef = output.portReference;
    auto srcInstance = network->instances[srcInstanceName];
    Sample sample = { interval, srcInstance->get_output(srcValRef) };
    Signal signal{ SignalType::Values, output };
    results[signal].emplace_back(sample);

    OptionalOutput optionalHermite = outIt->second;
    if (!optionalHermite.exists)
    {
      continue;
    }
    Output hermite = optionalHermite.output;
    fmi2ValueReference& hermiteValRef = hermite.portReference;
    Sample hermiteSample = { interval, srcInstance->get_hermite(srcValRef, hermiteValRef, h) };
    Signal hermiteSignal{ SignalType::Values, hermite };
    results[hermiteSignal].emplace_back(hermiteSample);
  }
}

void fmi::record_output_defects(Results& results, const MonitoredOutputs& monitoredOutputs) {
  for (auto it = monitoredOutputs.cbegin(); it != monitoredOutputs.cend(); it++) {
    Signal outputSignal = { Values, it->first };
    OptionalOutput hermiteOptional = it->second;
    if (!hermiteOptional.exists)
    {
      continue;
    }
    Signal hermiteSignal = { Values, hermiteOptional.output };
    Signal defectSignal = { Defect, it->first };

    Sample outputSample = results[outputSignal].back();
    Sample hermiteSample = results[hermiteSignal].back();

    OutputDerivatives defectDerivatives = signal::calculate_difference(hermiteSample.derivatives, outputSample.derivatives);
    Sample defectSample = { outputSample.interval, defectDerivatives };

    results[defectSignal].emplace_back(defectSample);
  }
}

void fmi::record_connection_defects(Results& results, const Connections& connections) {
  for (auto it = connections.cbegin(); it != connections.cend(); it++) {
    Signal outputSignal = { Values, it->second };
    Signal inputSignal = { Values, it->first };
    Signal defectSignal = { Defect, it->first };
    Sample outputSample = signal::shift_to_beginning(results[outputSignal].back());
    Sample inputSample = results[inputSignal].back();
    Derivatives diffDer = signal::calculate_difference(inputSample.derivatives, outputSample.derivatives);
    Sample defect = { inputSample.interval, diffDer };
    results[defectSignal].emplace_back(defect);
  }
}

bool Port::operator<(const Port &other)  const {
  return instanceName == other.instanceName ? portReference < other.portReference : instanceName < other.instanceName;
}

ConstantStepSize::ConstantStepSize(double h0, double tEnd) :_h0(h0),  _h(h0), _tEnd(tEnd), _t(0.) {}

double ConstantStepSize::initial() {
  _t += _h;
  return _h0;
}

double ConstantStepSize::next(double eps) {
  _t += _h;
  return _h;
}

bool ConstantStepSize::empty() {
  return _t >= _tEnd;
}

StepSequence::StepSequence(const Samples& samples) {
  for (auto it = samples.cbegin(); it != samples.cend(); it++) {
    double h = it->interval.end - it->interval.begin;
    _stepSequence.push(h);
  }
}

double StepSequence::initial() {
  return next(0.);
}

double StepSequence::next(double eps) {
  double hNext = _stepSequence.front();
  _stepSequence.pop();
  return hNext;
}

bool StepSequence::empty() {
  return _stepSequence.empty();
}

VariableStepSize::VariableStepSize(double h0, double tEnd, double tol)
: _h0(h0), _t(0.), _h(h0), _tEnd(tEnd), _I(std::log(h0)) ,_KP(0.13), _KI(1. / 15), _tol(tol), _phimax(2) { }

double VariableStepSize::initial() {
  _t += _h0;
  return _h0;
}

bool VariableStepSize::empty() {
  return _t >= _tEnd;
}

double VariableStepSize::next(double eps) {
  eps = eps > 1e-16 ?  eps : 1e-16;
  double e = std::log(_tol) - log(eps);
  double Ip = _I + _KI * e;
  double Hp = std::exp(Ip + _KP * e);
  _h = Hp > _phimax * _h ? _phimax * _h : Hp;
  _I = Ip + std::log(_h) - std::log(Hp);
  _t += _h;
  return _h;
};


double fmi::calculate_output_defects(const Results& results, const MonitoredOutputs& defectOutputs) {
  double eps = 0.;
  for (auto it = defectOutputs.cbegin(); it != defectOutputs.cend(); it++) {
    if (!it->second.exists)
    {
        continue;
    }
    Signal signal = { Defect, it->first };
    eps = std::max(eps, signal::rms(results.at(signal).back()));
  }
  return eps;
}

double fmi::calculate_connection_defects(const Results& results, const Connections& connections) {
  double eps = 0.;
  for (auto it = connections.cbegin(); it != connections.cend(); it++) {
    Signal signal = { Defect, it->first };
    eps = std::max(eps, signal::rms(results.at(signal).back()));
  }
  return eps;
}

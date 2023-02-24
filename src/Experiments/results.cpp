#include "results.h"
#include "signal.h"

#include <memory>
#include <numeric>
#include <cmath>
#include <fstream>

using namespace fmi;
using namespace fmi::results;
using namespace std;


void fmi::results::store_results(const Results& results, const char* csvFileName) {
  ofstream csvFile(csvFileName);
  const char* delimiter = ",";
  for (auto it = results.cbegin(); it != results.cend(); ++it) {
    const Signal& signal = it->first;
    const Port& port = signal.port;
    const Samples& samples = it->second;
    size_t maxDerivativeOrder = samples[0].derivatives.values.size() - 1;
    csvFile << "maxDerivativeOrder" << delimiter << port.instanceName.c_str() << delimiter << port.portReference << delimiter << maxDerivativeOrder << endl;
    csvFile << "steps" << delimiter << port.instanceName.c_str() << delimiter << port.portReference;
    for (auto itVal = samples.cbegin(); itVal != samples.cend(); ++itVal) {
      csvFile << delimiter << itVal->interval.end - itVal->interval.begin;
    }
    std::string signalTypeStr;
    switch (signal.type) {
    case SignalType::Values:
      signalTypeStr = "values";
      break;
    case SignalType::Defect:
      signalTypeStr = "defect";
      break;
    }
    csvFile << endl << signalTypeStr.c_str() << delimiter << port.instanceName.c_str() << delimiter << port.portReference;
    for (auto itVal = samples.cbegin(); itVal != samples.cend(); ++itVal) {
      const Sample shifted = signal::shift_to_beginning(*itVal);
      for (size_t derivative = 0; derivative <= maxDerivativeOrder; derivative++) {
        csvFile << delimiter << shifted.derivatives.values[derivative];
      }
    }
    csvFile << endl;
  }
}

double fmi::results::calculate_error(const Samples& signalX, const Samples& signalY) {
  return signal::rms(signal::calculate_difference(signalX, signalY));
}

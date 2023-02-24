#ifndef FMI_RESULTS_H
#define FMI_RESULTS_H

#include "fmi.h"


namespace fmi {
  namespace results {

    void store_results(const Results& results, const char* csvFileName);

    double calculate_error(const Samples& signalX, const Samples& signalY);
  }
}

#endif // FMI_RESULTS_H

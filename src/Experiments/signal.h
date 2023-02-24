#pragma once

#include "fmi.h"

namespace fmi {
  namespace signal {
    std::vector<fmi2Real> shift_derivatives(const std::vector<fmi2Real>& values, double h);

    Sample shift_to_beginning(const Sample& sample);

    double integrate(const Sample& sample);

    Derivatives multiply(const Derivatives& sampleX, const Derivatives& sampleY);
    
    Derivatives square(const Derivatives& sample);

    double rms(const Sample& sample);

    double rms(const Samples& samples);

    Samples calculate_difference(const Samples& signalX, const Samples& signalY);

    Derivatives calculate_difference(const Derivatives& signalX, const Derivatives& signalY);

    double eval(const Derivatives& signal, double t);
  }
}

#include "fmi.h"
#include "signal.h"

#include <memory>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <fstream>
#include <cassert>

using namespace std;
using namespace fmi;

// static vector<fmi2Real> facts({ 1, 1, 2, 6, 24 });
static vector<fmi2Real> facts({ 1 });

static inline fmi2Real factorial(size_t n) {
  if (facts.size() <= n) {
    for (size_t i = facts.size(); i <= n; i++) {
      facts.push_back(facts[i - 1] * i);
    }
  }
  return facts[n];
}

Derivatives fmi::signal::calculate_difference(const Derivatives& signalX, const Derivatives& signalY) {
  size_t minSize = min(signalX.values.size(), signalY.values.size());
  size_t maxSize = max(signalX.values.size(), signalY.values.size());
  Derivatives vals = { std::vector<fmi2Real>(maxSize, 0.), signalX.t };
  assert(fabs(signalX.t - signalY.t) < 1e-6);
  for (size_t i = 0; i < signalX.values.size(); i++)
  {
    vals.values[i] += signalX.values[i];
  }
  for (size_t i = 0; i < signalY.values.size(); i++)
  {
    vals.values[i] -= signalY.values[i];
  }
  return vals;
}

// difference = signalX - signalY
Samples fmi::signal::calculate_difference(const Samples& signalX, const Samples& signalY) {
  Samples difference;
  Samples::const_iterator itX = signalX.cbegin(), itY = signalY.cbegin();
  double t0 = max(itX->interval.begin, itY->interval.begin);
  while (itX != signalX.cend() && itY != signalY.cend()) {
    double tX = itX->interval.end, tY = itY->interval.end;
    double t = min(tX, tY);
    Interval interval = { t0, t };
    Sample sample = {
      interval,
      calculate_difference(itX->derivatives, itY->derivatives)
    };
    difference.emplace_back(sample);
    t0 = t;
    if (tX >= tY) {
      itY++;
    }
    if (tY >= tX) {
      itX++;
    }
  }
  return difference;
}

std::vector<fmi2Real> fmi::signal::shift_derivatives(const std::vector<fmi2Real>& values, double h) {
  std::vector<fmi2Real> shifted = values;
  size_t n = values.size() - 1;
  for (size_t m = 0; m <= n; m++) {
    shifted[m] = 0;
    fmi2Real hPow = 1;
    for (size_t k = 0; k <= n - m; k++) {
      shifted[m] += values[k + m] * hPow / factorial(k);
      hPow *= h;
    }
  }
  return shifted;
}

Sample fmi::signal::shift_to_beginning(const Sample& sample) {
  Sample shifted;
  shifted.derivatives.t = sample.interval.begin;
  shifted.interval = sample.interval;
  double h = shifted.derivatives.t - sample.derivatives.t;
  shifted.derivatives.values = shift_derivatives(sample.derivatives.values, h);
  return shifted;
}

double fmi::signal::integrate(const Sample& sample) {
  double h = sample.interval.end - sample.interval.begin;
  const Sample shifted = shift_to_beginning(sample);
  double integral = 0;
  const std::vector<fmi2Real>& values = shifted.derivatives.values;
  double hpown = 1;

  for (size_t i = 0; i < sample.derivatives.values.size(); i++) {
    hpown *= h / (i + 1);
    integral += values[i] * hpown;
  }
  return integral;
}

double fmi::signal::eval(const Derivatives& derivatives, double t) {
  double h = t - derivatives.t;
  double value = 0.;
  double hpown = 1;
  for (size_t i = 0; i < derivatives.values.size(); i++) {
    value += derivatives.values[i] * hpown;
    hpown *= h / (i + 1);
  }
  return value;
}

Derivatives fmi::signal::multiply(const Derivatives& sampleX, const Derivatives& sampleY) {
  size_t numDerivatives = sampleX.values.size() + sampleX.values.size() - 1;
  Derivatives multiplied = { std::vector<fmi2Real>(numDerivatives, 0.), sampleX.t };
  assert(sampleX.t == sampleY.t);
  for (size_t iX = 0; iX < sampleX.values.size(); iX++) {
    double coeffX = sampleX.values[iX] / factorial(iX);
    for (size_t iY = 0; iY < sampleY.values.size(); iY++) {
      double coeffY = sampleX.values[iY] / factorial(iY);
      size_t iM = iX + iY;
      multiplied.values[iM] += coeffX * coeffY * factorial(iM);
    }
  }
  return multiplied;
}

Derivatives fmi::signal::square(const Derivatives& sample) {
  return multiply(sample, sample);
}

double fmi::signal::rms(const Sample& sample) {
  Sample squared = sample;
  squared.derivatives = square(sample.derivatives);
  return sqrt(integrate(squared) / (sample.interval.end - sample.interval.begin));
}

double fmi::signal::rms(const Samples& samples) {
  double sumSquared = 0;
  for (auto it = samples.cbegin(); it != samples.cend(); it++) {
    Sample squared = *it;
    squared.derivatives = square(it->derivatives);
    sumSquared += integrate(squared);
  }
  return sqrt(sumSquared / (samples.back().interval.end - samples.front().interval.begin));
}

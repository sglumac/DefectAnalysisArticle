#include <fmi2TypesPlatform.h>
#include <fmi2Functions.h>
#include <fmi2FunctionTypes.h>

#include <stdio.h>
#include <string.h>
#include <math.h>

#define FMU_GUID "{3063ef2f-10cf-48fd-898e-833e2d41174a}"
#define MAX_DERIVATIVE_ORDER 3

#define ANALYTIC 0
#define EULER 1

typedef struct {
  fmi2Char* instanceName;
  const fmi2CallbackFunctions* environmentCallbacks;
  fmi2Real value;
  fmi2Real time, inputTime;
  fmi2Real u[MAX_DERIVATIVE_ORDER + 1], y[MAX_DERIVATIVE_ORDER + 1];
  fmi2Real yhMid;
  fmi2Integer factorials[MAX_DERIVATIVE_ORDER + 1];
  fmi2Real x, x0;
  fmi2Real eulerDx;
  fmi2Boolean initialized;
  fmi2Boolean outputCalculated;
  fmi2Integer solver;
} Component;

const char* fmi2GetTypesPlatform(void)
{
  return fmi2TypesPlatform;
}

const char* fmi2GetVersion(void)
{
  return fmi2Version;
}

fmi2Status  fmi2SetDebugLogging(fmi2Component c, fmi2Boolean b1, size_t n, const fmi2String categories[])
{
  return fmi2Warning;
}

static void calculate_analytic_output(Component* component) {
  fmi2Real h = component->time - component->inputTime;
  component->y[0] = 0.5 * component->x;
  fmi2Real dx = component->x;
  for (size_t i = 1; i <= MAX_DERIVATIVE_ORDER; i++) {
    dx = -dx;
    fmi2Real hpow = 1;
    fmi2Integer num = 1, fact = 1;
    for (size_t j = i - 1; j <= MAX_DERIVATIVE_ORDER; j++) {
      dx += component->u[j] * hpow / fact;
      hpow *= h;
      fact *= num++;
      component->y[i] = 0.5 * dx;
    }
  }
}

static void calculate_output(Component* component) {
  switch (component->solver)
  {
    case ANALYTIC:
    {
      calculate_analytic_output(component);
      break;
    }
    case EULER:
    {
      component->y[0] = 0.5 * component->x;
      component->y[1] = 0.5 * component->eulerDx;
      for (size_t i = 2; i <= MAX_DERIVATIVE_ORDER; i++) {
        component->y[i] = 0.;
      }
      break;
    }
  }
  component->outputCalculated = 1;
}

fmi2Status fmi2GetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i])
    {
    case 0:
      value[i] = component->u[0];
      break;
    case 1:
      if (!component->outputCalculated) {
        calculate_output(component);
      }
      value[i] = component->y[0];
      break;
    case 2:
      value[i] = component->x0;
      break;
    case 4:
      value[i] = component->yhMid;
      break;
    default:
      return fmi2Error;
    }
  }
  return fmi2OK;
}

fmi2Status fmi2GetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i]) {
    case 3:
      value[i] = component->solver;
    default:
      return fmi2Error;
    }
  }
  return fmi2OK;
}

fmi2Status fmi2GetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Boolean value[])
{
  if (nvr != 0)
  {
    return fmi2Error;
  }
  else
  {
    return fmi2OK;
  }
}

fmi2Status fmi2GetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2String value[])
{
  if (nvr != 0)
  {
    return fmi2Error;
  }
  else
  {
    return fmi2OK;
  }
}

fmi2Status fmi2SetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Real value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i])
    {
    case 0:
      component->u[0] = value[i];
      component->inputTime = component->time;
      component->outputCalculated = 0;
      break;
    case 2:
      if (!component->initialized) {
        component->x0 = value[i];
        component->x = component->x0;
        calculate_output(component);
      }
      else {
        return fmi2Error;
      }
    default:
      return fmi2Error;
    }
  }

  return fmi2OK;
}

fmi2Status fmi2SetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i]) {
    case 3:
      if (!component->initialized) {
        switch (value[i]) {
        case ANALYTIC:
        case EULER:
          component->solver = value[i];
          component->outputCalculated = 0;
          break;
        default:
          return fmi2Error;
        }        
      }
      else {
        return fmi2Error;
      }
    default:
      return fmi2Error;
    }
  }
  return fmi2OK;
}

fmi2Status fmi2SetBoolean(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Boolean value[])
{
  if (nvr != 0)
  {
    return fmi2Error;
  }
  else
  {
    return fmi2OK;
  }
}

fmi2Status fmi2SetString(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2String value[])
{
  if (nvr != 0)
  {
    return fmi2Error;
  }
  else
  {
    return fmi2OK;
  }
}

fmi2Component fmi2Instantiate(fmi2String instanceName, fmi2Type fmuType, fmi2String fmuGuid, fmi2String fmuResourceLocation,
  const fmi2CallbackFunctions* functions, fmi2Boolean visible, fmi2Boolean loggingOn)
{
  if (functions == NULL) // unable to log anything, this model expects that logging is enabled
  {
    return NULL;
  }
  else if (instanceName == NULL || !strlen(instanceName)) // check if instance name is valid, should not be modified
  {
    functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error", "Invalid instance name!");
    return NULL;
  }
  else if (strcmp(fmuGuid, FMU_GUID) != 0) // check if GUID is valid, should not be modified
  {
    functions->logger(functions->componentEnvironment, instanceName, fmi2Error, "error", "Invalid GUID of a class!");
    return NULL;
  }
  else
  {
    // create model instance
    Component *component = (Component *)malloc(sizeof(Component));
    //initialize time
    component->time = 0.0;
    component->inputTime = component->time;
    // model variables initialization
    component->value = 3.0;

    component->instanceName = strdup(instanceName);
    component->environmentCallbacks = functions;

    component->environmentCallbacks->logger(component->environmentCallbacks->componentEnvironment, component->instanceName, fmi2OK, "", "Hello from fmi2Instantiate!");

    component->solver = ANALYTIC;

    component->u[0] = 0;
    component->x0 = 3;
    component->x = component->x0;
    component->eulerDx = -component->x0;

    for (size_t i = 1; i <= MAX_DERIVATIVE_ORDER; i++) {
      component->u[i] = 0;
    }
    calculate_output(component);

    component->initialized = fmi2False;

    component->factorials[0] = 1;
    for (size_t i = 1; i <= MAX_DERIVATIVE_ORDER; i++) {
      component->factorials[i] = component->factorials[i - 1] * i;
    }

    return component;
  }
}

fmi2Status fmi2SetupExperiment(fmi2Component c, fmi2Boolean tolDefined, fmi2Real tol, fmi2Real statTime, fmi2Boolean stopDefined, fmi2Real stopTime)
{
  Component* s = (Component*)c;

  s->environmentCallbacks->logger(s->environmentCallbacks->componentEnvironment, s->instanceName, fmi2OK, "", "Hello from fmi2SetupExperiment!");

  return fmi2OK;
}

fmi2Status fmi2EnterInitializationMode(fmi2Component c)
{
  if (c == NULL) 
  {
    return fmi2Error;
  }
  else
  {
    Component* s = (Component*)c;

    s->environmentCallbacks->logger(s->environmentCallbacks->componentEnvironment, s->instanceName, fmi2OK, "", "Hello from fmi2EnterInitializationMode!");

    return fmi2OK;
  }
}

fmi2Status fmi2ExitInitializationMode(fmi2Component c)
{
  if (c == NULL)
  {
    return fmi2Error;
  }
  else
  {
    Component* s = (Component*)c;

    s->environmentCallbacks->logger(s->environmentCallbacks->componentEnvironment, s->instanceName, fmi2OK, "", "Hello from fmi2ExitInitializationMode!");

    s->initialized = fmi2True;

    return fmi2OK;
  }
}

static void shift_input_derivatives(Component* component)
{
  fmi2Real shifted[MAX_DERIVATIVE_ORDER + 1];
  size_t n = MAX_DERIVATIVE_ORDER;
  fmi2Real h = component->time - component->inputTime;
  for (size_t m = 0; m <= n; m++) {
    shifted[m] = 0;
    fmi2Real hPow = 1;
    for (size_t k = 0; k <= n - m; k++) {
      shifted[m] += component->u[k + m] * hPow / component->factorials[k];
      hPow *= h;
    }
  }
  for (size_t m = 0; m <= n; m++) {
    component->u[m] = shifted[m];
  }
  component->inputTime = component->time;
}

static void analytic_step(Component* component, double h)
{
  double pxp[MAX_DERIVATIVE_ORDER + 1];
  pxp[MAX_DERIVATIVE_ORDER] = component->u[MAX_DERIVATIVE_ORDER] / component->factorials[MAX_DERIVATIVE_ORDER];
  for (int i = MAX_DERIVATIVE_ORDER - 1; i >= 0; i--)
  {
    pxp[i] = component->u[i] / component->factorials[i] - (i + 1) * pxp[i + 1];
  }
  double ch = component->x - pxp[0];
  double xh = ch * exp(-h);
  component->x = xh;
  double hp = 1;
  for (size_t i = 0; i <= MAX_DERIVATIVE_ORDER; i++)
  {
    component->x += pxp[i] * hp;
    hp *= h;
  }
}

fmi2Status fmi2DoStep(fmi2Component c, fmi2Real t, fmi2Real h,  fmi2Boolean noSet)
{
  Component* component = c;

  if (component->solver == ANALYTIC) {
    analytic_step(component, h / 2);
    component->time += h / 2;
    component->yhMid = 0.5 * component->x;
    shift_input_derivatives(component);
    analytic_step(component, h / 2);
    component->time += h / 2;
  }
  else if (component->solver == EULER) {
    component->eulerDx = -component->x + component->u[0];
    component->yhMid = 0.5 * component->x + 0.5 * component->eulerDx * h / 2;
    component->x += component->eulerDx * h;
    component->time += h;
  }
  component->outputCalculated = 0;

  shift_input_derivatives(component);

  return fmi2OK;
}

fmi2Status fmi2Reset(fmi2Component c)
{
  return fmi2OK;
}

fmi2Status fmi2Terminate(fmi2Component c)
{
  Component* s = (Component*)c;

  s->environmentCallbacks->logger(s->environmentCallbacks->componentEnvironment, s->instanceName, fmi2OK, "", "Hello from fmi2Terminate!");

  return fmi2OK;
}

void fmi2FreeInstance(fmi2Component c)
{
  Component* s = (Component*)c;

  s->environmentCallbacks->logger(s->environmentCallbacks->componentEnvironment, s->instanceName, fmi2OK, "", "Hello from fmi2FreeInstance!");
  
  free(s->instanceName);
  free(s);
}

fmi2Status fmi2GetFMUstate(fmi2Component c, fmi2FMUstate* s)
{
  return fmi2Error;
}

fmi2Status fmi2SetFMUstate(fmi2Component c, fmi2FMUstate s)
{
  return fmi2Error;
}

fmi2Status fmi2FreeFMUstate(fmi2Component c, fmi2FMUstate* s)
{
  return fmi2Error;
}

fmi2Status fmi2SerializedFMUstateSize(fmi2Component c, fmi2FMUstate s, size_t* n)
{
  return fmi2Error;
}

fmi2Status fmi2SerializeFMUstate(fmi2Component c, fmi2FMUstate s, fmi2Byte xs[], size_t y)
{
  return fmi2Error;
}

fmi2Status fmi2DeSerializeFMUstate(fmi2Component c, const fmi2Byte xs[], size_t n, fmi2FMUstate* s)
{
  return fmi2Error;
}

fmi2Status fmi2GetDirectionalDerivative(fmi2Component c, const fmi2ValueReference vrs[], size_t nvr,
  const fmi2ValueReference vrs2[], size_t n,
  const fmi2Real someReals1[], fmi2Real someReals2[])
{
  return fmi2Error;
}

fmi2Status fmi2SetRealInputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], const fmi2Real value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i])
    {
    case 0:
      component->u[order[i]] = value[i];
      component->inputTime = component->time;
      component->outputCalculated = 0;
      break;
    default:
      return fmi2Error;
    }
  }

  return fmi2OK;
}

fmi2Status fmi2GetRealOutputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], fmi2Real value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i])
    {
    case 1:
      if (!component->outputCalculated) {
        calculate_output(component);
      }
      value[i] = component->y[order[i]];
      break;
    default:
      return fmi2Error;
    }
  }
  return fmi2OK;
}

fmi2Status fmi2CancelStep(fmi2Component c)
{
  return fmi2Error;
}

fmi2Status fmi2GetStatus(fmi2Component c, const fmi2StatusKind k, fmi2Status* s)
{
  return fmi2Error;
}

fmi2Status fmi2GetRealStatus(fmi2Component c, const fmi2StatusKind k, fmi2Real* r)
{
  return fmi2Error;
}

fmi2Status fmi2GetIntegerStatus(fmi2Component c, const fmi2StatusKind k, fmi2Integer* i)
{
  return fmi2Error;
}

fmi2Status fmi2GetBooleanStatus(fmi2Component c, const fmi2StatusKind k, fmi2Boolean* b)
{
  return fmi2Error;
}
fmi2Status fmi2GetStringStatus(fmi2Component c, const fmi2StatusKind k, fmi2String* str)
{
  return fmi2Error;
}

#include <fmi2TypesPlatform.h>
#include <fmi2Functions.h>
#include <fmi2FunctionTypes.h>

#include <stdio.h>
#include <string.h>
#include <math.h>

#define FMU_GUID "{799f52bd-d05b-4e92-8325-81d2bf6f7b09}"
#define MAX_DERIVATIVE_ORDER 10

typedef struct {
  fmi2Char* instanceName;
  const fmi2CallbackFunctions* environmentCallbacks;
  fmi2Real value;
  fmi2Real time;
  fmi2Real y1[MAX_DERIVATIVE_ORDER + 1], y2[MAX_DERIVATIVE_ORDER + 1];
  fmi2Real y1h[MAX_DERIVATIVE_ORDER + 1], y2h[MAX_DERIVATIVE_ORDER + 1];
  fmi2Real x10, x20;
  fmi2Real x1, x2;
  fmi2Boolean initialized;
} Component;

static void calculate_output(Component* component, fmi2Real* y1, fmi2Real* y2) {
  fmi2Real dx1 = component->x1, dx2 = component->x2;
  y1[0] = 0.5 * component->x1;
  y2[0] = 0.5 * component->x2;
  for (size_t i = 1; i <= MAX_DERIVATIVE_ORDER; i++) {
    dx1 = -dx1 + y2[i - 1];
    dx2 = -dx2 + y1[i - 1];
    y1[i] = 0.5 * dx1;
    y2[i] = 0.5 * dx2;
  }
}

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

fmi2Status fmi2GetReal(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Real value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i])
    {
    case 0:
      value[i] = component->y1[0];
      break;
    case 1:
      value[i] = component->y2[0];
      break;
    case 2:
      value[i] = component->x10;
      break;
    case 3:
      value[i] = component->x20;
      break;
    case 4:
      value[i] = component->y1h[0];
      break;
    case 5:
      value[i] = component->y1h[0];
      break;
    default:
      return fmi2Error;
    }
  }
  return fmi2OK;
}

fmi2Status fmi2GetInteger(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, fmi2Integer value[])
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
    case 2:
      if (!component->initialized) {
        component->x10 = value[i];
        component->x1 = component->x10;
        calculate_output(component, component->y1, component->y2);
        calculate_output(component, component->y1h, component->y2h);
      }
      else {
        return fmi2Error;
      }
    case 3:
      if (!component->initialized) {
        component->x20 = value[i];
        component->x2 = component->x20;
        calculate_output(component, component->y1, component->y2);
        calculate_output(component, component->y1h, component->y2h);
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
  if (nvr != 0)
  {
    return fmi2Error;
  }
  else
  {
    return fmi2OK;
  }
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
    // model variables initialization
    component->value = 3.0;

    component->instanceName = strdup(instanceName);
    component->environmentCallbacks = functions;

    component->environmentCallbacks->logger(component->environmentCallbacks->componentEnvironment, component->instanceName, fmi2OK, "", "Hello from fmi2Instantiate!");

    component->x10 = 3;
    component->x1 = component->x10;

    component->x20 = 3;
    component->x2 = component->x20;

    calculate_output(component, component->y1, component->y2);
    calculate_output(component, component->y1h, component->y2h);

    component->initialized = fmi2False;

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

static void step(Component* component, double h)
{
  fmi2Real C11 = 0.5 * component->x1 - 0.5 * component->x2;
  fmi2Real C12 = 0.5 * component->x1 + 0.5 * component->x2;
  fmi2Real C21 = 0.5 * component->x2 - 0.5 * component->x1;
  fmi2Real C22 = 0.5 * component->x2 + 0.5 * component->x1;

  component->x1 = C11 * exp(-1.5 * h) + C12 * exp(-0.5 * h);
  component->x2 = C21 * exp(-1.5 * h) + C22 * exp(-0.5 * h);
}


fmi2Status fmi2DoStep(fmi2Component c, fmi2Real t, fmi2Real h, fmi2Boolean noSet)
{
  Component* component = c;

  step(component, h / 2);
  calculate_output(component, component->y1h, component->y2h);

  step(component, h / 2);
  calculate_output(component, component->y1, component->y2);

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

fmi2Status fmi2SetRealInputDerivatives(fmi2Component c, const fmi2ValueReference vrs[], size_t n, const fmi2Integer sd[], const fmi2Real fdg[])
{
  return fmi2Error;
}

fmi2Status fmi2GetRealOutputDerivatives(fmi2Component c, const fmi2ValueReference vr[], size_t nvr, const fmi2Integer order[], fmi2Real value[])
{
  Component* component = c;
  size_t i;

  for (i = 0; i < nvr; i++)
  {
    switch (vr[i])
    {
    case 0:
      value[i] = component->y1[order[i]];
      break;
    case 1:
      value[i] = component->y2[order[i]];
      break;
    case 4:
      value[i] = component->y1h[order[i]];
      break;
    case 5:
      value[i] = component->y2h[order[i]];
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

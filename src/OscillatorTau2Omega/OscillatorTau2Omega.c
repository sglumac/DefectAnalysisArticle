/*
 * MIT License
 * 
 * Copyright (c) 2017 Slaven Glumac
 * 
 * Single mass force-to-displacement oscillator.
 */
#include <fmi2Functions.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <cvode/cvode_dense.h>
#include <cvode/cvode_diag.h>
#include <sundials/sundials_types.h>

#define MAX_DERIVATIVE_ORDER 2
#define NUMBER_OF_REALS 9
#define NUMBER_OF_INTEGERS 0
#define NUMBER_OF_BOOLEANS 0
#define NUMBER_OF_STRINGS 0
#define NUMBER_OF_VARS 9

const fmi2ValueReference ivrs[] = {0, 1, 2, 3, 4, 5, 1, 6, 7, 8};

#define vr_tauOther 0
#define _tauOther(d) r(vr_tauOther,d)
#define vr_omegaThis 1
#define _omegaThis(d) r(vr_omegaThis,d)
#define _J r(2,0)
#define _c r(3,0)
#define _d r(4,0)
#define _phiThis0 r(5,0)
#define _omegaThis0 r(6,0)
#define _solver r(7, 0)
#define vr_omegaThisHMid 8
#define _omegaThisHMid r(vr_omegaThisHMid,0)
#define _phiThis(d) r(9, d)

#define USES_CVODE (_solver < 0.5)

#define NUMBER_OF_STATES 2
#define _phiThisF y[0]
#define _omegaThisF y[1]
#define _dphiThisF dy[0]
#define _domegaThisF dy[1]
#define Jac(i,j) DENSE_ELEM(J,i,j)

#define _ycvode _internal.ycvode
#define _cvode _internal.cvode
#define _yE _internal.yE
struct Internal
{
    N_Vector ycvode;
    void* cvode;
    double yE[NUMBER_OF_STATES];
};

#include <template.h>

static void f(fmi2Component component, double* y, double* dy, double tauOther)
{
    _dphiThisF = _omegaThisF;
    _domegaThisF = -_c / _J * _phiThisF - _d / _J * _omegaThisF + 1. / _J * tauOther;
}

static int fcvode(realtype t, N_Vector ycvode, N_Vector dycvode, void *user_data)
{
    fmi2Component component = user_data;
    double y[NUMBER_OF_STATES], dy[NUMBER_OF_STATES];
    for (size_t i = 0; i < NUMBER_OF_STATES; i++)
    {
        y[i] = NV_Ith_S(ycvode,i);
    }

    fmi2Real tauOther = interp(component, vr_tauOther, t);
    f(component, y, dy, tauOther);

    for (size_t i = 0; i < NUMBER_OF_STATES; i++)
    {
        NV_Ith_S(dycvode,i) = dy[i];
    }
    return CV_SUCCESS;
}

static void update_cvode_state_derivatives(fmi2Component component) {
    for (size_t i = 1; i <= MAX_DERIVATIVE_ORDER; i++)
    {
        _phiThis(i) = _omegaThis(i-1);
        _omegaThis(i) = -_c / _J * _phiThis(i-1) - _d / _J * _omegaThis(i-1) + 1. / _J * _tauOther(i-1);
    }
}

static void update_euler_state_derivatives(fmi2Component component) {
    _phiThis(1) = _omegaThis(0);
    _omegaThis(1) = -_c / _J * _phiThis(0) - _d / _J * _omegaThis(0) + 1. / _J * _tauOther(0);
    _phiThis(2) = 0.;
    _omegaThis(2) = 0.;
}

static int Jacobian(long int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    fmi2Component component = user_data;

    Jac(0,0) = 0.;
    Jac(0,1) = 1.;
    Jac(1,0) = -_c / _J;
    Jac(1,1) = -_d / _J;

    return CV_SUCCESS;
}

void InstantiateInternal(fmi2Component component)
{
    _cvode = CVodeCreate(CV_BDF, CV_NEWTON);
    _ycvode = N_VNew_Serial(NUMBER_OF_STATES);
    _solver = 0;
}

void FreeInternal(fmi2Component component)
{
    N_VDestroy_Serial(_ycvode);
    CVodeFree(&_cvode);
}

fmi2Status InitializeIntegrator(fmi2Component component)
{
    realtype reltol = 0.;
    realtype abstol = 1e-8;
    log(fmi2OK, "Hello from InitializeIntegrator!");
    if (CVodeInit(_cvode, fcvode, _t, _ycvode) != CV_SUCCESS)
    {
        return fmi2Error;
    }
    if (CVodeSetUserData(_cvode, component) != CV_SUCCESS)
    {
        return fmi2Error;
    }
    if (CVodeSStolerances(_cvode, reltol, abstol) != CV_SUCCESS)
    {
      return fmi2Error;
    }
    if (CVDense(_cvode, NUMBER_OF_STATES) != CV_SUCCESS)
    {
        return fmi2Error;
    }
    if (CVDlsSetDenseJacFn(_cvode, Jacobian) != CV_SUCCESS)
    {
        return fmi2Error;
    }
    return fmi2OK;
}

void StartInitialization(fmi2Component component)
{
    _tauOther(0) = 0.;
    _J = 1.;
    _c = 1.;
    _d = 1.;
    _phiThis0 = 0.1;
    _omegaThis0 = 0.1;
}

fmi2Status FinishInitialization(fmi2Component component)
{
    double y[NUMBER_OF_STATES];
    _phiThisF = _phiThis0;
    _phiThis(0) = _phiThis0;
    _omegaThisF = _omegaThis0;
	_omegaThis(0) = _omegaThis0;
    for (size_t i = 0; i < NUMBER_OF_STATES; i++)
    {
        NV_Ith_S(_ycvode,i) = y[i];
        _yE[i] = y[i];
    }
    return InitializeIntegrator(component);
}

fmi2Status StateUpdate(fmi2Component component, fmi2Real h)
{
    double y[NUMBER_OF_STATES];
    if (USES_CVODE)
    {
        realtype tReached;
        if (CVode(_cvode, _t + h / 2, _ycvode, &tReached, CV_NORMAL) != CV_SUCCESS)
        {
            log(fmi2Error, "The integration failed!");
            return fmi2Error;
        }
        shift_input_derivatives(component, vr_tauOther, _t + h / 2);
        _omegaThisHMid = NV_Ith_S(_ycvode,1);
        if (CVode(_cvode, _t + h, _ycvode, &tReached, CV_NORMAL) != CV_SUCCESS)
        {
            log(fmi2Error, "The integration failed!");
            return fmi2Error;
        }
        for (size_t i = 0; i < NUMBER_OF_STATES; i++)
        {
                y[i] = NV_Ith_S(_ycvode,i);
        }
        shift_input_derivatives(component, vr_tauOther, _t + h);
        _phiThis(0) = _phiThisF;
        _omegaThis(0) = _omegaThisF;
        update_cvode_state_derivatives(component);
    }
    else
    {
        update_euler_state_derivatives(component);
        double dyE[NUMBER_OF_STATES];
        double tauOther = interp(component, vr_tauOther, _t);
        f(component, _yE, dyE, tauOther);
        euler(_yE, dyE, NUMBER_OF_STATES, h / 2);
        _omegaThisHMid = _yE[1];
        euler(_yE, dyE, NUMBER_OF_STATES, h / 2);
        for (size_t i = 0; i < NUMBER_OF_STATES; i++)
        {
            y[i] = _yE[i];
        }
        shift_input_derivatives(component, vr_tauOther, _t + h);
        _phiThis(0) = _phiThisF;
        _omegaThis(0) = _omegaThisF;
    }
    return fmi2OK;
}

fmi2Status OutputUpdate(fmi2Component component)
{
	return fmi2OK;
}

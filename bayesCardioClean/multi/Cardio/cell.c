#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "geometry.h"
#include "FEM.h"
#include "../Amazigh/printFile.h"
#include "../Amazigh/mathVector.h"
#include "../Amazigh/solvers.h"
#include "../Amazigh/plot2D.h"
#include "../Amazigh/supportC.h"


/*
There are a total of 8 entries in the algebraic variable array.
There are a total of 3 entries in each of the rate and state variable arrays.
There are a total of 22 entries in the constant variable array.
*/
/*
* VOI is time in component environment (ms).
* STATES[0] is u in component membrane (dimensionless).
* CONSTANTS[0] is Cm in component membrane (uF_per_cm2).
* ALGEBRAIC[0] is Vm in component membrane (mV).
* CONSTANTS[1] is V_0 in component membrane (mV).
* CONSTANTS[2] is V_fi in component membrane (mV).
* ALGEBRAIC[2] is J_fi in component fast_inward_current (per_ms).
* ALGEBRAIC[4] is J_so in component slow_outward_current (per_ms).
* ALGEBRAIC[6] is J_si in component slow_inward_current (per_ms).
* ALGEBRAIC[7] is Istim in component stimulus_protocol (per_ms).
* ALGEBRAIC[1] is p in component p (dimensionless).
* CONSTANTS[3] is u_c in component p (dimensionless).
* ALGEBRAIC[3] is q in component q (dimensionless).
* CONSTANTS[4] is u_v in component q (dimensionless).
* CONSTANTS[21] is tau_d in component fast_inward_current (ms).
* CONSTANTS[5] is g_fi_max in component fast_inward_current (mS_per_cm2).
* STATES[1] is v in component fast_inward_current_v_gate (dimensionless).
* ALGEBRAIC[5] is tau_v_minus in component fast_inward_current_v_gate (ms).
* CONSTANTS[6] is tau_v1_minus in component fast_inward_current_v_gate (ms).
* CONSTANTS[7] is tau_v2_minus in component fast_inward_current_v_gate (ms).
* CONSTANTS[8] is tau_v_plus in component fast_inward_current_v_gate (ms).
* CONSTANTS[9] is tau_0 in component slow_outward_current (ms).
* CONSTANTS[10] is tau_r in component slow_outward_current (ms).
* CONSTANTS[11] is tau_si in component slow_inward_current (ms).
* CONSTANTS[12] is u_csi in component slow_inward_current (dimensionless).
* CONSTANTS[13] is k in component slow_inward_current (dimensionless).
* STATES[2] is w in component slow_inward_current_w_gate (dimensionless).
* CONSTANTS[14] is tau_w_minus in component slow_inward_current_w_gate (ms).
* CONSTANTS[15] is tau_w_plus in component slow_inward_current_w_gate (ms).
* CONSTANTS[16] is IstimStart in component stimulus_protocol (ms).
* CONSTANTS[17] is IstimEnd in component stimulus_protocol (ms).
* CONSTANTS[18] is IstimAmplitude in component stimulus_protocol (per_ms).
* CONSTANTS[19] is IstimPeriod in component stimulus_protocol (ms).
* CONSTANTS[20] is IstimPulseDuration in component stimulus_protocol (ms).
* RATES[0] is d/dt u in component membrane (dimensionless).
* RATES[1] is d/dt v in component fast_inward_current_v_gate (dimensionless).
* RATES[2] is d/dt w in component slow_inward_current_w_gate (dimensionless).
*/
void
initConsts_FemtomKarma(double* CONSTANTS, double* RATES, double *STATES)
{
	STATES[0] = 0;
	CONSTANTS[0] = 1;
	CONSTANTS[1] = -85;
	CONSTANTS[2] = 15;
	CONSTANTS[3] = 0.13;
	CONSTANTS[4] = 0.04;
	CONSTANTS[5] = 4;
	STATES[1] = 1;
	CONSTANTS[6] = 1250;
	CONSTANTS[7] = 19.6;
	CONSTANTS[8] = 3.33;
	CONSTANTS[9] = 12.5;
	CONSTANTS[10] = 33.33;
	CONSTANTS[11] = 29;
	CONSTANTS[12] = 0.85;
	CONSTANTS[13] = 10;
	STATES[2] = 1;
	CONSTANTS[14] = 41;
	CONSTANTS[15] = 870;
	CONSTANTS[16] = 10;
	CONSTANTS[17] = 50000;
	CONSTANTS[18] = -0.2;
	CONSTANTS[19] = 1000;
	CONSTANTS[20] = 1;
	CONSTANTS[21] = CONSTANTS[0] / CONSTANTS[5];
}
void
computeRates_FemtomKarma(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
	ALGEBRAIC[1] = (STATES[0]<CONSTANTS[3] ? 0.00000 : 1.00000);
	RATES[2] = ((1.00000 - ALGEBRAIC[1])*(1.00000 - STATES[2])) / CONSTANTS[14] - (ALGEBRAIC[1] * STATES[2]) / CONSTANTS[15];
	ALGEBRAIC[3] = (STATES[0]<CONSTANTS[4] ? 0.00000 : 1.00000);
	ALGEBRAIC[5] = ALGEBRAIC[3] * CONSTANTS[6] + (1.00000 - ALGEBRAIC[3])*CONSTANTS[7];
	RATES[1] = ((1.00000 - ALGEBRAIC[1])*(1.00000 - STATES[1])) / ALGEBRAIC[5] - (ALGEBRAIC[1] * STATES[1]) / CONSTANTS[8];
	ALGEBRAIC[2] = (-STATES[1] * ALGEBRAIC[1] * (1.00000 - STATES[0])*(STATES[0] - CONSTANTS[3])) / CONSTANTS[21];
	ALGEBRAIC[4] = (STATES[0] * (1.00000 - ALGEBRAIC[1])) / CONSTANTS[9] + ALGEBRAIC[1] / CONSTANTS[10];
	ALGEBRAIC[6] = (-STATES[2] * (1.00000 + tanh(CONSTANTS[13] * (STATES[0] - CONSTANTS[12])))) / (2.00000*CONSTANTS[11]);
	ALGEBRAIC[7] = (VOI >= CONSTANTS[16] && VOI <= CONSTANTS[17] && (VOI - CONSTANTS[16]) - floor((VOI - CONSTANTS[16]) / CONSTANTS[19])*CONSTANTS[19] <= CONSTANTS[20] ? CONSTANTS[18] : 0.00000);
	//ALGEBRAIC[7] = 0;
	RATES[0] = -(ALGEBRAIC[2] + ALGEBRAIC[4] + ALGEBRAIC[6] + ALGEBRAIC[7]);
}
void
computeVariables_FemtomKarma(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
	ALGEBRAIC[1] = (STATES[0]<CONSTANTS[3] ? 0.00000 : 1.00000);
	ALGEBRAIC[3] = (STATES[0]<CONSTANTS[4] ? 0.00000 : 1.00000);
	ALGEBRAIC[5] = ALGEBRAIC[3] * CONSTANTS[6] + (1.00000 - ALGEBRAIC[3])*CONSTANTS[7];
	ALGEBRAIC[2] = (-STATES[1] * ALGEBRAIC[1] * (1.00000 - STATES[0])*(STATES[0] - CONSTANTS[3])) / CONSTANTS[21];
	ALGEBRAIC[4] = (STATES[0] * (1.00000 - ALGEBRAIC[1])) / CONSTANTS[9] + ALGEBRAIC[1] / CONSTANTS[10];
	ALGEBRAIC[6] = (-STATES[2] * (1.00000 + tanh(CONSTANTS[13] * (STATES[0] - CONSTANTS[12])))) / (2.00000*CONSTANTS[11]);
	ALGEBRAIC[7] = (VOI >= CONSTANTS[16] && VOI <= CONSTANTS[17] && (VOI - CONSTANTS[16]) - floor((VOI - CONSTANTS[16]) / CONSTANTS[19])*CONSTANTS[19] <= CONSTANTS[20] ? CONSTANTS[18] : 0.00000);
	ALGEBRAIC[0] = CONSTANTS[1] + STATES[0] * (CONSTANTS[2] - CONSTANTS[1]);
}

void run_FemtomKarma()
{
	double* ALGEBRAIC = newDoubleV(8);
	double* CONSTANTS = newDoubleV(22);
	double* RATES = newDoubleV(3);
	double* STATES = newDoubleV(3);

	initConsts_FemtomKarma(CONSTANTS, RATES, STATES);

	double dt = 0.125;
	double tmax = 2000;
	double tinit = 0;
	int steps = (int)(tmax - tinit) / dt;

	double** Values = newDouble(steps, 3);
	
	char** Names[3];
	Names[0] = "U";
	Names[1] = "V";
	Names[2] = "W";

	

	double t = 0;
	for (int i = 0; i < steps; i++)
	{
		

		double Istim = 0;
		//computeVariables(t, CONSTANTS, RATES, STATES, ALGEBRAIC);
		computeRates_FemtomKarma(t, CONSTANTS, RATES, STATES, ALGEBRAIC);

		eulerODE(RATES, STATES, 3, dt);

		for (int j = 0; j < 3; j++)
		{
			Values[i][j] = STATES[j];
		}
		
		t += dt;
	
	}

	plot2DComparisonNamesNew(Values, 3, steps, "Results", "Graph.html", Names,"Amplitude");

	freeMatrix(Values, steps);
	
	free(CONSTANTS);
	free(RATES);
	free(STATES);
	free(ALGEBRAIC);
	
}

void run_FemtomKarmaRK2()
{
	double* ALGEBRAIC = newDoubleV(8);
	double* CONSTANTS = newDoubleV(22);
	double* RATES = newDoubleV(3);
	double* STATES = newDoubleV(3);

	initConsts_FemtomKarma(CONSTANTS, RATES, STATES);

	double dt = 0.125;
	double tmax = 2000;
	double tinit = 0;
	int steps = (int)(tmax - tinit) / dt;

	double** Values = newDouble(steps, 3);

	char** Names[3];
	Names[0] = "U";
	Names[1] = "V";
	Names[2] = "W";



	double t = 0;
	for (int i = 0; i < steps; i++)
	{
		
		double* k1 = newDoubleV(3);
		computeRates_FemtomKarma(t, CONSTANTS, k1, STATES, ALGEBRAIC);

		double* k2 = newDoubleV(3);
		double* tempk1 = multConstantV(k1, dt / 2,3);
		double* statesk1 = addVector(STATES, tempk1,3);
		computeRates_FemtomKarma(t+dt/2, CONSTANTS, k2, statesk1, ALGEBRAIC);


		for (int j = 0; j < 3; j++)
		{
			STATES[j] += k2[j] * dt;
			Values[i][j] = STATES[j];
		}

		free(k1);
		free(k2);
		free(tempk1);
		free(statesk1);

		t += dt;
	}

	plot2DComparisonNamesNew(Values, 3, steps, "Results", "Graph.html", Names, "Amplitude");

	freeMatrix(Values, steps);

	free(CONSTANTS);
	free(RATES);
	free(STATES);
	free(ALGEBRAIC);

}

/*
There are a total of 58 entries in the algebraic variable array.
There are a total of 12 entries in each of the rate and state variable arrays.
There are a total of 60 entries in the constant variable array.
*/
/*
* VOI is time in component environment (ms).
* STATES[0] is V in component membrane (mV).
* CONSTANTS[0] is R in component membrane (gas_constant_units).
* CONSTANTS[1] is T in component membrane (kelvin).
* CONSTANTS[2] is F in component membrane (faradays_constant_units).
* ALGEBRAIC[53] is dV_dt in component membrane (mV_per_ms).
* CONSTANTS[3] is Cm in component membrane (uF_per_mm2).
* ALGEBRAIC[3] is I_st in component membrane (uA_per_mm2).
* ALGEBRAIC[14] is i_Na in component fast_sodium_current (uA_per_mm2).
* ALGEBRAIC[26] is i_Ca_L in component L_type_Ca_channel (uA_per_mm2).
* ALGEBRAIC[29] is i_K in component time_dependent_potassium_current (uA_per_mm2).
* ALGEBRAIC[52] is i_NaCa in component Na_Ca_exchanger (uA_per_mm2).
* ALGEBRAIC[34] is i_K1 in component time_independent_potassium_current (uA_per_mm2).
* ALGEBRAIC[37] is i_Kp in component plateau_potassium_current (uA_per_mm2).
* ALGEBRAIC[38] is i_p_Ca in component sarcolemmal_calcium_pump (uA_per_mm2).
* ALGEBRAIC[40] is i_Na_b in component sodium_background_current (uA_per_mm2).
* ALGEBRAIC[42] is i_Ca_b in component calcium_background_current (uA_per_mm2).
* ALGEBRAIC[44] is i_NaK in component sodium_potassium_pump (uA_per_mm2).
* ALGEBRAIC[51] is i_ns_Ca in component non_specific_calcium_activated_current (uA_per_mm2).
* CONSTANTS[4] is stimPeriod in component membrane (dimensionless).
* CONSTANTS[5] is stimDuration in component membrane (dimensionless).
* CONSTANTS[6] is stimCurrent in component membrane (dimensionless).
* ALGEBRAIC[10] is E_Na in component fast_sodium_current (mV).
* CONSTANTS[7] is g_Na in component fast_sodium_current (mS_per_mm2).
* STATES[1] is Nai in component ionic_concentrations (mM).
* CONSTANTS[8] is Nao in component ionic_concentrations (mM).
* STATES[2] is m in component fast_sodium_current_m_gate (dimensionless).
* STATES[3] is h in component fast_sodium_current_h_gate (dimensionless).
* STATES[4] is j in component fast_sodium_current_j_gate (dimensionless).
* ALGEBRAIC[0] is alpha_m in component fast_sodium_current_m_gate (per_ms).
* ALGEBRAIC[7] is beta_m in component fast_sodium_current_m_gate (per_ms).
* ALGEBRAIC[1] is alpha_h in component fast_sodium_current_h_gate (per_ms).
* ALGEBRAIC[8] is beta_h in component fast_sodium_current_h_gate (per_ms).
* ALGEBRAIC[2] is alpha_j in component fast_sodium_current_j_gate (per_ms).
* ALGEBRAIC[9] is beta_j in component fast_sodium_current_j_gate (per_ms).
* ALGEBRAIC[23] is i_CaCa in component L_type_Ca_channel (uA_per_mm2).
* ALGEBRAIC[25] is i_CaK in component L_type_Ca_channel (uA_per_mm2).
* ALGEBRAIC[24] is i_CaNa in component L_type_Ca_channel (uA_per_mm2).
* CONSTANTS[9] is gamma_Nai in component L_type_Ca_channel (dimensionless).
* CONSTANTS[10] is gamma_Nao in component L_type_Ca_channel (dimensionless).
* CONSTANTS[11] is gamma_Ki in component L_type_Ca_channel (dimensionless).
* CONSTANTS[12] is gamma_Ko in component L_type_Ca_channel (dimensionless).
* ALGEBRAIC[17] is I_CaCa in component L_type_Ca_channel (uA_per_mm2).
* ALGEBRAIC[21] is I_CaK in component L_type_Ca_channel (uA_per_mm2).
* ALGEBRAIC[20] is I_CaNa in component L_type_Ca_channel (uA_per_mm2).
* CONSTANTS[13] is P_Ca in component L_type_Ca_channel (mm_per_ms).
* CONSTANTS[14] is P_Na in component L_type_Ca_channel (mm_per_ms).
* CONSTANTS[15] is P_K in component L_type_Ca_channel (mm_per_ms).
* CONSTANTS[16] is gamma_Cai in component L_type_Ca_channel (dimensionless).
* CONSTANTS[17] is gamma_Cao in component L_type_Ca_channel (dimensionless).
* STATES[5] is Cai in component ionic_concentrations (mM).
* CONSTANTS[18] is Cao in component ionic_concentrations (mM).
* CONSTANTS[19] is Ko in component ionic_concentrations (mM).
* STATES[6] is Ki in component ionic_concentrations (mM).
* STATES[7] is d in component L_type_Ca_channel_d_gate (dimensionless).
* STATES[8] is f in component L_type_Ca_channel_f_gate (dimensionless).
* ALGEBRAIC[22] is f_Ca in component L_type_Ca_channel_f_Ca_gate (dimensionless).
* ALGEBRAIC[15] is alpha_d in component L_type_Ca_channel_d_gate (per_ms).
* ALGEBRAIC[18] is beta_d in component L_type_Ca_channel_d_gate (per_ms).
* ALGEBRAIC[4] is d_infinity in component L_type_Ca_channel_d_gate (dimensionless).
* ALGEBRAIC[11] is tau_d in component L_type_Ca_channel_d_gate (ms).
* ALGEBRAIC[16] is alpha_f in component L_type_Ca_channel_f_gate (per_ms).
* ALGEBRAIC[19] is beta_f in component L_type_Ca_channel_f_gate (per_ms).
* ALGEBRAIC[5] is f_infinity in component L_type_Ca_channel_f_gate (dimensionless).
* ALGEBRAIC[12] is tau_f in component L_type_Ca_channel_f_gate (ms).
* CONSTANTS[20] is Km_Ca in component L_type_Ca_channel_f_Ca_gate (mM).
* CONSTANTS[21] is g_K_max in component time_dependent_potassium_current (mS_per_mm2).
* CONSTANTS[54] is g_K in component time_dependent_potassium_current (mS_per_mm2).
* ALGEBRAIC[27] is E_K in component time_dependent_potassium_current (mV).
* CONSTANTS[22] is PR_NaK in component time_dependent_potassium_current (dimensionless).
* STATES[9] is X in component time_dependent_potassium_current_X_gate (dimensionless).
* ALGEBRAIC[28] is Xi in component time_dependent_potassium_current_Xi_gate (dimensionless).
* ALGEBRAIC[6] is alpha_X in component time_dependent_potassium_current_X_gate (per_ms).
* ALGEBRAIC[13] is beta_X in component time_dependent_potassium_current_X_gate (per_ms).
* ALGEBRAIC[30] is E_K1 in component time_independent_potassium_current (mV).
* CONSTANTS[23] is g_K1_max in component time_independent_potassium_current (mS_per_mm2).
* CONSTANTS[55] is g_K1 in component time_independent_potassium_current (mS_per_mm2).
* ALGEBRAIC[33] is K1_infinity in component time_independent_potassium_current_K1_gate (dimensionless).
* ALGEBRAIC[31] is alpha_K1 in component time_independent_potassium_current_K1_gate (per_ms).
* ALGEBRAIC[32] is beta_K1 in component time_independent_potassium_current_K1_gate (per_ms).
* ALGEBRAIC[35] is E_Kp in component plateau_potassium_current (mV).
* CONSTANTS[24] is g_Kp in component plateau_potassium_current (mS_per_mm2).
* ALGEBRAIC[36] is Kp in component plateau_potassium_current (dimensionless).
* CONSTANTS[25] is K_mpCa in component sarcolemmal_calcium_pump (mM).
* CONSTANTS[26] is I_pCa in component sarcolemmal_calcium_pump (uA_per_mm2).
* CONSTANTS[27] is g_Nab in component sodium_background_current (mS_per_mm2).
* ALGEBRAIC[39] is E_NaN in component sodium_background_current (mV).
* CONSTANTS[28] is g_Cab in component calcium_background_current (mS_per_mm2).
* ALGEBRAIC[41] is E_CaN in component calcium_background_current (mV).
* CONSTANTS[29] is I_NaK in component sodium_potassium_pump (uA_per_mm2).
* ALGEBRAIC[43] is f_NaK in component sodium_potassium_pump (dimensionless).
* CONSTANTS[30] is K_mNai in component sodium_potassium_pump (mM).
* CONSTANTS[31] is K_mKo in component sodium_potassium_pump (mM).
* CONSTANTS[56] is sigma in component sodium_potassium_pump (dimensionless).
* ALGEBRAIC[48] is i_ns_Na in component non_specific_calcium_activated_current (uA_per_mm2).
* ALGEBRAIC[50] is i_ns_K in component non_specific_calcium_activated_current (uA_per_mm2).
* CONSTANTS[32] is P_ns_Ca in component non_specific_calcium_activated_current (mm_per_ms).
* ALGEBRAIC[47] is I_ns_Na in component non_specific_calcium_activated_current (uA_per_mm2).
* ALGEBRAIC[49] is I_ns_K in component non_specific_calcium_activated_current (uA_per_mm2).
* CONSTANTS[33] is K_m_ns_Ca in component non_specific_calcium_activated_current (mM).
* ALGEBRAIC[46] is Vns in component non_specific_calcium_activated_current (mV).
* ALGEBRAIC[45] is EnsCa in component non_specific_calcium_activated_current (mV).
* CONSTANTS[34] is K_NaCa in component Na_Ca_exchanger (uA_per_mm2).
* CONSTANTS[35] is K_mNa in component Na_Ca_exchanger (mM).
* CONSTANTS[36] is K_mCa in component Na_Ca_exchanger (mM).
* CONSTANTS[37] is K_sat in component Na_Ca_exchanger (dimensionless).
* CONSTANTS[38] is eta in component Na_Ca_exchanger (dimensionless).
* ALGEBRAIC[54] is i_rel in component calcium_fluxes_in_the_SR (mM_per_ms).
* ALGEBRAIC[55] is i_up in component calcium_fluxes_in_the_SR (mM_per_ms).
* ALGEBRAIC[56] is i_leak in component calcium_fluxes_in_the_SR (mM_per_ms).
* ALGEBRAIC[57] is i_tr in component calcium_fluxes_in_the_SR (mM_per_ms).
* CONSTANTS[59] is G_rel in component calcium_fluxes_in_the_SR (per_ms).
* CONSTANTS[57] is G_rel_peak in component calcium_fluxes_in_the_SR (per_ms).
* CONSTANTS[39] is G_rel_max in component calcium_fluxes_in_the_SR (per_ms).
* CONSTANTS[40] is tau_on in component calcium_fluxes_in_the_SR (ms).
* CONSTANTS[41] is tau_off in component calcium_fluxes_in_the_SR (ms).
* CONSTANTS[42] is t_CICR in component calcium_fluxes_in_the_SR (ms).
* CONSTANTS[43] is tau_tr in component calcium_fluxes_in_the_SR (ms).
* CONSTANTS[44] is K_mrel in component calcium_fluxes_in_the_SR (mM).
* CONSTANTS[45] is K_mup in component calcium_fluxes_in_the_SR (mM).
* CONSTANTS[58] is K_leak in component calcium_fluxes_in_the_SR (per_ms).
* CONSTANTS[46] is I_up in component calcium_fluxes_in_the_SR (mM_per_ms).
* CONSTANTS[47] is Ca_NSR_max in component calcium_fluxes_in_the_SR (mM).
* CONSTANTS[48] is delta_Ca_i2 in component calcium_fluxes_in_the_SR (mM).
* CONSTANTS[49] is delta_Ca_ith in component calcium_fluxes_in_the_SR (mM).
* STATES[10] is Ca_JSR in component ionic_concentrations (mM).
* STATES[11] is Ca_NSR in component ionic_concentrations (mM).
* CONSTANTS[50] is Am in component ionic_concentrations (per_mm).
* CONSTANTS[51] is V_myo in component ionic_concentrations (dimensionless).
* CONSTANTS[52] is V_JSR in component ionic_concentrations (dimensionless).
* CONSTANTS[53] is V_NSR in component ionic_concentrations (dimensionless).
* RATES[0] is d/dt V in component membrane (mV).
* RATES[2] is d/dt m in component fast_sodium_current_m_gate (dimensionless).
* RATES[3] is d/dt h in component fast_sodium_current_h_gate (dimensionless).
* RATES[4] is d/dt j in component fast_sodium_current_j_gate (dimensionless).
* RATES[7] is d/dt d in component L_type_Ca_channel_d_gate (dimensionless).
* RATES[8] is d/dt f in component L_type_Ca_channel_f_gate (dimensionless).
* RATES[9] is d/dt X in component time_dependent_potassium_current_X_gate (dimensionless).
* RATES[1] is d/dt Nai in component ionic_concentrations (mM).
* RATES[5] is d/dt Cai in component ionic_concentrations (mM).
* RATES[6] is d/dt Ki in component ionic_concentrations (mM).
* RATES[10] is d/dt Ca_JSR in component ionic_concentrations (mM).
* RATES[11] is d/dt Ca_NSR in component ionic_concentrations (mM).
*/
void
initConsts_LuoRudy(double* CONSTANTS, double* RATES, double *STATES)
{
	STATES[0] = -84.624;
	CONSTANTS[0] = 8.3145e3;
	CONSTANTS[1] = 310.0;
	CONSTANTS[2] = 96845.0;
	CONSTANTS[3] = 0.01;
	CONSTANTS[4] = 1e3;
	CONSTANTS[5] = 0.5;
	CONSTANTS[6] = 0.5;
	CONSTANTS[7] = 0.16;
	STATES[1] = 10.0;
	CONSTANTS[8] = 140.0;
	STATES[2] = 0.0;
	STATES[3] = 1.0;
	STATES[4] = 1.0;
	CONSTANTS[9] = 0.75;
	CONSTANTS[10] = 0.75;
	CONSTANTS[11] = 0.75;
	CONSTANTS[12] = 0.75;
	CONSTANTS[13] = 5.4e-6;
	CONSTANTS[14] = 6.75e-9;
	CONSTANTS[15] = 1.93e-9;
	CONSTANTS[16] = 1.0;
	CONSTANTS[17] = 0.34;
	STATES[5] = 0.12e-3;
	CONSTANTS[18] = 1.8;
	CONSTANTS[19] = 5.4;
	STATES[6] = 145.0;
	STATES[7] = 0.0;
	STATES[8] = 1.0;
	CONSTANTS[20] = 0.6e-3;
	CONSTANTS[21] = 2.82e-3;
	CONSTANTS[22] = 0.01833;
	STATES[9] = 0.0;
	CONSTANTS[23] = 7.5e-3;
	CONSTANTS[24] = 1.83e-4;
	CONSTANTS[25] = 0.5e-3;
	CONSTANTS[26] = 1.15e-2;
	CONSTANTS[27] = 1.41e-5;
	CONSTANTS[28] = 3.016e-5;
	CONSTANTS[29] = 1.5e-2;
	CONSTANTS[30] = 10.0;
	CONSTANTS[31] = 1.5;
	CONSTANTS[32] = 1.75e-9;
	CONSTANTS[33] = 1.2e-3;
	CONSTANTS[34] = 20.0;
	CONSTANTS[35] = 87.5;
	CONSTANTS[36] = 1.38;
	CONSTANTS[37] = 0.1;
	CONSTANTS[38] = 0.35;
	CONSTANTS[39] = 60.0;
	CONSTANTS[40] = 2.0;
	CONSTANTS[41] = 2.0;
	CONSTANTS[42] = 0.0;
	CONSTANTS[43] = 180.0;
	CONSTANTS[44] = 0.8e-3;
	CONSTANTS[45] = 0.92e-3;
	CONSTANTS[46] = 0.005;
	CONSTANTS[47] = 15.0;
	CONSTANTS[48] = 0.0;
	CONSTANTS[49] = 0.18e-3;
	STATES[10] = 1.8;
	STATES[11] = 1.8;
	CONSTANTS[50] = 200;
	CONSTANTS[51] = 0.68;
	CONSTANTS[52] = 0.0048;
	CONSTANTS[53] = 0.0552;
	CONSTANTS[54] = CONSTANTS[21] * pow((CONSTANTS[19] / 5.40000), 1.0 / 2);
	CONSTANTS[55] = CONSTANTS[23] * pow((CONSTANTS[19] / 5.40000), 1.0 / 2);
	CONSTANTS[56] = (1.00000 / 7.00000)*(exp(CONSTANTS[8] / 67.3000) - 1.00000);
	CONSTANTS[57] = (CONSTANTS[48]<CONSTANTS[49] ? 0.00000 : CONSTANTS[39]);
	CONSTANTS[58] = CONSTANTS[46] / CONSTANTS[47];
	CONSTANTS[59] = CONSTANTS[57] * ((CONSTANTS[48] - CONSTANTS[49]) / ((CONSTANTS[44] + CONSTANTS[48]) - CONSTANTS[49]))*(1.00000 - exp(-(CONSTANTS[42] / CONSTANTS[40])))*exp(-(CONSTANTS[42] / CONSTANTS[41]));
}
void
computeRates_LuoRudy(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
	ALGEBRAIC[0] = (0.320000*(STATES[0] + 47.1300)) / (1.00000 - exp(-0.100000*(STATES[0] + 47.1300)));
	ALGEBRAIC[7] = 0.0800000*exp(-STATES[0] / 11.0000);
	RATES[2] = ALGEBRAIC[0] * (1.00000 - STATES[2]) - ALGEBRAIC[7] * STATES[2];
	ALGEBRAIC[1] = (STATES[0]<-40.0000 ? 0.135000*exp((80.0000 + STATES[0]) / -6.80000) : 0.00000);
	ALGEBRAIC[8] = (STATES[0]<-40.0000 ? 3.56000*exp(0.0790000*STATES[0]) + 310000.*exp(0.350000*STATES[0]) : 1.00000 / (0.130000*(1.00000 + exp((STATES[0] + 10.6600) / -11.1000))));
	RATES[3] = ALGEBRAIC[1] * (1.00000 - STATES[3]) - ALGEBRAIC[8] * STATES[3];
	ALGEBRAIC[2] = (STATES[0]<-40.0000 ? (-127140.*exp(0.244400*STATES[0]) - 3.47400e-05*exp(-0.0439100*STATES[0]))*((STATES[0] + 37.7800) / (1.00000 + exp(0.311000*(STATES[0] + 79.2300)))) : 0.00000);
	ALGEBRAIC[9] = (STATES[0]<-40.0000 ? (0.121200*exp(-0.0105200*STATES[0])) / (1.00000 + exp(-0.137800*(STATES[0] + 40.1400))) : (0.300000*exp(-2.53500e-07*STATES[0])) / (1.00000 + exp(-0.100000*(STATES[0] + 32.0000))));
	RATES[4] = ALGEBRAIC[2] * (1.00000 - STATES[4]) - ALGEBRAIC[9] * STATES[4];
	ALGEBRAIC[6] = (7.19000e-05*(STATES[0] + 30.0000)) / (1.00000 - exp(-0.148000*(STATES[0] + 30.0000)));
	ALGEBRAIC[13] = (0.000131000*(STATES[0] + 30.0000)) / (-1.00000 + exp(0.0687000*(STATES[0] + 30.0000)));
	RATES[9] = ALGEBRAIC[6] * (1.00000 - STATES[9]) - ALGEBRAIC[13] * STATES[9];
	ALGEBRAIC[4] = 1.00000 / (1.00000 + exp(-((STATES[0] + 10.0000) / 6.24000)));
	ALGEBRAIC[11] = ALGEBRAIC[4] * ((1.00000 - exp(-((STATES[0] + 10.0000) / 6.24000))) / (0.0350000*(STATES[0] + 10.0000)));
	ALGEBRAIC[15] = ALGEBRAIC[4] / ALGEBRAIC[11];
	ALGEBRAIC[18] = (1.00000 - ALGEBRAIC[4]) / ALGEBRAIC[11];
	RATES[7] = ALGEBRAIC[15] * (1.00000 - STATES[7]) - ALGEBRAIC[18] * STATES[7];
	ALGEBRAIC[5] = 1.00000 / (1.00000 + exp((STATES[0] + 35.0600) / 8.60000)) + 0.600000 / (1.00000 + exp((50.0000 - STATES[0]) / 20.0000));
	ALGEBRAIC[12] = 1.00000 / (0.0197000*exp(-pow(0.0337000*(STATES[0] + 10.0000), 2.00000)) + 0.0200000);
	ALGEBRAIC[16] = ALGEBRAIC[5] / ALGEBRAIC[12];
	ALGEBRAIC[19] = (1.00000 - ALGEBRAIC[5]) / ALGEBRAIC[12];
	RATES[8] = ALGEBRAIC[16] * (1.00000 - STATES[8]) - ALGEBRAIC[19] * STATES[8];
	ALGEBRAIC[27] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log((CONSTANTS[19] + CONSTANTS[22] * CONSTANTS[8]) / (STATES[6] + CONSTANTS[22] * STATES[1]));
	ALGEBRAIC[28] = 1.00000 / (1.00000 + exp((STATES[0] - 56.2600) / 32.1000));
	ALGEBRAIC[29] = CONSTANTS[54] * pow(STATES[9], 2.00000)*ALGEBRAIC[28] * (STATES[0] - ALGEBRAIC[27]);
	ALGEBRAIC[30] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log(CONSTANTS[19] / STATES[6]);
	ALGEBRAIC[31] = 1.02000 / (1.00000 + exp(0.238500*((STATES[0] - ALGEBRAIC[30]) - 59.2150)));
	ALGEBRAIC[32] = (0.491240*exp(0.0803200*((STATES[0] + 5.47600) - ALGEBRAIC[30])) + exp(0.0617500*(STATES[0] - (ALGEBRAIC[30] + 594.310)))) / (1.00000 + exp(-0.514300*((STATES[0] - ALGEBRAIC[30]) + 4.75300)));
	ALGEBRAIC[33] = ALGEBRAIC[31] / (ALGEBRAIC[31] + ALGEBRAIC[32]);
	ALGEBRAIC[34] = CONSTANTS[55] * ALGEBRAIC[33] * (STATES[0] - ALGEBRAIC[30]);
	ALGEBRAIC[35] = ALGEBRAIC[30];
	ALGEBRAIC[36] = 1.00000 / (1.00000 + exp((7.48800 - STATES[0]) / 5.98000));
	ALGEBRAIC[37] = CONSTANTS[24] * ALGEBRAIC[36] * (STATES[0] - ALGEBRAIC[35]);
	ALGEBRAIC[43] = 1.00000 / ((1.00000 + 0.124500*exp(-0.100000*((STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])))) + 0.0365000*CONSTANTS[56] * exp(-((STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1]))));
	ALGEBRAIC[44] = CONSTANTS[29] * ALGEBRAIC[43] * (1.00000 / (1.00000 + pow(CONSTANTS[30] / STATES[1], 1.50000)))*(CONSTANTS[19] / (CONSTANTS[19] + CONSTANTS[31]));
	ALGEBRAIC[21] = CONSTANTS[15] * pow(1.00000, 2.00000)*((STATES[0] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[11] * STATES[6] * exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[12] * CONSTANTS[19]) / (exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[22] = 1.00000 / (1.00000 + pow(STATES[5] / CONSTANTS[20], 2.00000));
	ALGEBRAIC[25] = STATES[7] * STATES[8] * ALGEBRAIC[22] * ALGEBRAIC[21];
	ALGEBRAIC[45] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log((CONSTANTS[19] + CONSTANTS[8]) / (STATES[6] + STATES[1]));
	ALGEBRAIC[46] = STATES[0] - ALGEBRAIC[45];
	ALGEBRAIC[49] = CONSTANTS[32] * pow(1.00000, 2.00000)*((ALGEBRAIC[46] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[11] * STATES[6] * exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[12] * CONSTANTS[19]) / (exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[50] = ALGEBRAIC[49] * (1.00000 / (1.00000 + pow(CONSTANTS[33] / STATES[5], 3.00000)));
	RATES[6] = -(ALGEBRAIC[25] + ALGEBRAIC[29] + ALGEBRAIC[34] + ALGEBRAIC[37] + ALGEBRAIC[50] + -(ALGEBRAIC[44] * 2.00000))*(CONSTANTS[50] / (CONSTANTS[51] * CONSTANTS[2]));
	ALGEBRAIC[10] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log(CONSTANTS[8] / STATES[1]);
	ALGEBRAIC[14] = CONSTANTS[7] * pow(STATES[2], 3.00000)*STATES[3] * STATES[4] * (STATES[0] - ALGEBRAIC[10]);
	ALGEBRAIC[52] = CONSTANTS[34] * (1.00000 / (pow(CONSTANTS[35], 3.00000) + pow(CONSTANTS[8], 3.00000)))*(1.00000 / (CONSTANTS[36] + CONSTANTS[18]))*(1.00000 / (1.00000 + CONSTANTS[37] * exp((CONSTANTS[38] - 1.00000)*STATES[0] * (CONSTANTS[2] / (CONSTANTS[0] * CONSTANTS[1])))))*(exp(CONSTANTS[38] * STATES[0] * (CONSTANTS[2] / (CONSTANTS[0] * CONSTANTS[1])))*pow(STATES[1], 3.00000)*CONSTANTS[18] - exp((CONSTANTS[38] - 1.00000)*STATES[0] * (CONSTANTS[2] / (CONSTANTS[0] * CONSTANTS[1])))*pow(CONSTANTS[8], 3.00000)*STATES[5]);
	ALGEBRAIC[39] = ALGEBRAIC[10];
	ALGEBRAIC[40] = CONSTANTS[27] * (STATES[0] - ALGEBRAIC[39]);
	ALGEBRAIC[20] = CONSTANTS[14] * pow(1.00000, 2.00000)*((STATES[0] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[9] * STATES[1] * exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[10] * CONSTANTS[8]) / (exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[24] = STATES[7] * STATES[8] * ALGEBRAIC[22] * ALGEBRAIC[20];
	ALGEBRAIC[47] = CONSTANTS[32] * pow(1.00000, 2.00000)*((ALGEBRAIC[46] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[9] * STATES[1] * exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[10] * CONSTANTS[8]) / (exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[48] = ALGEBRAIC[47] * (1.00000 / (1.00000 + pow(CONSTANTS[33] / STATES[5], 3.00000)));
	RATES[1] = -(ALGEBRAIC[14] + ALGEBRAIC[24] + ALGEBRAIC[40] + ALGEBRAIC[48] + ALGEBRAIC[52] * 3.00000 + ALGEBRAIC[44] * 3.00000)*(CONSTANTS[50] / (CONSTANTS[51] * CONSTANTS[2]));
	ALGEBRAIC[3] = ((int)(VOI) % (int)(CONSTANTS[4])<CONSTANTS[5] ? CONSTANTS[6] : 0.00000);
	ALGEBRAIC[17] = CONSTANTS[13] * pow(2.00000, 2.00000)*((STATES[0] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[16] * STATES[5] * exp((2.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[17] * CONSTANTS[18]) / (exp((2.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[23] = STATES[7] * STATES[8] * ALGEBRAIC[22] * ALGEBRAIC[17];
	ALGEBRAIC[26] = ALGEBRAIC[23] + ALGEBRAIC[25] + ALGEBRAIC[24];
	ALGEBRAIC[38] = CONSTANTS[26] * (STATES[5] / (CONSTANTS[25] + STATES[5]));
	ALGEBRAIC[41] = ((CONSTANTS[0] * CONSTANTS[1]) / (2.00000*CONSTANTS[2]))*log(CONSTANTS[18] / STATES[5]);
	ALGEBRAIC[42] = CONSTANTS[28] * (STATES[0] - ALGEBRAIC[41]);
	ALGEBRAIC[51] = ALGEBRAIC[48] + ALGEBRAIC[50];
	ALGEBRAIC[53] = (ALGEBRAIC[3] - (ALGEBRAIC[14] + ALGEBRAIC[26] + ALGEBRAIC[29] + ALGEBRAIC[34] + ALGEBRAIC[37] + ALGEBRAIC[52] + ALGEBRAIC[38] + ALGEBRAIC[40] + ALGEBRAIC[42] + ALGEBRAIC[44] + ALGEBRAIC[51])) / CONSTANTS[3];
	RATES[0] = ALGEBRAIC[53];
	ALGEBRAIC[54] = CONSTANTS[59] * (STATES[10] - STATES[5]);
	ALGEBRAIC[55] = CONSTANTS[46] * (STATES[5] / (STATES[5] + CONSTANTS[45]));
	ALGEBRAIC[56] = CONSTANTS[58] * STATES[11];
	RATES[5] = -((ALGEBRAIC[23] + ALGEBRAIC[38] + ALGEBRAIC[42]) - ALGEBRAIC[52])*(CONSTANTS[50] / (2.00000*CONSTANTS[51] * CONSTANTS[2])) + ALGEBRAIC[54] * (CONSTANTS[52] / CONSTANTS[51]) + (ALGEBRAIC[56] - ALGEBRAIC[55])*(CONSTANTS[53] / CONSTANTS[51]);
	ALGEBRAIC[57] = (STATES[11] - STATES[10]) / CONSTANTS[43];
	RATES[10] = -(ALGEBRAIC[54] - ALGEBRAIC[57] * (CONSTANTS[53] / CONSTANTS[52]));
	RATES[11] = -((ALGEBRAIC[56] + ALGEBRAIC[57]) - ALGEBRAIC[55]);
}
void
computeVariables_LuoRudy(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
	ALGEBRAIC[0] = (0.320000*(STATES[0] + 47.1300)) / (1.00000 - exp(-0.100000*(STATES[0] + 47.1300)));
	ALGEBRAIC[7] = 0.0800000*exp(-STATES[0] / 11.0000);
	ALGEBRAIC[1] = (STATES[0]<-40.0000 ? 0.135000*exp((80.0000 + STATES[0]) / -6.80000) : 0.00000);
	ALGEBRAIC[8] = (STATES[0]<-40.0000 ? 3.56000*exp(0.0790000*STATES[0]) + 310000.*exp(0.350000*STATES[0]) : 1.00000 / (0.130000*(1.00000 + exp((STATES[0] + 10.6600) / -11.1000))));
	ALGEBRAIC[2] = (STATES[0]<-40.0000 ? (-127140.*exp(0.244400*STATES[0]) - 3.47400e-05*exp(-0.0439100*STATES[0]))*((STATES[0] + 37.7800) / (1.00000 + exp(0.311000*(STATES[0] + 79.2300)))) : 0.00000);
	ALGEBRAIC[9] = (STATES[0]<-40.0000 ? (0.121200*exp(-0.0105200*STATES[0])) / (1.00000 + exp(-0.137800*(STATES[0] + 40.1400))) : (0.300000*exp(-2.53500e-07*STATES[0])) / (1.00000 + exp(-0.100000*(STATES[0] + 32.0000))));
	ALGEBRAIC[6] = (7.19000e-05*(STATES[0] + 30.0000)) / (1.00000 - exp(-0.148000*(STATES[0] + 30.0000)));
	ALGEBRAIC[13] = (0.000131000*(STATES[0] + 30.0000)) / (-1.00000 + exp(0.0687000*(STATES[0] + 30.0000)));
	ALGEBRAIC[4] = 1.00000 / (1.00000 + exp(-((STATES[0] + 10.0000) / 6.24000)));
	ALGEBRAIC[11] = ALGEBRAIC[4] * ((1.00000 - exp(-((STATES[0] + 10.0000) / 6.24000))) / (0.0350000*(STATES[0] + 10.0000)));
	ALGEBRAIC[15] = ALGEBRAIC[4] / ALGEBRAIC[11];
	ALGEBRAIC[18] = (1.00000 - ALGEBRAIC[4]) / ALGEBRAIC[11];
	ALGEBRAIC[5] = 1.00000 / (1.00000 + exp((STATES[0] + 35.0600) / 8.60000)) + 0.600000 / (1.00000 + exp((50.0000 - STATES[0]) / 20.0000));
	ALGEBRAIC[12] = 1.00000 / (0.0197000*exp(-pow(0.0337000*(STATES[0] + 10.0000), 2.00000)) + 0.0200000);
	ALGEBRAIC[16] = ALGEBRAIC[5] / ALGEBRAIC[12];
	ALGEBRAIC[19] = (1.00000 - ALGEBRAIC[5]) / ALGEBRAIC[12];
	ALGEBRAIC[27] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log((CONSTANTS[19] + CONSTANTS[22] * CONSTANTS[8]) / (STATES[6] + CONSTANTS[22] * STATES[1]));
	ALGEBRAIC[28] = 1.00000 / (1.00000 + exp((STATES[0] - 56.2600) / 32.1000));
	ALGEBRAIC[29] = CONSTANTS[54] * pow(STATES[9], 2.00000)*ALGEBRAIC[28] * (STATES[0] - ALGEBRAIC[27]);
	ALGEBRAIC[30] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log(CONSTANTS[19] / STATES[6]);
	ALGEBRAIC[31] = 1.02000 / (1.00000 + exp(0.238500*((STATES[0] - ALGEBRAIC[30]) - 59.2150)));
	ALGEBRAIC[32] = (0.491240*exp(0.0803200*((STATES[0] + 5.47600) - ALGEBRAIC[30])) + exp(0.0617500*(STATES[0] - (ALGEBRAIC[30] + 594.310)))) / (1.00000 + exp(-0.514300*((STATES[0] - ALGEBRAIC[30]) + 4.75300)));
	ALGEBRAIC[33] = ALGEBRAIC[31] / (ALGEBRAIC[31] + ALGEBRAIC[32]);
	ALGEBRAIC[34] = CONSTANTS[55] * ALGEBRAIC[33] * (STATES[0] - ALGEBRAIC[30]);
	ALGEBRAIC[35] = ALGEBRAIC[30];
	ALGEBRAIC[36] = 1.00000 / (1.00000 + exp((7.48800 - STATES[0]) / 5.98000));
	ALGEBRAIC[37] = CONSTANTS[24] * ALGEBRAIC[36] * (STATES[0] - ALGEBRAIC[35]);
	ALGEBRAIC[43] = 1.00000 / ((1.00000 + 0.124500*exp(-0.100000*((STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])))) + 0.0365000*CONSTANTS[56] * exp(-((STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1]))));
	ALGEBRAIC[44] = CONSTANTS[29] * ALGEBRAIC[43] * (1.00000 / (1.00000 + pow(CONSTANTS[30] / STATES[1], 1.50000)))*(CONSTANTS[19] / (CONSTANTS[19] + CONSTANTS[31]));
	ALGEBRAIC[21] = CONSTANTS[15] * pow(1.00000, 2.00000)*((STATES[0] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[11] * STATES[6] * exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[12] * CONSTANTS[19]) / (exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[22] = 1.00000 / (1.00000 + pow(STATES[5] / CONSTANTS[20], 2.00000));
	ALGEBRAIC[25] = STATES[7] * STATES[8] * ALGEBRAIC[22] * ALGEBRAIC[21];
	ALGEBRAIC[45] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log((CONSTANTS[19] + CONSTANTS[8]) / (STATES[6] + STATES[1]));
	ALGEBRAIC[46] = STATES[0] - ALGEBRAIC[45];
	ALGEBRAIC[49] = CONSTANTS[32] * pow(1.00000, 2.00000)*((ALGEBRAIC[46] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[11] * STATES[6] * exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[12] * CONSTANTS[19]) / (exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[50] = ALGEBRAIC[49] * (1.00000 / (1.00000 + pow(CONSTANTS[33] / STATES[5], 3.00000)));
	ALGEBRAIC[10] = ((CONSTANTS[0] * CONSTANTS[1]) / CONSTANTS[2])*log(CONSTANTS[8] / STATES[1]);
	ALGEBRAIC[14] = CONSTANTS[7] * pow(STATES[2], 3.00000)*STATES[3] * STATES[4] * (STATES[0] - ALGEBRAIC[10]);
	ALGEBRAIC[52] = CONSTANTS[34] * (1.00000 / (pow(CONSTANTS[35], 3.00000) + pow(CONSTANTS[8], 3.00000)))*(1.00000 / (CONSTANTS[36] + CONSTANTS[18]))*(1.00000 / (1.00000 + CONSTANTS[37] * exp((CONSTANTS[38] - 1.00000)*STATES[0] * (CONSTANTS[2] / (CONSTANTS[0] * CONSTANTS[1])))))*(exp(CONSTANTS[38] * STATES[0] * (CONSTANTS[2] / (CONSTANTS[0] * CONSTANTS[1])))*pow(STATES[1], 3.00000)*CONSTANTS[18] - exp((CONSTANTS[38] - 1.00000)*STATES[0] * (CONSTANTS[2] / (CONSTANTS[0] * CONSTANTS[1])))*pow(CONSTANTS[8], 3.00000)*STATES[5]);
	ALGEBRAIC[39] = ALGEBRAIC[10];
	ALGEBRAIC[40] = CONSTANTS[27] * (STATES[0] - ALGEBRAIC[39]);
	ALGEBRAIC[20] = CONSTANTS[14] * pow(1.00000, 2.00000)*((STATES[0] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[9] * STATES[1] * exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[10] * CONSTANTS[8]) / (exp((1.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[24] = STATES[7] * STATES[8] * ALGEBRAIC[22] * ALGEBRAIC[20];
	ALGEBRAIC[47] = CONSTANTS[32] * pow(1.00000, 2.00000)*((ALGEBRAIC[46] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[9] * STATES[1] * exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[10] * CONSTANTS[8]) / (exp((1.00000*ALGEBRAIC[46] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[48] = ALGEBRAIC[47] * (1.00000 / (1.00000 + pow(CONSTANTS[33] / STATES[5], 3.00000)));
	ALGEBRAIC[3] = ((int)(VOI) % (int)(CONSTANTS[4])<CONSTANTS[5] ? CONSTANTS[6] : 0.00000);
	ALGEBRAIC[17] = CONSTANTS[13] * pow(2.00000, 2.00000)*((STATES[0] * pow(CONSTANTS[2], 2.00000)) / (CONSTANTS[0] * CONSTANTS[1]))*((CONSTANTS[16] * STATES[5] * exp((2.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - CONSTANTS[17] * CONSTANTS[18]) / (exp((2.00000*STATES[0] * CONSTANTS[2]) / (CONSTANTS[0] * CONSTANTS[1])) - 1.00000));
	ALGEBRAIC[23] = STATES[7] * STATES[8] * ALGEBRAIC[22] * ALGEBRAIC[17];
	ALGEBRAIC[26] = ALGEBRAIC[23] + ALGEBRAIC[25] + ALGEBRAIC[24];
	ALGEBRAIC[38] = CONSTANTS[26] * (STATES[5] / (CONSTANTS[25] + STATES[5]));
	ALGEBRAIC[41] = ((CONSTANTS[0] * CONSTANTS[1]) / (2.00000*CONSTANTS[2]))*log(CONSTANTS[18] / STATES[5]);
	ALGEBRAIC[42] = CONSTANTS[28] * (STATES[0] - ALGEBRAIC[41]);
	ALGEBRAIC[51] = ALGEBRAIC[48] + ALGEBRAIC[50];
	ALGEBRAIC[53] = (ALGEBRAIC[3] - (ALGEBRAIC[14] + ALGEBRAIC[26] + ALGEBRAIC[29] + ALGEBRAIC[34] + ALGEBRAIC[37] + ALGEBRAIC[52] + ALGEBRAIC[38] + ALGEBRAIC[40] + ALGEBRAIC[42] + ALGEBRAIC[44] + ALGEBRAIC[51])) / CONSTANTS[3];
	ALGEBRAIC[54] = CONSTANTS[59] * (STATES[10] - STATES[5]);
	ALGEBRAIC[55] = CONSTANTS[46] * (STATES[5] / (STATES[5] + CONSTANTS[45]));
	ALGEBRAIC[56] = CONSTANTS[58] * STATES[11];
	ALGEBRAIC[57] = (STATES[11] - STATES[10]) / CONSTANTS[43];
}

void run_LuoRudyRK2()
{

	double* ALGEBRAIC = newDoubleV(58);
	double* CONSTANTS = newDoubleV(60);
	double* RATES = newDoubleV(12);
	double* STATES = newDoubleV(12);

	initConsts_LuoRudy(CONSTANTS, RATES, STATES);

	double dt = 0.0075;
	double tmax = 500;
	double tinit = 0;
	int steps = (int)(tmax - tinit) / dt;

	double** Values = newDouble(steps, 12);

	char** Names[12];
	Names[0] = "0";
	Names[1] = "1";
	Names[2] = "2";
	Names[3] = "3";
	Names[4] = "4";
	Names[5] = "5";
	Names[6] = "6";
	Names[7] = "7";
	Names[8] = "8";
	Names[9] = "9";
	Names[10] = "10";
	Names[11] = "11";


	double t = 0;
	for (int i = 0; i < steps; i++)
	{

		double* k1 = newDoubleV(12);
		computeRates_LuoRudy(t, CONSTANTS, k1, STATES, ALGEBRAIC);

		double* k2 = newDoubleV(12);
		double* tempk1 = multConstantV(k1, dt / 2, 12);
		double* statesk1 = addVector(STATES, tempk1, 12);
		computeRates_LuoRudy(t + dt / 2, CONSTANTS, k2, statesk1, ALGEBRAIC);


		for (int j = 0; j < 12; j++)
		{
			STATES[j] += k2[j] * dt;
			Values[i][j] = STATES[j];
		}

		free(k1);
		free(k2);
		free(tempk1);
		free(statesk1);

		t += dt;
	}

	plot2DComparisonNamesNew(Values, 12, steps, "Results", "Graph.html", Names, "Amplitude");

	freeMatrix(Values, steps);

	free(CONSTANTS);
	free(RATES);
	free(STATES);
	free(ALGEBRAIC);

}

void run_LuoRudyRK4()
{

	double* ALGEBRAIC = newDoubleV(58);
	double* CONSTANTS = newDoubleV(60);
	double* RATES = newDoubleV(12);
	double* STATES = newDoubleV(12);

	initConsts_LuoRudy(CONSTANTS, RATES, STATES);

	double dt = 0.01;
	double tmax = 500;
	double tinit = 0;
	int steps = (int)(tmax - tinit) / dt;

	double** Values = newDouble(steps, 12);

	char** Names[12];
	Names[0] = "0";
	Names[1] = "1";
	Names[2] = "2";
	Names[3] = "3";
	Names[4] = "4";
	Names[5] = "5";
	Names[6] = "6";
	Names[7] = "7";
	Names[8] = "8";
	Names[9] = "9";
	Names[10] = "10";
	Names[11] = "11";


	double t = 0;
	for (int i = 0; i < steps; i++)
	{

		double* k1 = newDoubleV(12);
		computeRates_LuoRudy(t, CONSTANTS, k1, STATES, ALGEBRAIC);

		double* k2 = newDoubleV(12);
		double* tempk1 = multConstantV(k1, dt / 2, 12);
		double* statesk1 = addVector(STATES, tempk1, 12);
		computeRates_LuoRudy(t + dt / 2, CONSTANTS, k2, statesk1, ALGEBRAIC);

		double* k3 = newDoubleV(12);
		double* tempk2 = multConstantV(k2, dt / 2, 12);
		double* statesk2 = addVector(STATES, tempk2, 12);
		computeRates_LuoRudy(t + dt / 2, CONSTANTS, k3, statesk2, ALGEBRAIC);

		double* k4 = newDoubleV(12);
		double* tempk3 = multConstantV(k3, dt, 12);
		double* statesk3 = addVector(STATES, tempk3, 12);
		computeRates_LuoRudy(t + dt, CONSTANTS, k4, statesk3, ALGEBRAIC);

		for (int j = 0; j < 12; j++)
		{
			STATES[j] += dt / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
			Values[i][j] = STATES[j];
		}

		free(k1);
		
		free(k2);
		free(tempk1);
		free(statesk1);

	
		free(k3);
		free(tempk2);
		free(statesk2);

		
		free(k4);
		free(tempk3);
		free(statesk3);

		t += dt;
	}

	plot2DComparisonNamesNew(Values, 12, steps, "Results", "Graph.html", Names, "Amplitude");

	freeMatrix(Values, steps);

	free(CONSTANTS);
	free(RATES);
	free(STATES);
	free(ALGEBRAIC);

}

/*
There are a total of 10 entries in the algebraic variable array.
There are a total of 4 entries in each of the rate and state variable arrays.
There are a total of 8 entries in the constant variable array.
*/
/*
* VOI is time in component environment (millisecond).
* STATES[0] is V in component membrane (millivolt).
* CONSTANTS[0] is E_R in component membrane (millivolt).
* CONSTANTS[1] is Cm in component membrane (microF_per_cm2).
* ALGEBRAIC[4] is i_Na in component sodium_channel (microA_per_cm2).
* ALGEBRAIC[8] is i_K in component potassium_channel (microA_per_cm2).
* ALGEBRAIC[9] is i_L in component leakage_current (microA_per_cm2).
* ALGEBRAIC[0] is i_Stim in component membrane (microA_per_cm2).
* CONSTANTS[2] is g_Na in component sodium_channel (milliS_per_cm2).
* CONSTANTS[5] is E_Na in component sodium_channel (millivolt).
* STATES[1] is m in component sodium_channel_m_gate (dimensionless).
* STATES[2] is h in component sodium_channel_h_gate (dimensionless).
* ALGEBRAIC[1] is alpha_m in component sodium_channel_m_gate (per_millisecond).
* ALGEBRAIC[5] is beta_m in component sodium_channel_m_gate (per_millisecond).
* ALGEBRAIC[2] is alpha_h in component sodium_channel_h_gate (per_millisecond).
* ALGEBRAIC[6] is beta_h in component sodium_channel_h_gate (per_millisecond).
* CONSTANTS[3] is g_K in component potassium_channel (milliS_per_cm2).
* CONSTANTS[6] is E_K in component potassium_channel (millivolt).
* STATES[3] is n in component potassium_channel_n_gate (dimensionless).
* ALGEBRAIC[3] is alpha_n in component potassium_channel_n_gate (per_millisecond).
* ALGEBRAIC[7] is beta_n in component potassium_channel_n_gate (per_millisecond).
* CONSTANTS[4] is g_L in component leakage_current (milliS_per_cm2).
* CONSTANTS[7] is E_L in component leakage_current (millivolt).
* RATES[0] is d/dt V in component membrane (millivolt).
* RATES[1] is d/dt m in component sodium_channel_m_gate (dimensionless).
* RATES[2] is d/dt h in component sodium_channel_h_gate (dimensionless).
* RATES[3] is d/dt n in component potassium_channel_n_gate (dimensionless).
*/
void
initConsts_Hodgkin(double* CONSTANTS, double* RATES, double *STATES)
{
	STATES[0] = -75;
	CONSTANTS[0] = -75;
	CONSTANTS[1] = 1;
	CONSTANTS[2] = 120;
	STATES[1] = 0.05;
	STATES[2] = 0.6;
	CONSTANTS[3] = 36;
	STATES[3] = 0.325;
	CONSTANTS[4] = 0.3;
	CONSTANTS[5] = CONSTANTS[0] + 115.000;
	CONSTANTS[6] = CONSTANTS[0] - 12.0000;
	CONSTANTS[7] = CONSTANTS[0] + 10.6130;
}
void
computeRates_Hodgkin(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
	ALGEBRAIC[1] = (-0.100000*(STATES[0] + 50.0000)) / (exp(-(STATES[0] + 50.0000) / 10.0000) - 1.00000);
	ALGEBRAIC[5] = 4.00000*exp(-(STATES[0] + 75.0000) / 18.0000);
	RATES[1] = ALGEBRAIC[1] * (1.00000 - STATES[1]) - ALGEBRAIC[5] * STATES[1];
	ALGEBRAIC[2] = 0.0700000*exp(-(STATES[0] + 75.0000) / 20.0000);
	ALGEBRAIC[6] = 1.00000 / (exp(-(STATES[0] + 45.0000) / 10.0000) + 1.00000);
	RATES[2] = ALGEBRAIC[2] * (1.00000 - STATES[2]) - ALGEBRAIC[6] * STATES[2];
	ALGEBRAIC[3] = (-0.0100000*(STATES[0] + 65.0000)) / (exp(-(STATES[0] + 65.0000) / 10.0000) - 1.00000);
	ALGEBRAIC[7] = 0.125000*exp((STATES[0] + 75.0000) / 80.0000);
	RATES[3] = ALGEBRAIC[3] * (1.00000 - STATES[3]) - ALGEBRAIC[7] * STATES[3];
	ALGEBRAIC[4] = CONSTANTS[2] * pow(STATES[1], 3.00000)*STATES[2] * (STATES[0] - CONSTANTS[5]);
	ALGEBRAIC[8] = CONSTANTS[3] * pow(STATES[3], 4.00000)*(STATES[0] - CONSTANTS[6]);
	ALGEBRAIC[9] = CONSTANTS[4] * (STATES[0] - CONSTANTS[7]);
	ALGEBRAIC[0] = (VOI >= 10.0000&&VOI <= 10.5000 ? 20.0000 : 0.00000);
	RATES[0] = -(-ALGEBRAIC[0] + ALGEBRAIC[4] + ALGEBRAIC[8] + ALGEBRAIC[9]) / CONSTANTS[1];
}
void
computeVariables_Hodgkin(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
	ALGEBRAIC[1] = (-0.100000*(STATES[0] + 50.0000)) / (exp(-(STATES[0] + 50.0000) / 10.0000) - 1.00000);
	ALGEBRAIC[5] = 4.00000*exp(-(STATES[0] + 75.0000) / 18.0000);
	ALGEBRAIC[2] = 0.0700000*exp(-(STATES[0] + 75.0000) / 20.0000);
	ALGEBRAIC[6] = 1.00000 / (exp(-(STATES[0] + 45.0000) / 10.0000) + 1.00000);
	ALGEBRAIC[3] = (-0.0100000*(STATES[0] + 65.0000)) / (exp(-(STATES[0] + 65.0000) / 10.0000) - 1.00000);
	ALGEBRAIC[7] = 0.125000*exp((STATES[0] + 75.0000) / 80.0000);
	ALGEBRAIC[4] = CONSTANTS[2] * pow(STATES[1], 3.00000)*STATES[2] * (STATES[0] - CONSTANTS[5]);
	ALGEBRAIC[8] = CONSTANTS[3] * pow(STATES[3], 4.00000)*(STATES[0] - CONSTANTS[6]);
	ALGEBRAIC[9] = CONSTANTS[4] * (STATES[0] - CONSTANTS[7]);
	ALGEBRAIC[0] = (VOI >= 10.0000&&VOI <= 10.5000 ? 20.0000 : 0.00000);
}

void run_Hodgkin()
{
	double* ALGEBRAIC = newDoubleV(10);
	double* CONSTANTS = newDoubleV(8);
	double* RATES = newDoubleV(4);
	double* STATES = newDoubleV(4);

	initConsts_Hodgkin(CONSTANTS, RATES, STATES);

	double dt = 0.00125;
	double tmax = 50;
	double tinit = 0;
	int steps = (int)(tmax - tinit) / dt;

	double** Values = newDouble(steps, 4);

	char** Names[4];
	Names[0] = "0";
	Names[1] = "1";
	Names[2] = "2";
	Names[3] = "3";



	double t = 0;
	for (int i = 0; i < steps; i++)
	{


		double Istim = 0;

		computeRates_Hodgkin(t, CONSTANTS, RATES, STATES, ALGEBRAIC);

		eulerODE(RATES, STATES, 4, dt);
		computeVariables_Hodgkin(t, CONSTANTS, RATES, STATES, ALGEBRAIC);
		for (int j = 0; j < 4; j++)
		{
			Values[i][j] = STATES[j];
		}

		t += dt;

	}

	plot2DComparisonNamesNew(Values, 4, steps, "Results", "Graph.html", Names, "Amplitude");

	freeMatrix(Values, steps);

	free(CONSTANTS);
	free(RATES);
	free(STATES);
	free(ALGEBRAIC);

}

/*
There are a total of 8 entries in the algebraic variable array.
There are a total of 3 entries in each of the rate and state variable arrays.
There are a total of 22 entries in the constant variable array.
*/
/*
* VOI is time in component environment (ms).
* STATES[0] is u in component membrane (dimensionless).
* CONSTANTS[0] is Cm in component membrane (uF_per_cm2).
* ALGEBRAIC[0] is Vm in component membrane (mV).
* CONSTANTS[1] is V_0 in component membrane (mV).
* CONSTANTS[2] is V_fi in component membrane (mV).
* ALGEBRAIC[2] is J_fi in component fast_inward_current (per_ms).
* ALGEBRAIC[4] is J_so in component slow_outward_current (per_ms).
* ALGEBRAIC[6] is J_si in component slow_inward_current (per_ms).
* ALGEBRAIC[7] is Istim in component stimulus_protocol (per_ms).
* ALGEBRAIC[1] is p in component p (dimensionless).
* CONSTANTS[3] is u_c in component p (dimensionless).
* ALGEBRAIC[3] is q in component q (dimensionless).
* CONSTANTS[4] is u_v in component q (dimensionless).
* CONSTANTS[21] is tau_d in component fast_inward_current (ms).
* CONSTANTS[5] is g_fi_max in component fast_inward_current (mS_per_cm2).
* STATES[1] is v in component fast_inward_current_v_gate (dimensionless).
* ALGEBRAIC[5] is tau_v_minus in component fast_inward_current_v_gate (ms).
* CONSTANTS[6] is tau_v1_minus in component fast_inward_current_v_gate (ms).
* CONSTANTS[7] is tau_v2_minus in component fast_inward_current_v_gate (ms).
* CONSTANTS[8] is tau_v_plus in component fast_inward_current_v_gate (ms).
* CONSTANTS[9] is tau_0 in component slow_outward_current (ms).
* CONSTANTS[10] is tau_r in component slow_outward_current (ms).
* CONSTANTS[11] is tau_si in component slow_inward_current (ms).
* CONSTANTS[12] is u_csi in component slow_inward_current (dimensionless).
* CONSTANTS[13] is k in component slow_inward_current (dimensionless).
* STATES[2] is w in component slow_inward_current_w_gate (dimensionless).
* CONSTANTS[14] is tau_w_minus in component slow_inward_current_w_gate (ms).
* CONSTANTS[15] is tau_w_plus in component slow_inward_current_w_gate (ms).
* CONSTANTS[16] is IstimStart in component stimulus_protocol (ms).
* CONSTANTS[17] is IstimEnd in component stimulus_protocol (ms).
* CONSTANTS[18] is IstimAmplitude in component stimulus_protocol (per_ms).
* CONSTANTS[19] is IstimPeriod in component stimulus_protocol (ms).
* CONSTANTS[20] is IstimPulseDuration in component stimulus_protocol (ms).
* RATES[0] is d/dt u in component membrane (dimensionless).
* RATES[1] is d/dt v in component fast_inward_current_v_gate (dimensionless).
* RATES[2] is d/dt w in component slow_inward_current_w_gate (dimensionless).
*/
void
initConsts_FemtomKarmaMod(double* CONSTANTS, double* RATES, double *STATES)
{
	STATES[0] = 0;
	CONSTANTS[0] = 1;
	CONSTANTS[1] = -85;
	CONSTANTS[2] = 15;
	CONSTANTS[3] = 0.13;
	CONSTANTS[4] = 0.04;
	CONSTANTS[5] = 4;
	STATES[1] = 1;
	CONSTANTS[6] = 1250;
	CONSTANTS[7] = 19.6;
	CONSTANTS[8] = 3.33;
	CONSTANTS[9] = 12.5;
	CONSTANTS[10] = 33.33;
	CONSTANTS[11] = 29;
	CONSTANTS[12] = 0.85;
	CONSTANTS[13] = 10;
	STATES[2] = 1;
	CONSTANTS[14] = 41;
	CONSTANTS[15] = 870;
	CONSTANTS[16] = 10;
	CONSTANTS[17] = 50000;
	CONSTANTS[18] = -0.2;
	CONSTANTS[19] = 1000;
	CONSTANTS[20] = 1;
	CONSTANTS[21] = CONSTANTS[0] / CONSTANTS[5];
}
void
computeRates_FemtomKarmaMod(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, double Istim)
{
	ALGEBRAIC[1] = (STATES[0]<CONSTANTS[3] ? 0.00000 : 1.00000);
	RATES[2] = ((1.00000 - ALGEBRAIC[1])*(1.00000 - STATES[2])) / CONSTANTS[14] - (ALGEBRAIC[1] * STATES[2]) / CONSTANTS[15];
	ALGEBRAIC[3] = (STATES[0]<CONSTANTS[4] ? 0.00000 : 1.00000);
	ALGEBRAIC[5] = ALGEBRAIC[3] * CONSTANTS[6] + (1.00000 - ALGEBRAIC[3])*CONSTANTS[7];
	RATES[1] = ((1.00000 - ALGEBRAIC[1])*(1.00000 - STATES[1])) / ALGEBRAIC[5] - (ALGEBRAIC[1] * STATES[1]) / CONSTANTS[8];
	ALGEBRAIC[2] = (-STATES[1] * ALGEBRAIC[1] * (1.00000 - STATES[0])*(STATES[0] - CONSTANTS[3])) / CONSTANTS[21];
	ALGEBRAIC[4] = (STATES[0] * (1.00000 - ALGEBRAIC[1])) / CONSTANTS[9] + ALGEBRAIC[1] / CONSTANTS[10];
	ALGEBRAIC[6] = (-STATES[2] * (1.00000 + tanh(CONSTANTS[13] * (STATES[0] - CONSTANTS[12])))) / (2.00000*CONSTANTS[11]);
	//ALGEBRAIC[7] = (VOI >= CONSTANTS[16] && VOI <= CONSTANTS[17] && (VOI - CONSTANTS[16]) - floor((VOI - CONSTANTS[16]) / CONSTANTS[19])*CONSTANTS[19] <= CONSTANTS[20] ? CONSTANTS[18] : 0.00000);
	ALGEBRAIC[7] = Istim;
	RATES[0] = -(ALGEBRAIC[2] + ALGEBRAIC[4] + ALGEBRAIC[6] + ALGEBRAIC[7]);
}
void
computeVariables_FemtomKarmaMod(double VOI, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC)
{
	ALGEBRAIC[1] = (STATES[0]<CONSTANTS[3] ? 0.00000 : 1.00000);
	ALGEBRAIC[3] = (STATES[0]<CONSTANTS[4] ? 0.00000 : 1.00000);
	ALGEBRAIC[5] = ALGEBRAIC[3] * CONSTANTS[6] + (1.00000 - ALGEBRAIC[3])*CONSTANTS[7];
	ALGEBRAIC[2] = (-STATES[1] * ALGEBRAIC[1] * (1.00000 - STATES[0])*(STATES[0] - CONSTANTS[3])) / CONSTANTS[21];
	ALGEBRAIC[4] = (STATES[0] * (1.00000 - ALGEBRAIC[1])) / CONSTANTS[9] + ALGEBRAIC[1] / CONSTANTS[10];
	ALGEBRAIC[6] = (-STATES[2] * (1.00000 + tanh(CONSTANTS[13] * (STATES[0] - CONSTANTS[12])))) / (2.00000*CONSTANTS[11]);
	ALGEBRAIC[7] = (VOI >= CONSTANTS[16] && VOI <= CONSTANTS[17] && (VOI - CONSTANTS[16]) - floor((VOI - CONSTANTS[16]) / CONSTANTS[19])*CONSTANTS[19] <= CONSTANTS[20] ? CONSTANTS[18] : 0.00000);
	ALGEBRAIC[0] = CONSTANTS[1] + STATES[0] * (CONSTANTS[2] - CONSTANTS[1]);
}

void run_FemtomKarmaMod()
{
	double* ALGEBRAIC = newDoubleV(8);
	double* CONSTANTS = newDoubleV(22);
	double* RATES = newDoubleV(3);
	double* STATES = newDoubleV(3);

	initConsts_FemtomKarmaMod(CONSTANTS, RATES, STATES);

	double dt = 0.25;
	double tmax = 2000;
	double tinit = 0;
	int steps = (int)(tmax - tinit) / dt;

	double** Values = newDouble(steps, 3);

	char** Names[3];
	Names[0] = "U";
	Names[1] = "V";
	Names[2] = "W";



	double t = 0;
	for (int i = 0; i < steps; i++)
	{
		double Istim = 0;
		if (t>100 && t<101)
		Istim = -0.2;

		computeRates_FemtomKarmaMod(t, CONSTANTS, RATES, STATES, ALGEBRAIC, Istim);

		eulerODE(RATES, STATES, 3, dt);

		for (int j = 0; j < 3; j++)
		{
			Values[i][j] = STATES[j];
		}

		t += dt;

	}

	plot2DComparisonNamesNew(Values, 3, steps, "Results", "Graph.html", Names, "Amplitude");

	freeMatrix(Values, steps);

	free(CONSTANTS);
	free(RATES);
	free(STATES);
	free(ALGEBRAIC);

}
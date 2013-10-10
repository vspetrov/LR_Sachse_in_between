#include "Sachse_fibroblast.h"


double Vf_function(double V, double O_shkr, double Iext)
{
	double I_kir, I_shkr, I_b;
	double O_kir = 1./(a_kir+exp(b_kir*(V-E_k)*FbyRT));
	I_kir = G_kir*O_kir*sqrt(K_0/1000.)*(V-E_k);

	I_shkr = P_shkr*O_shkr*V*F*FbyRT*(K_i-K_0*exp(-V*FbyRT))/(1.-exp(-V*FbyRT));

	I_b = G_b*(V-E_b);

    return -1./Cf*(I_kir+I_b+I_shkr)+Iext;
}

double C0_function(double C0, double C1, double V)
{
	double kv, k_v;
	kv = kv0*exp(V*zv*FbyRT);
	k_v = k_v0*exp(V*z_v*FbyRT);

	return -4.*kv*C0+k_v*C1;
}

double C1_function(double C0, double C1, double C2, double V)
{
	double kv, k_v;
	kv = kv0*exp(V*zv*FbyRT);
	k_v = k_v0*exp(V*z_v*FbyRT);

	return 4.*kv*C0-(3.*kv+k_v)*C1+2.*k_v*C2;
}

double C2_function(double C1, double C2, double C3, double V)
{
	double kv, k_v;
	kv = kv0*exp(V*zv*FbyRT);
	k_v = k_v0*exp(V*z_v*FbyRT);

	return 3.*kv*C1-(2.*kv+2.*k_v)*C2+3.*k_v*C3;
}

double C3_function(double C2, double C3, double C4, double V)
{
	double kv, k_v;
	kv = kv0*exp(V*zv*FbyRT);
	k_v = k_v0*exp(V*z_v*FbyRT);

	return 2.*kv*C2-(kv+3.*k_v)*C3+4.*k_v*C4;
}

double C4_function(double C3, double C4, double O_shkr, double V)
{
	double kv, k_v;
	kv = kv0*exp(V*zv*FbyRT);
	k_v = k_v0*exp(V*z_v*FbyRT);

	return kv*C3-(ko+4.*k_v)*C4+k_o*O_shkr;
}

double O_function(double C4, double O_shkr)
{
	return ko*C4-k_o*O_shkr;
}

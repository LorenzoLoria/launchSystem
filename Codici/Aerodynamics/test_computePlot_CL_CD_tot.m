% test_computePlot_CL_CD_tot.m / test_computePlot_CL_CD_tot_2.m
clear;
clc;

%% dati

vector.M = [0 : 0.01 : 5];
vector.M_cases = [0.1 0.8 2 7];
vector.alpha_deg = [0 : 0.01 : 15];
vector.alpha_deg_cases = [3];



bodyGeom.l = 3;
bodyGeom.d = 0.3;
bodyGeom.Lnose = 0.6;
bodyGeom.Aref = pi*(bodyGeom.d/2)^2; 
bodyGeom.Anose = bodyGeom.Aref;
bodyGeom.Abase = bodyGeom.Aref;
bodyGeom.Aexit = [];  
bodyGeom.phi = deg2rad(0);
bodyGeom.Ab = bodyGeom.Aref;
bodyGeom.Ap = bodyGeom.l * bodyGeom.d;



cr   = 0.35;                      % root chord [m]
ct   = 0.0;                       % tip chord [m] (triangolo puro)
s    = 0.20;                      % semispan [m]

finsGeom.Nfins = 1;
finsGeom.be = s;                  % span equivalente ~ semispan reale
finsGeom.Se = 0.5 * cr * s;       % area della fin [m^2]
finsGeom.cmac = (2/3) * cr;       % mean aerodynamic chord [m]
finsGeom.delta_le = 45;
finsGeom.lambda_le = 0;
finsGeom.b = 2 * cr;
finsGeom.tmac = 0.08 * finsGeom.cmac;    



bodyInfo.isPowered = false;
bodyInfo.a_sub = 0;
bodyInfo.b_sub = 1;
bodyInfo.Cdn = 1.2;



finsInfo.rho = 1.225;
finsInfo.vsound = 340;



gasProp.gamma = 1.4;



%% test

[CL_tot_alpha, CD_tot_mach] = computePlot_CL_CD_tot(vector, bodyGeom, finsGeom, bodyInfo, finsInfo);

[CL_tot_alpha2, CD_tot_mach2] = computePlot_CL_CD_tot_2(vector, bodyGeom, finsGeom, bodyInfo, finsInfo, gasProp);
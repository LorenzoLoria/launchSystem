% test_computePlot_CL_CD_tot.m
clear;
clc;

% dati

vector.M = [0 : 0.01 : 5];
vector.M_cases = [0.8 1.3 2 5];
vector.alpha_deg = [0 : 0.01 : 15];
vector.alpha_cases = [0 10 30 50];

bodyGeom.L = 3;
bodyGeom.l = bodyGeom.L;
bodyGeom.d = 0.3;
bodyGeom.Lnose = 0.6;
bodyGeom.Aref = pi*(bodyGeom.d/2)^2; 
bodyGeom.Anose = bodyGeom.Aref;
bodyGeom.Abase = bodyGeom.Aref;
bodyGeom.Aexit = [];  
bodyGeom.phi = deg2rad(0);
bodyGeom.Ab = bodyGeom.Aref;
bodyGeom.Ap = bodyGeom.L * bodyGeom.d;

finsGeom.Nfins = 1;
finsGeom.Se = 0.5 * 0.35 * 0.20; 
finsGeom.cmac = (2/3) * 0.35; 
finsGeom.delta_le = 45;
finsGeom.lambda_le = 0;
finsGeom.b = 2 * 0.35;
finsGeom.tmac = 0.08 * finsGeom.cmac;    

bodyInfo.isPowered = false;
bodyInfo.a_sub = 0;
bodyInfo.b_sub = 1;
bodyInfo.Cdn = 1.2;

finsInfo.rho = 1.225;
finsInfo.vsound = 340;



[CL_tot_alpha, CD_tot_mach] = computePlot_CL_CD_tot(vector, bodyGeom, finsGeom, bodyInfo, finsInfo);
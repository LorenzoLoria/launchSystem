% test_computePlot_CL_CD_tot.m / test_computePlot_CL_CD_tot_2.m
%% ======================= Prova Runnata ================================== 

clear all;
clc;
close all

% ==========================  DATI ========================================

addpath(genpath('..\..\'))

[mission, settings] = dataStructGlobal;

mission.structure.alphaQmax = 0;

% launcher = [nStages, nMotore1, nMotore2, nMotore3, %massa1, %massa2,
% %massa3];
launcher = [2,2,3,4,0.459952176990556, 0.753370531158904, 0.634795741885559];

for i = 1:launcher(1)
    configuration.stage{i}.engine = mission.engines{launcher(1+i)};
end

[mer,staging,configuration] = initialMassEstimation(mission,configuration,settings,launcher);
thrustDataGA = load('thrustdataVecTraj.mat','xGATraj');

thrustDataVecFMC(:,:,1) = [0.902082365568723	1.480898931628005
                            0.999984156345040	23.253294859564580
                            0.900002678098914	52.979241033086943
                            0.900000000000007	59.571815331701984
                            0.903941814015555	55.058714159781090];


thrustDataVecFMC(:,:,2 ) = [0.400917809388214	65.122710138507202
                            0.964494359624014	79.658359202140389
                            0.975968800776448	91.801043507018605
                            0.992714640706230	89.085172454227390
                            0.993244065056187	99.345740209598944];



[timeCollocation, stateCollocation] = totalTrajectoryGlobalGA(launcher,configuration,mission,settings,thrustDataVecFMC);



%% dati

vector.M = [0 : 0.01 : 25];
%vector.M = [2];
vector.M_cases = [0.1 0.5 1.1 2 5];
%vector.M_cases = [2];
vector.alpha_deg = [0 : 0.01 : 15];
%vector.alpha_deg = [3];
vector.alpha_deg_cases = [0];



bodyGeom.l = configuration.geometry.totalLength;
bodyGeom.d = mission.capsule.radius * 2;
bodyGeom.Lnose = mission.capsule.height;
bodyGeom.Aref = pi*(bodyGeom.d/2)^2; 
bodyGeom.Anose = pi*(mission.capsule.radius/1.8)^2; 
bodyGeom.Abase = bodyGeom.Aref;
bodyGeom.Aexit = mission.engines{2}.effAreaZero;  
bodyGeom.phi = deg2rad(0);
bodyGeom.Ab = bodyGeom.Aref;
bodyGeom.Ap = bodyGeom.l * bodyGeom.d;



cr   = 1.81;                       % root chord [m]
ct   = 0.45;                       % tip chord [m] (triangolo puro)     

finsGeom.Nfins = 4;
finsGeom.be = 1.81;                  % span equivalente ~ semispan reale
finsGeom.Se = 2.05;                  % area della fin [m^2]
finsGeom.cmac = 1.27;                % mean aerodynamic chord [m]
finsGeom.delta_le = 30; 
finsGeom.lambda_le = 36.9;
finsGeom.b = finsGeom.be;
finsGeom.tmac = 0.10;    




% ============================ Fins section ===============================
% Hexagonal
%
% --> Sezione root
% delta_le = 40 deg (angolo totale)
% tratto1 = 0.1374
% tratto2 = 1.5352
% tratto3 = 0.1374
% 
% --> Sezione tip
% delta_le = 40 deg (angolo totale)
% tratto1 = 0.1374
% tratto2 = 0.752
% tratto3 = 0.1374




bodyInfo.isPowered = true;
bodyInfo.a_sub = 0;
bodyInfo.b_sub = 1;
bodyInfo.Cdn = 1.2;


% vel = stateCollocation(4:6,:,1:end-1)-stateCollocation(4:6,1,1);
% vel = mission.target.Rfinal* vel(1:3,:);
% absVel = sqrt ( vel(1,:).^2+ vel(2,:).^2 + vel(3,:).^2 );
% 
% pos = stateCollocation(1:3,:,1:end-1);
% pos = pos(1:3,:);
% absH = sqrt ( pos(1,:).^2+ pos(2,:).^2 + pos(3,:).^2 )-mission.environment.rEarth;
% 
% rhoVec = mission.environment.rhoFun(absH);
% 
% % q = 0.5*rhoVec.*(absVel).^2;
% 
% stageNumber = 2;
% alpha = 0;
% 
% finsInfo.vsound = mission.aerodynamics.soundspeedFun(absH);
% 
% finsInfo.rho = rhoVec;


finsInfo.vsound = 340;
finsInfo.rho = 1.225;
gasProp.gamma = 1.4;



%% test

%[CL_tot_alpha, CD_tot_mach] = computePlot_CL_CD_tot(vector, bodyGeom, finsGeom, bodyInfo, finsInfo);

[CL_tot_alpha2, CD_tot_mach2] = computePlot_CL_CD_tot_2(vector, bodyGeom, finsGeom, bodyInfo, finsInfo, gasProp, settings);
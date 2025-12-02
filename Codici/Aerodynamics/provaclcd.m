
mission.aerodynamics.bodyGeom.l = 3;
mission.aerodynamics.bodyGeom.d = 0.3;
mission.aerodynamics.bodyGeom.Lnose = 0.6;
mission.aerodynamics.bodyGeom.Aref = pi*(mission.aerodynamics.bodyGeom.d/2)^2; 
mission.aerodynamics.bodyGeom.Anose = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Abase = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Aexit = 0;  
mission.aerodynamics.bodyGeom.phi = deg2rad(0);
mission.aerodynamics.bodyGeom.Ab = mission.aerodynamics.bodyGeom.Aref;
mission.aerodynamics.bodyGeom.Ap = mission.aerodynamics.bodyGeom.l * mission.aerodynamics.bodyGeom.d;



cr   = 0.35;                      % root chord [m]
ct   = 0.0;                       % tip chord [m] (triangolo puro)
s    = 0.20;                      % semispan [m]

mission.aerodynamics.finsGeom.Nfins = 1;
mission.aerodynamics.finsGeom.be = s;                  % span equivalente ~ semispan reale
mission.aerodynamics.finsGeom.Se = 0.5 * cr * s;       % area della fin [m^2]
mission.aerodynamics.finsGeom.cmac = (2/3) * cr;       % mean aerodynamic chord [m]
mission.aerodynamics.finsGeom.delta_le = 45;
mission.aerodynamics.finsGeom.lambda_le = 0;
mission.aerodynamics.finsGeom.b = 2 * cr;
mission.aerodynamics.finsGeom.tmac = 0.08 * mission.aerodynamics.finsGeom.cmac;    



mission.aerodynamics.bodyInfo.isPowered = false;
mission.aerodynamics.bodyInfo.a_sub = 0;
mission.aerodynamics.bodyInfo.b_sub = 1;
mission.aerodynamics.bodyInfo.Cdn = 1.2;



mission.aerodynamics.finsInfo.rho = 1.225;
mission.aerodynamics.finsInfo.vsound = 340;


[CL,CD] = CLCDcomputation(2,deg2rad(3),283220,1,mission);


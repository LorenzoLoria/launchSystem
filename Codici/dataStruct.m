% data script
function [mission] = dataStruct

mission = struct();

mission.launcher.engines{1}.name = 'Merlin1D' ;
mission.launcher.engines{2}.name = 'Raptor' ;
mission.launcher.engines{3}.name = 'RS-25' ;

mission.launcher.engines{1}.isp = 283;
mission.launcher.engines{1}.thrust = 854*1e3;
mission.launcher.engines{1}.weight = 470;
mission.launcher.engines{1}.OF = 2.36;
mission.launcher.engines{1}.oxDens = 1143;
mission.launcher.engines{1}.fuelDens = 835;
 
mission.launcher.engines{2}.isp = 327;
mission.launcher.engines{2}.thrust = 2750*1e3;
mission.launcher.engines{2}.weight = 1525;
mission.launcher.engines{2}.OF = 3.6;
mission.launcher.engines{2}.oxDens = 1143;
mission.launcher.engines{2}.fuelDens = 423;

mission.launcher.engines{3}.isp = 366;
mission.launcher.engines{3}.thrust = 1860*1e3;
mission.launcher.engines{3}.weight = 3526;
mission.launcher.engines{3}.OF = 6;
mission.launcher.engines{3}.oxDens = 1143;
mission.launcher.engines{3}.fuelDens = 71;

mission.capsule.weigth = 8600;
mission.capsule.Area = 3.7^2*pi;
mission.capsule.supersonicCD = 1.23;
mission.capsule.subsonicCD = 0.45;

mission.environment.altRange = (-1000:100:1000000);

warning('off', 'all'); 
for i=1:length(mission.environment.altRange)
    [T,rho] = atmosnrlmsise00(mission.environment.altRange(i), 45, 120, 2020,50,1300);
    rhoVal(i) = rho(6);
    Tval(i) = T(1,2) ; 
end
warning('on', 'all');   
mission.environment.rho = rhoVal;
mission.environment.T = Tval;
mission.environment.rEarth = 6371e3;
mission.environment.g0 = 9.81;
mission.environment.GM = 398600.4e9;
end
function [mer,staging] = initialMassEstimation(mission,opt,settings,launcher)

nStages = launcher(1);
lOverD = 5;
    avionics =350; % 10 * (M0)^0.361; the previous was overyly conservative as per il libro che ho ciulato
    wiring = 10;

    if nStages == 1
        mer.stage{1}.avionics = avionics;
        mer.stage{1}.wiring = wiring;
        mer.stage{1}.battery = 20;
        mer.stage{1}.payloadAttach = 0.0775*mission.capsule.weigth + 50;
    elseif nStages == 2
        mer.stage{1}.avionics = 0.2*avionics;
        mer.stage{2}.avionics = 0.8*avionics;
        
        mer.stage{1}.wiring = 0.2*wiring;
        mer.stage{2}.wiring = 0.8*wiring;

        mer.stage{1}.battery = 20;
        mer.stage{2}.battery = 8;

        mer.stage{2}.payloadAttach = 0.0775*mission.capsule.weigth + 50;        
    else
        mer.stage{1}.avionics = 0.1*avionics;
        mer.stage{2}.avionics = 0.1*avionics;
        mer.stage{3}.avionics = 0.8*avionics;
    
        mer.stage{1}.wiring = 0.1*wiring;
        mer.stage{2}.wiring = 0.1*wiring;
        mer.stage{3}.wiring = 0.8*wiring;

        mer.stage{1}.battery = 20;
        mer.stage{2}.battery = 8;
        mer.stage{3}.battery = 8;

        mer.stage{3}.payloadAttach = 0.0775*mission.capsule.weigth + 50;        
    end
    
    mStack = mission.capsule.weigth;

    for i = nStages:-1:1
        couple = opt.stage{i}.engine.couple;

        nEngine = 0;
        twRatio = 0;
        mStage = 0;

        if i == 1
            thrustValue = opt.stage{i}.engine.thrustZero;
        else
            thrustValue = opt.stage{i}.engine.thrustVacum;
        end



        while twRatio < 1
            mStack = mStack-mStage;
            nEngine = nEngine+1;
            mProp = nEngine * thrustValue/9.81 * launcher(3+i);

            if strcmp(couple, 'RP1-LOX')

                mFuel = 1/(1+opt.stage{i}.engine.OF)*mProp;
                mOx = opt.stage{i}.engine.OF/(1+opt.stage{i}.engine.OF)*mProp;

                volumeOx = mOx/opt.stage{i}.engine.oxDens;
                volumeFuel = mFuel/opt.stage{i}.engine.fuelDens;

                volumeTankLOX = volumeOx*1.055;
                volumeTankFuel = volumeFuel*1.055;
                pressurantMass = 340000*(volumeTankLOX + volumeTankFuel)/(8314/4*200);
                pressurantTankMass = 2*pressurantMass;
                fun1 =@(x) [ volumeTankFuel-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); volumeTankLOX-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
                x01 = [2;10;10];
                
                sol1 = fsolve(fun1,x01,settings.fsolveOptions);

                areaTankLOX = 2*sol1(1)*pi*sol1(3) + 4*pi*sol1(1)^2;
                areaTankFuel = 2*sol1(1)*pi*sol1(2) + 4*pi*sol1(1)^2;

                mer.stage{i}.tankMassOx = 12.2 * volumeTankLOX + 255.2;
                mer.stage{i}.tankMassFuel= 12.16 * volumeTankFuel + 255.2;
                mer.stage{i}.pressurant = pressurantTankMass + pressurantMass;
                mer.stage{i}.cryoInsuOx = 1.12 * areaTankLOX;
                mer.stage{i}.cryoInsuFuel = 0*areaTankFuel;
                mer.stage{i}.thrustFrame = 2.55e-4*nEngine * thrustValue;
                mer.stage{i}.engineWeight = nEngine*opt.stage{i}.engine.weight;
                mer.stage{i}.tvc = 0.12*nEngine*opt.stage{i}.engine.weight;
                mer.stage{i}.interStage = 13.7*sol1(1)^2*pi;
            elseif strcmp(couple, 'CH4-LOX')

                mFuel = 1/(1+opt.stage{i}.engine.OF)*mProp;
                mOx = opt.stage{i}.engine.OF/(1+opt.stage{i}.engine.OF)*mProp;

                volumeOx = mOx/opt.stage{i}.engine.oxDens;
                volumeFuel = mFuel/opt.stage{i}.engine.fuelDens;

                volumeTankLOX = volumeOx*1.055;
                volumeTankFuel = volumeFuel*1.055;

                pressurantMass = 340000*(volumeTankLOX + volumeTankFuel)/(8314/4*200);
                pressurantTankMass = 2*pressurantMass;

                fun1 =@(x) [ volumeTankFuel-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); volumeTankLOX-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
                x01 = [2;10;10];
                sol1 = fsolve(fun1,x01,settings.fsolveOptions);

                areaTankLOX = 2*sol1(1)*pi*sol1(3) + 4*pi*sol1(1)^2;
                areaTankFuel = 2*sol1(1)*pi*sol1(2) + 4*pi*sol1(1)^2;

                mer.stage{i}.tankMassOx = 12.2 * volumeTankLOX + 255.2;
                mer.stage{i}.tankMassFuel = 12.2 * volumeTankFuel + 255.2;
                mer.stage{i}.pressurant = pressurantTankMass + pressurantMass;
                mer.stage{i}.cryoInsuOx = 1.12 * areaTankLOX;
                mer.stage{i}.cryoInsuFuel = 1.12 * areaTankFuel;
                mer.stage{i}.thrustFrame = 2.55e-4 *nEngine*thrustValue;
                mer.stage{i}.engineWeight = nEngine*opt.stage{i}.engine.weight;
                mer.stage{i}.tvc = 0.12*nEngine*opt.stage{i}.engine.weight;   
                mer.stage{i}.interStage = 13.7*sol1(1)^2*pi;
            else

             
                mFuel = 1/(1+opt.stage{i}.engine.OF)*mProp;
                mOx = opt.stage{i}.engine.OF/(1+opt.stage{i}.engine.OF)*mProp;

                volumeOx = mOx/opt.stage{i}.engine.oxDens;
                volumeFuel = mFuel/opt.stage{i}.engine.fuelDens;

                volumeTankLOX = volumeOx*1.055;
                volumeTankFuel = volumeFuel*1.055;

                pressurantMass = 340000*(volumeTankLOX + volumeTankFuel)/(8314/4*200);
                pressurantTankMass = 2*pressurantMass;

                fun1 =@(x) [ volumeTankFuel-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); volumeTankLOX-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
                x01 = [2;10;10];
                sol1 = fsolve(fun1,x01,settings.fsolveOptions);

                areaTankLOX = 2*sol1(1)*pi*sol1(3) + 4*pi*sol1(1)^2;
                areaTankFuel = 2*sol1(1)*pi*sol1(2) + 4*pi*sol1(1)^2;

                mer.stage{i}.tankMassOx = 12.2 * volumeTankLOX + 255.2;
                mer.stage{i}.tankMassFuel = 9.08 * volumeTankFuel + 100.09;
                mer.stage{i}.pressurant = pressurantTankMass + pressurantMass;
                mer.stage{i}.cryoInsuOx = 1.12 * areaTankLOX;
                mer.stage{i}.cryoInsuFuel = 2.88 * areaTankFuel;
                mer.stage{i}.thrustFrame = 2.55e-4*nEngine*thrustValue;
                mer.stage{i}.engineWeight = nEngine*opt.stage{i}.engine.weight;
                mer.stage{i}.tvc = 0.12*nEngine*opt.stage{i}.engine.weight;
                mer.stage{i}.interStage = 13.7*sol1(1)^2*pi;                
            end

            mStructure = mer.stage{i}.tankMassOx + mer.stage{i}.tankMassFuel + ...
                         mer.stage{i}.cryoInsuOx + mer.stage{i}.cryoInsuFuel + ...
                         mer.stage{i}.thrustFrame + mer.stage{i}.engineWeight + ...
                         mer.stage{i}.avionics + mer.stage{i}.wiring + mer.stage{i}.tvc +  ...
                         mer.stage{i}.battery + mer.stage{i}.pressurant + mer.stage{i}.interStage;

            mStage = mStructure + mProp ;
            mStack = mStack + mStage;
            twRatio = nEngine*thrustValue/mStack/9.81;
        end

    staging{i}.mStage = mStage;
    staging{i}.mProp = mProp;
    staging{i}.mStruct = mStructure;

    end



end
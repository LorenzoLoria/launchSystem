function [mer,staging] = initialMassEstimation(mission,opt)

nStages = opt.nStages;
lOverD = 5;
    avionics =350; % 10 * (M0)^0.361; the previous was overyly conservative as per il libro che ho ciulato
    wiring = 10;

    if nStages == 1
        mer.stage{1}.avionics = avionics;
        mer.stage{1}.wiring = wiring;
    elseif nStages == 2
        mer.stage{1}.avionics = 0.2*avionics;
        mer.stage{2}.avionics = 0.8*avionics;
    
        mer.stage{1}.wiring = 0.2*wiring;
        mer.stage{2}.wiring = 0.8*wiring;
    else
        mer.stage{1}.avionics = 0.1*avionics;
        mer.stage{2}.avionics = 0.1*avionics;
        mer.stage{3}.avionics = 0.8*avionics;
    
        mer.stage{1}.wiring = 0.1*wiring;
        mer.stage{2}.wiring = 0.1*wiring;
        mer.stage{3}.wiring = 0.8*wiring;
    end
    
    mStack = mission.capsule.weigth;

    for i = nStages:-1:1
        couple = opt.stage{i}.engine.couple;

        nEngine = 0;
        twRatio = 0;
        mStage = 0;
        while twRatio < 1
            mStack = mStack-mStage;
            nEngine = nEngine+1;
            mProp = nEngine * opt.stage{i}.engine.thrust/9.81 * opt.stage{i}.percentage;

            if strcmp(couple, 'RP1-LOX')

                mFuel = 1/(1+opt.stage{i}.engine.OF)*mProp;
                mOx = opt.stage{i}.engine.OF/(1+opt.stage{i}.engine.OF)*mProp;

                volumeOx = mOx/opt.stage{i}.engine.oxDens;
                volumeFuel = mFuel/opt.stage{i}.engine.fuelDens;

                volumeTankLOX = volumeOx*1.055;
                volumeTankFuel = volumeFuel*1.055;

                fun1 =@(x) [ volumeTankFuel-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); volumeTankLOX-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
                x01 = [2;10;10];
                sol1 = fsolve(fun1,x01);

                areaTankLOX = sol1(1)^2*pi;

                mer.stage{i}.tankMassOx = 12.2 * volumeTankLOX + 255.2;
                mer.stage{i}.tankMassFuel= 12.16 * volumeTankFuel;
                mer.stage{i}.cryoInsuOx = 1.12 * areaTankLOX;
                mer.stage{i}.cryoInsuFuel = 0;
                mer.stage{i}.thrustFrame = 2.55e-4*nEngine * opt.stage{i}.engine.thrust;
                mer.stage{i}.engineWeight = nEngine*opt.stage{i}.engine.weight;

            elseif strcmp(couple, 'CH4-LOX')

                mFuel = 1/(1+opt.stage{i}.engine.OF)*mProp;
                mOx = opt.stage{i}.engine.OF/(1+opt.stage{i}.engine.OF)*mProp;

                volumeOx = mOx/opt.stage{i}.engine.oxDens;
                volumeFuel = mFuel/opt.stage{i}.engine.fuelDens;

                volumeTankLOX = volumeOx*1.055;
                volumeTankFuel = volumeFuel*1.055;

                fun1 =@(x) [ volumeTankFuel-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); volumeTankLOX-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
                x01 = [2;10;10];
                sol1 = fsolve(fun1,x01);

                areaTankLOX = sol1(1)^2*pi;
                areaTankFuel = areaTankLOX;

                mer.stage{i}.tankMassOx = 12.2 * volumeTankLOX + 255.2;
                mer.stage{i}.tankMassFuel = 12.2 * volumeTankFuel + 255.2;
                mer.stage{i}.cryoInsuOx = 1.12 * areaTankLOX;
                mer.stage{i}.cryoInsuFuel = 1.12 * areaTankFuel;
                mer.stage{i}.thrustFrame = 2.55e-4 *nEngine*opt.stage{i}.engine.thrust;
                mer.stage{i}.engineWeight = nEngine*opt.stage{i}.engine.weight;

            else


                mFuel = 1/(1+opt.stage{i}.engine.OF)*mProp;
                mOx = opt.stage{i}.engine.OF/(1+opt.stage{i}.engine.OF)*mProp;

                volumeOx = mOx/opt.stage{i}.engine.oxDens;
                volumeFuel = mFuel/opt.stage{i}.engine.fuelDens;

                volumeTankLOX = volumeOx*1.055;
                volumeTankFuel = volumeFuel*1.055;

                fun1 =@(x) [ volumeTankFuel-4/3*pi*x(1)^3 - pi*x(1)^2*x(2); volumeTankLOX-4/3*pi*x(1)^3 - pi*x(1)^2*x(3); 0.8*lOverD- (4*x(1) + x(2)+ x(3))/(2*x(1)) ];
                x01 = [2;10;10];
                sol1 = fsolve(fun1,x01);

                areaTankLOX = sol1(1)^2*pi;
                areaTankFuel = areaTankLOX;

                mer.stage{i}.tankMassOx = 12.2 * volumeTankLOX + 255.2;
                mer.stage{i}.tankMassFuel = 9.08 * volumeTankFuel + 100.09;
                mer.stage{i}.cryoInsuOx = 1.12 * areaTankLOX;
                mer.stage{i}.cryoInsuFuel = 2.88 * areaTankFuel;
                mer.stage{i}.thrustFrame = 2.55e-4*nEngine*opt.stage{i}.engine.thrust;
                mer.stage{i}.engineWeight = nEngine*opt.stage{i}.engine.weight;

            end

            mStructure = mer.stage{i}.tankMassOx + mer.stage{i}.tankMassFuel + ...
                         mer.stage{i}.cryoInsuOx + mer.stage{i}.cryoInsuFuel + ...
                         mer.stage{i}.thrustFrame + mer.stage{i}.engineWeight + ...
                         mer.stage{i}.avionics + mer.stage{i}.wiring ;

            mStage = mStructure + mProp ;
            mStack = mStack + mStage;
            twRatio = nEngine*opt.stage{i}.engine.thrust/mStack/9.81;
        end

    staging{i}.mStage = mStage;
    staging{i}.mProp = mProp;
    staging{i}.mStruct = mStructure;

    end



end
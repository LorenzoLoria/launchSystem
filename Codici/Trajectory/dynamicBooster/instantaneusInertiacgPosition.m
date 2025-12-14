function [Isym,Iy,mass,cg] = instantaneusInertiacgPosition(mission,configuration,mer,staging,stageNumber) 


fuelPercentage = 0.05;
tankMass = mer.stage{stageNumber}.tankMassOx + mer.stage{stageNumber}.tankMassOx + mer.stage{stageNumber}.cryoInsuOx + mer.stage{stageNumber}.cryoInsuFuel;
tankRadius = configuration.geometry.stage{stageNumber}.radius;
tankHeightTotal = configuration.stage{stageNumber}.fuelTankH + configuration.stage{stageNumber}.oxTankH;
xcgTanks = tankRadius + tankHeightTotal/2;
[Isym,Iy] = tankInertia(tankMass, tankRadius, tankHeightTotal);
propMass = fuelPercentage * staging{stageNumber}.mProp;
xcgProp =(1-fuelPercentage)*tankHeightTotal; %approximation
remainingStructuralMass = staging{stageNumber}.mStruct -  tankMass - mer.stage{stageNumber}.interStage;
xcgRemainings = tankHeightTotal + configuration.geometry.stage{stageNumber}.thrustFrame/2;
xCGvec = [xcgTanks,xcgProp,xcgRemainings];
massVec = [tankMass,propMass,remainingStructuralMass];

cg = centerOfGravity(massVec,xCGvec);
mass = sum(massVec);

Iy = Iy + tankMass*(cg-xcgTanks)^2 + propMass*(cg-xcgProp)^2 + remainingStructuralMass * (cg-xcgRemainings)^2;

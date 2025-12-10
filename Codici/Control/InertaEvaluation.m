function [inertia] = InertaEvaluation(mission, configuration, launcher, totalStageNumber )
% INERTIAEVALUATION evaluate the inertia matrix of launcher + capsule
% Capsule data (assumption of conical shape) -- > from Sforza Book
mCapsule = mission.capsule.weight ;
rCapsule = mission.capsule.radius ;
hCapsule = mission.capsule.height ;

inertia.capsule = 3 * mCapsule * hCapsule^2 / 80 * ( 1 + 4 * (rCapsule / hCapsule)^2) ;

% Stage number
stageNumber = launcher(1) ;
inertia.stage = 0 ;
inertia.interStage = 0 ;

% launcher Data
for ii = totalStageNumber : -1 : stageNumber
   hStageIter = configuration.geometry.stage{ii}.length;
   rStageIter = configuration.geometry.stage{ii}.radius;
   mStageIter = configuration.stage{ii}.mStage ;

   hInterstageIter = configuration.geometry.stage{ii}.interstage.lenght ;
   rInterstageIter = configuration.geometry.stage{ii-1}.radius;
   mInterstageIter = mer.stage{ii}.interStage ;
   inertia.stage(ii) = mStageIter * hStageIter^2 / 12 * (1 + 3 * (rStageIter / hStageIter)^2 ) ;
   inertia.stage = inertia.stage + inertia.stage(ii) ;

   inertia.interstage(ii) = 3 * mInterstageIter * (rStageIter^5 - rInterstageIter^5) / 10 / (rStageIter^3 - rInterstageIter^3)...
       - mInterstageIter * hInterstageIter^2 / 16 * ( (rStageIter^2 + 2 * ...
        rStageIter * rInterstageIter + 3 * rInterstageIter^2) / (rStageIter^2 + 2 * rStageIter * rInterstageIter + rInterstageIter^2)) ;
   inertia.interStage = inertia.interStage + inertia.interStage(ii) ;
end

inertia = inertia.capsule + inertia.stage + inertia.interStage ;

end
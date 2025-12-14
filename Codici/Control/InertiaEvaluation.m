function [inertia,Xcg] = InertiaEvaluation(mission, configuration, mer, stageNumber, m, launcher)
% INERTIAEVALUATION evaluate the inertia matrix of launcher + capsule

inertia = zeros(3, 3);
Xcg = computeXcgGlobal(mission, configuration, launcher, mer, m, stageNumber);

% CAPSULE data (assumption of conical shape) --> from Sforza Book
mCapsule = mission.capsule.weight;
rCapsule = mission.capsule.radius;
hCapsule = mission.capsule.height;
inertiaShiftCaps = Xcg - 3/4 * hCapsule; % cg capsule is at 3/4 h from the nose
r_com_capsule = [0; 0; inertiaShiftCaps];

% Local inertia matrix for capsule (about its COM)
Ixx_capsule = (3/80) * mCapsule * (4 * rCapsule^2 + hCapsule^2);
Izz_capsule = (3/10) * mCapsule * rCapsule^2;
I_local_capsule = diag([Ixx_capsule, Ixx_capsule, Izz_capsule ]) ;

inertia = inertia + translate_inertia(I_local_capsule, mCapsule, r_com_capsule) ;

%% STAGES AND INTERSTAGES

cumulative_height = hCapsule;

for ii = launcher(1):-1:stageNumber
    % Stage data
    hStageIter = configuration.geometry.stage{ii}.tanksLength;
    rStageIter = configuration.geometry.stage{ii}.radius;

    if ii == stageNumber 

        if stageNumber == 1
            mStageIter =  configuration.stage{ii}.mStage - (configuration.totalMass - m) ;
        else
            totalMassIter = configuration.totalMass - configuration.stage{stageNumber - 1}.mStage ;
            mStageIter =  configuration.stage{ii}.mStage - (totalMassIter - m) ;
        end

    else
        mStageIter = configuration.stage{ii}.mStage ;
    end

    hInterstageIter = configuration.geometry.stage{ii}.interstage.length;
    
    if ii == launcher(1)
        rBottom = configuration.geometry.stage{ii}.radius; % bottom = lower stage
        rTop = mission.capsule.radius;
    else
        rBottom = configuration.geometry.stage{ii}.radius; % bottom = lower stage
        rTop = configuration.geometry.stage{ii+1}.radius;
    end
    mInterstageIter = mer.stage{ii}.interStage;

    % Center of mass for frustum
    z_com_local = hInterstageIter - hInterstageIter * (rBottom^2 + 2*rBottom*rTop + 3*rTop^2) / ...
                  (4 * (rBottom^2 + rBottom*rTop + rTop^2));

    % Moments of inertia for frustum (approximation)
    Izz_interstage = (3/10) * mInterstageIter * (rBottom^5 - rTop^5) / (rBottom^3 - rTop^3);
    Ixx_interstage = (3/20) * mInterstageIter * (rBottom^5 - rTop^5) / (rBottom^3 - rTop^3);
    I_local_inter = diag([Ixx_interstage, Ixx_interstage, Izz_interstage]);

    z_interstage = Xcg - (cumulative_height + z_com_local);
    r_com_interstage = [0; 0; z_interstage];

    inertia = inertia + translate_inertia(I_local_inter, mInterstageIter, r_com_interstage);
    cumulative_height = cumulative_height + hInterstageIter;

    % Local inertia matrix for stage (cylinder)
    Ixx_stage = (1/12) * mStageIter * (3 * rStageIter^2 + hStageIter^2);
    Izz_stage = (1/2) * mStageIter * rStageIter^2;
    I_local_stage = diag([Ixx_stage, Ixx_stage, Izz_stage]);

    z_stage = Xcg - (cumulative_height + hStageIter / 2); % z-coordinate of the stage's COM
    r_com_stage = [0; 0; z_stage];

    inertia = inertia + translate_inertia(I_local_stage, mStageIter, r_com_stage);

    cumulative_height = cumulative_height + hStageIter;



end

end




function I_global = translate_inertia(I_local, m, r_com)
    x = r_com(1);
    y = r_com(2); 
    z = r_com(3);
    
    I_global = I_local + m * [ ...
        y^2 + z^2, -x*y, -x*z; ...
        -x*y, x^2 + z^2, -y*z; ...
        -x*z, -y*z, x^2 + y^2 ...
    ];
end

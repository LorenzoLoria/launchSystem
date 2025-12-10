function [inertia] = InertaEvaluation(mission, configuration, mer, launcher, totalStageNumber)
% INERTIAEVALUATION evaluate the inertia matrix of launcher + capsule

inertia = zeros(3, 3);

%% CAPSULE
% Capsule data (assumption of conical shape) --> from Sforza Book
mCapsule = mission.capsule.weight;
rCapsule = mission.capsule.radius;
hCapsule = mission.capsule.height;
z_capsule = hCapsule / 4; % COM of a cone is at h/4 from the base
r_com_capsule = [0; 0; z_capsule];


% Local inertia matrix for capsule (about its COM)
Ixx_capsule = (3/80) * mCapsule * (4 * rCapsule^2 + hCapsule^2);
Izz_capsule = (3/10) * mCapsule * rCapsule^2;

%% STAGES AND INTERSTAGES
stageNumber = launcher(1);

cumulative_height = hCapsule;

for ii = totalStageNumber:-1:stageNumber
    % Stage data
    hStageIter = configuration.geometry.stage{ii}.length;
    rStageIter = configuration.geometry.stage{ii}.radius;
    mStageIter = configuration.stage{ii}.mStage;

    % Local inertia matrix for stage (cylinder)
    Ixx_stage = (1/12) * mStageIter * (3 * rStageIter^2 + hStageIter^2);
    Izz_stage = (1/2) * mStageIter * rStageIter^2;
    I_local_stage = diag([Ixx_stage, Ixx_stage, Izz_stage]);

    z_stage = cumulative_height + hStageIter / 2; % z-coordinate of the stage's COM
    r_com_stage = [0; 0; z_stage];

    inertia = inertia + translate_inertia(I_local_stage, mStageIter, r_com_stage);

    cumulative_height = cumulative_height + hStageIter;

    % Interstage data
    if ii > 1
        hInterstageIter = configuration.geometry.stage{ii}.interstage.length;
        rBottom = configuration.geometry.stage{ii-1}.radius; % bottom = lower stage
        rTop = configuration.geometry.stage{ii}.radius;
        mInterstageIter = mer.stage{ii}.interStage;

        % Center of mass for frustum
        z_com_local = hInterstageIter * (rBottom^2 + 2*rBottom*rTop + 3*rTop^2) / ...
                      (4 * (rBottom^2 + rBottom*rTop + rTop^2));

        % Moments of inertia for frustum (approximation)
        Izz_interstage = (3/10) * mInterstageIter * (rBottom^5 - rTop^5) / (rBottom^3 - rTop^3);
        Ixx_interstage = (3/20) * mInterstageIter * (rBottom^5 - rTop^5) / (rBottom^3 - rTop^3);
        I_local_inter = diag([Ixx_interstage, Ixx_interstage, Izz_interstage]);

        z_interstage = cumulative_height + z_com_local;
        r_com_interstage = [0; 0; z_interstage];

        inertia = inertia + translate_inertia(I_local_inter, mInterstageIter, r_com_interstage);
        cumulative_height = cumulative_height + hInterstageIter;
    end
end

end

function I_global = translate_inertia(I_local, m, r_com)
    x = r_com(1); y = r_com(2); z = r_com(3);
    I_global = I_local + m * [ ...
        y^2 + z^2, -x*y, -x*z; ...
        -x*y, x^2 + z^2, -y*z; ...
        -x*z, -y*z, x^2 + y^2 ...
    ];
end

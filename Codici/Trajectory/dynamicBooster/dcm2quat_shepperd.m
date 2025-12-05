function q = dcm2quat_shepperd(R)
% DCM2QUAT_SHEPPERD Converts a Rotation Matrix (DCM) to a Quaternion
% using Shepperd's Algorithm for numerical stability.
%
%   q = dcm2quat_shepperd(R)
%
%   INPUTS:
%       R - 3x3 Rotation Matrix (Direction Cosine Matrix)
%           Transforms from Body frame to Inertial frame (or vice versa).
%           Assumed to be orthonormal.
%
%   OUTPUTS:
%       q - 1x4 Quaternion vector [w, x, y, z] (Scalar First convention)
%           Where q = w + xi + yj + zk
%
%   REFERENCES:
%       Shepperd, S.W., "Quaternion from Rotation Matrix," 
%       Journal of Guidance and Control, Vol. 1, No. 3, 1978.

    % Extract diagonal and trace
    r11 = R(1,1); r12 = R(1,2); r13 = R(1,3);
    r21 = R(2,1); r22 = R(2,2); r23 = R(2,3);
    r31 = R(3,1); r32 = R(3,2); r33 = R(3,3);
    
    tr = r11 + r22 + r33;

    % Shepperd's Logic: 
    % Find the maximum of the four possible denominators to ensure stability.
    % The four decision variables correspond to the squares of w, x, y, z.
    
    d0 = tr;             % Associated with w (scalar)
    d1 = r11;            % Associated with x
    d2 = r22;            % Associated with y
    d3 = r33;            % Associated with z
    
    % Find which case to use
    [~, best_index] = max([d0, d1, d2, d3]);
    
    switch best_index
        case 1 % Trace is the dominant factor (Standard Case)
            % Formula: 4w^2 = 1 + trace
            S = sqrt(1.0 + tr) * 2; % S = 4w
            w = 0.25 * S;
            x = (r32 - r23) / S;
            y = (r13 - r31) / S;
            z = (r21 - r12) / S;
            
        case 2 % R(1,1) is dominant
            % Formula: 4x^2 = 1 + r11 - r22 - r33
            S = sqrt(1.0 + r11 - r22 - r33) * 2; % S = 4x
            w = (r32 - r23) / S;
            x = 0.25 * S;
            y = (r12 + r21) / S;
            z = (r13 + r31) / S;
            
        case 3 % R(2,2) is dominant
            % Formula: 4y^2 = 1 - r11 + r22 - r33
            S = sqrt(1.0 + r22 - r11 - r33) * 2; % S = 4y
            w = (r13 - r31) / S;
            x = (r12 + r21) / S;
            y = 0.25 * S;
            z = (r23 + r32) / S;
            
        case 4 % R(3,3) is dominant
            % Formula: 4z^2 = 1 - r11 - r22 + r33
            S = sqrt(1.0 + r33 - r11 - r22) * 2; % S = 4z
            w = (r21 - r12) / S;
            x = (r13 + r31) / S;
            y = (r23 + r32) / S;
            z = 0.25 * S;
    end
    
    % Pack into vector
    q = [w, x, y, z];
    
    % Enforce unit norm (correction for numerical errors in S calculation)
    q = q / norm(q);

end
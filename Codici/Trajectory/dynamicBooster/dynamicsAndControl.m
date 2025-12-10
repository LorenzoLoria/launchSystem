function [dxdt] = dynamicsAndControl(t,x,targetPos,targetVel,gains,initialPos)
if ~iscolumn(x)
    x=x';
end
if ~iscolumn(gains)
    gains=gains';
end

posIR = x(1:3) ;
velIR = x(4:6) ;
qCurr = x(7:10) ;
omegaBR = x(11:13) ;




terminalHeight = 1 ;

if posIR(3) > terminalHeight

    heightWeight = posIR(3) / initialPos(3) ;

    kPos = heightWeight * gains(1:3) ;
    kVel = (1.7 - heightWeight) * gains(4:6) ;
    KpAtt = gains(7) ;
    KdAtt = gains(8) ;

else
    
    kPos = 0 ;
    kVel = 3 * gains(4:6) ; 
    KpAtt = gains(7) ;
    KdAtt = gains(8) ;

end





%Geometrical data
g = [0 0 -9.81]' ;
m = 20 ; 
I = diag([10 10 1]) ;
thrustArm = 2 ;
xCpCg = 2 ;
S = 4 * 1 ;
cd = 0.2 ; 
rho = 1.225 ; 

%Max thrust, gimbal and attitude
Tmax = 400 ;
deltaGimbalMax = 20 ;
maxTiltAngle = 20 ;

%Rotation matrices
Rbi = [ 1-2*(qCurr(3)^2+qCurr(4)^2),  2*(qCurr(2)*qCurr(3)-qCurr(4)*qCurr(1)),    2*(qCurr(2)*qCurr(4)+qCurr(3)*qCurr(1));
      2*(qCurr(2)*qCurr(3)+qCurr(4)*qCurr(1)),    1-2*(qCurr(2)^2+qCurr(4)^2),  2*(qCurr(3)*qCurr(4)-qCurr(2)*qCurr(1));
      2*(qCurr(2)*qCurr(4)-qCurr(3)*qCurr(1)),    2*(qCurr(3)*qCurr(4)+qCurr(2)*qCurr(1)),    1-2*(qCurr(2)^2+qCurr(3)^2) ];

Rib = Rbi' ;



%Guidance control law
kAcc = 0 ;
if posIR(3)<2
    kAcc = 5 ;
end

aCommand = kPos .* (targetPos - posIR) + kVel .* (targetVel - velIR) - g ;



% Desired attitude retrieval
if posIR(3) > terminalHeight
    zDes = aCommand / norm(aCommand) ;
else
    zDes = [0 0 1]' ;
end
% zDes = aCommand / norm(aCommand) ;

if acosd(zDes(3)/norm(zDes)) > maxTiltAngle
    zDes(3) = norm(zDes(1:2)) / tand(maxTiltAngle) ;
    zDes = zDes / norm(zDes) ;
end

yDes = [0 1 0]' ;
if abs(dot(zDes, yDes)) > 0.99 % Singularity check
        yDes = [1; 0; 0];
end
xDes = cross(yDes,zDes) ;
yDes = cross(zDes,xDes) ;

Rdes = [xDes yDes zDes] ;
qDes = dcm2quat_shepperd(Rdes)' ;

qDesMat =   [qDes(1) -qDes(2) -qDes(3) -qDes(4);...
             qDes(2) qDes(1) -qDes(4)  qDes(3);...
             qDes(3) qDes(4) qDes(1)  -qDes(2);...
             qDes(4) -qDes(3) qDes(2)  qDes(1)];

% Attitude control law
qCurrConj = [q(1) ; - q(2:4)];
qErrIR = qDesMat * qCurrConj ;

if qErrIR(1) < 0    % short path rotation
qErrIR = -qErrIR;
end

qErrBR = [qErrIR(1) ; IRFtoBRF * qErrIR(2:4)] ;

torqueCommand = KpAtt * qErrBR(2:4) - KdAtt * omegaBR ; 


% Drag computation

dragIR = - 0.5 * rho * norm(velIR) * velIR * S * cd ;
dragBR = Rib * dragIR ;

qDot = 0.5 * [-qCurr(2) -qCurr(3) -qCurr(4) ;...
              qCurr(1) -qCurr(4) qCurr(3) ;...
              qCurr(4) qCurr(1) -qCurr(2) ;...
              -qCurr(3) qCurr(2) qCurr(1)] ; 

qDot = qDot * omegaBR ;

% Finding gimbal angle and updating acceleration



% xCurr = Rbi(:,1) / norm(Rbi(:,1)) ;
% yCurr = Rbi(:,2) / norm(Rbi(:,2)) ;
% zCurr = Rbi(:,3) / norm(Rbi(:,3)) ;


thrustGuidanceBR = Rib * (m * aCommand) ;
thrustBR = [-torqueCommand(2)/thrustArm ; torqueCommand(1)/thrustArm ; thrustGuidanceBR(3)] ;

%Checks on maxThrust and gimbal angle
if acosd(thrustBR(3)/norm(thrustBR)) > deltaGimbalMax 
    thrustBR(3) = norm(thrustBR(1:2)) / tand(deltaGimbalMax) ;
end

if norm(thrustBR) > Tmax
    thrustBR = thrustBR / norm(thrustBR) * Tmax ;
end

torqueBR = cross([0 0 -thrustArm]', thrustBR) + cross([0 0 xCpCg]', dragBR);
omegaDot = inv(I) * (torqueBR - cross(omegaBR,I*omegaBR)) ;



thrustIR = Rbi * thrustBR ;

totalForceIR = thrustIR + dragIR ;

%Variables update: 
%velocity is updated as point mass with the actual thrust applied, not
%accelration commanded by guidance


dxdt(1:3) = velIR ;
dxdt(4:6) = totalForceIR / m + g; 
dxdt(7:10) = qDot ;
dxdt(11:13) = omegaDot ;

dxdt = dxdt' ;

if any(isnan(dxdt))
    keyboard
end



end
function [output] = findGains(x,targetPos,targetVel,y0,tSpan,r0)


gains = x ;
if ~iscolumn(gains)
    gains = gains' ;
end


options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@touchDown, 'MinStep', 1e-6);
[t, sol] = ode45(@(t,x) dynamicsAndControl(t,x,targetPos,targetVel,gains,r0), tSpan, y0,options);




output = sum(abs(sol(end,1:3)'-targetPos)) + sum(abs(sol(end,4:6)'-targetVel)) ;






end
function [Fcmd, Mcmd] = controller(err_pos, err_vel, err_theta, err_omega, KPos, KVel, KpAtt, KdAtt, m, g)


Fcmd = KPos * err_pos + KVel * err_vel + m * g;


Mcmd = KpAtt * err_theta + KdAtt * err_omega;

end

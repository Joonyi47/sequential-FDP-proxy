function [pos,vel] = checkBoundaries(pos,vel,maxvel,var_max,var_min)
global N Np

MAXLIM   = repmat(var_max, 1, Np)';
MINLIM   = repmat(var_min, 1, Np)';
MAXVEL   = repmat(maxvel(:)',Np,1);
MINVEL   = repmat(-maxvel(:)',Np,1);

vel(vel>MAXVEL) = MAXVEL(vel>MAXVEL);
vel(vel<MINVEL) = MINVEL(vel<MINVEL);
% vel(pos>MAXLIM) = (-1).*vel(pos>MAXLIM);
pos(pos>MAXLIM) = MAXLIM(pos>MAXLIM);
% vel(pos<MINLIM) = (-1).*vel(pos<MINLIM);
pos(pos<MINLIM) = MINLIM(pos<MINLIM);
end
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
pos(:,1:2*N)    = round(pos(:,1:2*N));

for j = 1:Np
    tmp = [];
    for i = 1:N
        tmp = [tmp, DefLimit(pos(j, 2*i-1:2*i))];
    end
    MINLIM(j,1:2*N) = tmp(1,:);
    MAXLIM(j,1:2*N) = tmp(2,:);
end
pos(pos>MAXLIM) = MAXLIM(pos>MAXLIM);
pos(pos<MINLIM) = MINLIM(pos<MINLIM);
end
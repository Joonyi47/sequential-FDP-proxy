function [lim] = DefLimit(pos)
global ACTIVE nx ny

active = reshape(ACTIVE, nx*ny, 2);
active = active(:,1).*active(:,2);
active = reshape(active, nx, ny)';

pos = reshape(pos, 2, 1);

tmp = find(active(pos(2), :) == 1);
minlim_x = tmp(1);
maxlim_x = tmp(end);

tmp = find(active(:, pos(1)) == 1);
minlim_y = tmp(1);
maxlim_y = tmp(end);

[tmpx, Ix] = min([minlim_x - pos(1), maxlim_x - pos(1)]);
[tmpy, Iy] = min([minlim_y - pos(2), maxlim_y - pos(2)]);
if Ix == 1
    if Iy == 1
        if tmpx <= tmpy
            minlim_y = pos(2);
        else
            minlim_x = pos(1);
        end
    else
        if tmpx <= tmpy
            maxlim_y = pos(2);
        else
            minlim_x = pos(1);
        end
    end
else
    if Iy == 1
        if tmpx <= tmpy
            minlim_y = pos(2);
        else
            maxlim_y = pos(1);
        end
    else
        if tmpx <= tmpy
            maxlim_y = pos(2);
        else
            maxlim_x = pos(1);
        end
    end
end

minlim = [minlim_x, minlim_y];
maxlim = [maxlim_x, maxlim_y];

lim = [minlim; maxlim];
end
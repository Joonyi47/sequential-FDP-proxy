function [viot] = constViolation(type)
%%% pos ´Â ¿­º¤ÅÍ = {x1, y1, x2, y2, ... }

global N Nstep dstep pmax ...
       opt_type

type = reshape(type, 4, N);

vio = zeros(N, Nstep);
for i = 1:N
    vio(i, (type(1,i) == 0) + ceil(Nstep * type(1,i))) = 1;
end

if sum((opt_type >= 1/3) .* vio(:,1)') == 0
    viot = 100;
else    
    vio = sum(vio, 1);
    if vio(1) == 0
        viot = 100;
    else
        vio = vio - 3;
        viot = sum(vio(1,vio > 0));
    end
end
end
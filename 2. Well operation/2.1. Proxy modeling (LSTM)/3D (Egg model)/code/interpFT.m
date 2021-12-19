function [output] = interpFT(input)
global pmax dstep Nstep tstep

Ndata = length(input);
output = cell(Ndata, 1);

for j = 1:Ndata
    for p = 1:size(input{1},1)
        
        tmp = interp1([0:pmax/tstep/Nstep:pmax/tstep], [0, input{j}(p,:)], [0:1:(pmax/tstep)], 'pchip');
        tmp2 = tmp(1,2:end);
        tmp2(tmp2<0) = 0;
        
        output{j}(p,:) = tmp2;
    end
end

end
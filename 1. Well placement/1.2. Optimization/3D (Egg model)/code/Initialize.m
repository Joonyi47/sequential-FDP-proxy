function [pos] = Initialize(type)
global N Np nx ny ACTIVE area

active = reshape(ACTIVE, nx*ny, 2);
active = active(:,1).*active(:,2);
active = reshape(reshape(active, nx, ny)', nx*ny,1);
idx_act = find(active == 1);
radius = sqrt(43560*area/pi);

thres = 0.5*radius^2*ones(Np,1);
thres(randperm(Np,3),1) = 0.001;

% while length(find(vio == 0)) < Np
    pos = []; vio = [];
    
    h = waitbar(0,'Please wait...');
    for i = 1:Np
        
        tmp3 = radius^2;
        
        while tmp3 > thres(i,1)
            
            tmp  = randperm(length(find(active == 1)), N);
            tmp  = idx_act(tmp);
            tmp1 = reshape([ceil(tmp'/ny); mod(tmp', ny) + (mod(tmp',ny) == 0)*ny], 1, 2*N);
            
            if isempty(type)
                tmp2 = [1, -1, -1 + 2*rand(1,N-2)];
            else
                tmp2 = type;
            end
            
            tmp3 = constViolation(tmp1, tmp2);
        end
        
        pos = [pos; [tmp1, tmp2]];
        vio = [vio; tmp3];
        
        waitbar(i / Np)
        
    end
%     disp(int2str(length(find(vio < 0.5*radius^2))));
% end
close(h)
end
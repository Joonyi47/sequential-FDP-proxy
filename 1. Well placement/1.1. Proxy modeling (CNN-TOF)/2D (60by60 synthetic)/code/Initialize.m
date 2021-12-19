function [pos] = Initialize(type)
global N Np nx ny area

radius = sqrt(43560*area/pi);

thres = 0.5*radius^2*ones(Np,1);
thres(randperm(Np,0.3*Np),1) = 0.001;


pos = []; vio = [];

h = waitbar(0,'Please wait...');
for i = 1:Np
    
    tmp3 = (1)*radius^2;
    
    while tmp3 > thres(i,1)
        
%         tmp = reshape([randi([1 nx], 1, N); randi([1 ny], 1, N)], 1, 2*N);
%         tmp = round((ny-1)*lhsdesign(Np,2*N)+1);
        tmp = randperm(nx*ny, N);
        tmp = reshape([ceil(tmp'/ny); mod(tmp', ny) + (mod(tmp',ny) == 0)*ny], 1, 2*N);

        if isempty(type)
            
%             tmp2 = [1, 1, round(rand(1, N-3)), 0];
            tmp2 = [1/3+2/3*rand(1), -1 + 2*rand(1,N-1)];
%             tmp2 = [1/3 + 2/3*rand(1,4), -1 + 2/3*rand(1,4), -1 + 2*rand(1,N-8)];
            tmp3 = constViolation(tmp, tmp2);
            
        else
            
            tmp2 = [];
            tmp3 = constViolation(tmp, type);
            
        end
                
    end
        
    pos = [pos; [tmp, tmp2]];
    vio = [vio; tmp3];
    
    waitbar(i / Np)
    
end
close(h)
end


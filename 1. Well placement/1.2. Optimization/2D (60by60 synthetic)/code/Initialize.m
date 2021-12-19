function [pos] = Initialize(type)
global N Np nx ny area

radius = sqrt(43560*area/pi);

thres = 0.5*radius^2*ones(Np,1);
thres(randperm(Np,0.3*Np),1) = 0.001;


pos = []; vio = [];

for i = 1:Np
    
    tmp3 = (1)*radius^2;
    
    while tmp3 > thres(i,1)
        
        tmp = randperm(nx*ny, N);
        tmp = reshape([ceil(tmp'/ny); mod(tmp', ny) + (mod(tmp',ny) == 0)*ny], 1, 2*N);

        if isempty(type)

            tmp2 = [1/3 + 2/3*rand(1,1), -1 + 2*rand(1,N-1)];
            tmp3 = constViolation(tmp, tmp2);
            
        else
            
            tmp2 = [];
            tmp3 = constViolation(tmp, type);
            
        end
                
    end
        
    pos = [pos; [tmp, tmp2]];
    vio = [vio; tmp3];
    
end

end


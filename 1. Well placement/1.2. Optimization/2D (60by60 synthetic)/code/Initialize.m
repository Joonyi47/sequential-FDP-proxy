function [pos] = Initialize(type)
global N Np nx ny area

radius = sqrt(43560*area/pi);

thres = 0.5*radius^2*ones(Np,1);
thres(randperm(Np,0.3*Np),1) = 0.001;


pos = []; vio = [];

for i = 1:Np
    
    violation = (1)*radius^2;
    
    while violation > thres(i,1)
        
        loc_idx = randperm(nx*ny, N);
        loc = reshape([ceil(loc_idx'/ny); mod(loc_idx', ny) + (mod(loc_idx',ny) == 0)*ny], 1, 2*N);

        if isempty(type)

            tp = [1/3 + 2/3*rand(1,1), -1 + 2*rand(1,N-1)];
            violation = constViolation(loc, tp);
            
        else
            
            tp = [];
            violation = constViolation(loc, type);
            
        end
                
    end
        
    pos = [pos; [loc, tp]];
    vio = [vio; violation];
    
end

end


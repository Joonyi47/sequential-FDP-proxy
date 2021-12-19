function [type] = Initialize(p)
global N Nstep Np

type = [];
vio  = [];

thres = 3*ones(Np,1);
thres(randperm(Np,0.5*Np),1) = 0.001;

h = waitbar(0,'Please wait...');
for i = 1:Np
    
    tmp2 = max(thres)+1;
    
    while tmp2 > thres(i,1)
        
        tmp = [rand(1,N)*ceil(N/3)/Nstep; rand(1,N); pi*rand(2,N)];
        tmp = reshape(tmp, 1, 4*N);
        
        tmp2 = constViolation(tmp);
        
    end
    
    type = [type; tmp];
    vio  = [vio; tmp2];
    
    waitbar(i / Np)
    
end

close(h)
end
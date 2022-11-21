function [operation] = Initialize(p)
global N Nstep Np

operation = [];
vio  = [];

thres = 3*ones(Np,1);
thres(randperm(Np,1*Np),1) = 0.001;

h = waitbar(0,'Please wait...');
for i = 1:Np
    
    violdation = max(thres)+1;
    
    while violdation > thres(i,1)
        
        opr = [rand(1,N)*ceil(N/3)/Nstep; rand(1,N); pi*rand(2,N)];
        opr = reshape(opr, 1, 4*N);
        
        violdation = constViolation(opr);
        
    end
    
    operation = [operation; opr];
    vio  = [vio; violdation];
    
    waitbar(i / Np)
    
end

close(h)
end


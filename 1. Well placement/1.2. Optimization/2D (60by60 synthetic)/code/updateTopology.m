function REP = updateTopology(REP) 
    global Np Ni
    % topology update if best fitness value not improved
    if REP.pos_fit >= REP.prior.pos_fit
        
        mat = zeros(Np, Np);
        p = 1-(1-1/Np)^Ni;
        for i = 1:Np
            for j = 1:Np
                if i == j
                    mat(i,j) = 1;
                elseif (rand < p)
                    mat(i,j) = 1;
                end
            end
        end
        REP.neighbor_mat = mat;
    end
end
function REP = updateRepository(REP,POS,POS_fit,POS_vio)
    REP.prior  = REP;
    % Domination between particles
    DOMINATED1 = checkDomination(POS_fit(:,1),POS_vio);
    SELECTED1  = find(~DOMINATED1 == 1, 1);
    TEMP       = [REP.pos;POS(SELECTED1,:)];
    TEMP_fit   = [REP.pos_fit;POS_fit(SELECTED1,:)];
    TEMP_vio   = [REP.pos_vio;POS_vio(SELECTED1,:)];
    DOMINATED2 = checkDomination(TEMP_fit(:,1), TEMP_vio);
    if ~all(DOMINATED2==0)
        REP.pos = TEMP(~DOMINATED2,:);
        REP.pos_fit = TEMP_fit(~DOMINATED2,:);
        REP.pos_vio = TEMP_vio(~DOMINATED2,:);
    end
    % Updating the grid
    REP        = updateTopology(REP);
end
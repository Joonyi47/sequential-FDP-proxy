
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                      2D 2P IMPES                         %%%%%%%
%%%%%%              2019. 07.  -   by Joonyi Kim                %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%} time dependent 행렬만 첫글자 대문자 %%%%%%


function [P, So] = IMPES(pos, type, perm, swof, parameters)

names = fieldnames(parameters);
for k = 1:length(names)
    eval([names{k} '=parameters.' names{k} ';' ]);
end

% global nx ny nz dx dy dz dt T nstep...
%     pi phi cf ps co ro soi boi voi ...
%     cw vwi bwi swi ...
%     s rw pwf_p pwf_i

% load sample.mat;
% load PERMX.mat;
% load RelPermTable.mat;

loc_p = reshape(pos(:,repelem(type,1,2) == 1), 2, length(find(type == 1)));
loc_i = reshape(pos(:,repelem(type,1,2) == 0), 2, length(find(type == 0)));

% nx = 20;
% ny = nx;
% nz = 1;
% 
% dx = 240;
% dy = 240;
% dz = 40;

perm = reshape(reshape(perm, nx, ny)', nx*ny, 1);
% perm = reshape(reshape(PERMX.ups1(:,1), nx, ny)', nx*ny, 1);
% phi  = 0.2*ones(length(perm),1);
% cf   = 3E-05;
% pi   = 3500; % psi
% ps   = 14.7;

% dt = 15;
% T  = 90;
% nstep = T/dt;

% %%% oil properties %%%
% co  = 9.98E-07;
% ro  = 0.2*dx;
% % ro  = sqrt(dx*dy/3.141592);
% soi = 0.75;
% boi = 1.01/(1 - co * (2000 - pi));
% voi = 3.0;
% 
% %%% water properties %%%
% cw  = 5E-07;
% vwi = 1;
% swi = 1-soi;
% bwi = 1;            % constant



%%%%% Initialize  %%%%%
cell_prod = ny.*(loc_p(1,:) - 1) + loc_p(2,:);
cell_inj  = ny.*(loc_i(1,:) - 1) + loc_i(2,:);

P   =  prei*ones(nx*ny, nstep);
So  =  soi*ones(nx*ny, nstep);
Sw  =  swi*ones(nx*ny, nstep);
% Sw(cell_inj,:) = 1-swi;
Ct  =  cf + co*So + cw*Sw;
Vo  =  voi*ones(nx*ny, nstep);
Vw  =  vwi*ones(nx*ny, nstep);
Bo  =  boi*ones(nx*ny, nstep);
Bw  =  bwi*ones(nx*ny, nstep);
Vp  =  dx*dy*dz*phi.*ones(nx*ny, nstep);

A   = zeros(nx*ny, nx*ny, nstep);
B   = zeros(nx*ny, 1, nstep);


% pwf_p = 1000;
% pwf_i = 6000;
% rw    = 0.5;  % ft
% s     = 0;    % skin factor



J_model = zeros(nx*ny,1);
% J_model(cell_prod) = 2*3.141592*0.00633*perm(cell_prod)*dz/(log(ro/rw) + s);
J_model(cell_prod) = 2*3.141592*1.959178;
% J_model(cell_inj) = 2*3.141592*0.00633*perm(cell_inj)*dz/(log(ro/rw) + s);
J_model(cell_inj) = 2*3.141592*1.959178;


%% main loop
for t = 1:nstep
    
    temp  = zeros(nx*ny,nx*ny + 2*ny-1);  % for A matrix
    temp2 = zeros(nx*ny, 1);              % for B matrix
    relperm_oil = CalculateRelPerm(swof, So(:, t), 'o');
    relperm_wat = CalculateRelPerm(swof, Sw(:, t), 'w');
    relperm_wat(cell_inj) = max(swof(:,2));

    %     lamda = (relperm_oil./Vo(:, t) + relperm_wat./Vw(:, t));
    
    %%%%%%%% A matrix %%%%%%%%%%%
    for cellnum = 1:nx*ny
        
        
        if ( cellnum - ny ) > 0   % West 여부
            I = P(cellnum,t) < P(cellnum - ny,t);
            kro = relperm_oil(cellnum - ny*I);
            krw = relperm_wat(cellnum - ny*I);
            perm_ = 2*perm(cellnum)*perm(cellnum-ny)/(perm(cellnum) + perm(cellnum-ny));
            bo  = ( Bo(cellnum - ny, t) + Bo(cellnum, t) ) / 2;
            bw  = ( Bw(cellnum - ny, t) + Bw(cellnum, t) ) / 2;
            vo  = ( Vo(cellnum - ny, t) + Vo(cellnum, t) ) / 2;
            vw  = ( Vw(cellnum - ny, t) + Vw(cellnum, t) ) / 2;
            aw   = Bo(cellnum, t) * 0.00633 * dy * dz / dx * perm_ * kro / bo / vo ...
                + Bw(cellnum, t) * 0.00633 * dy * dz / dx * perm_ * krw / bw / vw;
        else
            aw = 0;
        end
        temp(cellnum, cellnum) = -aw;
        
        if mod(cellnum - 1, ny) ~= 0  % North 여부
            I = P(cellnum,t) < P(cellnum - 1,t);
            kro = relperm_oil(cellnum-1*I);
            krw = relperm_wat(cellnum-1*I);
            perm_ = 2*perm(cellnum)*perm(cellnum-1)/(perm(cellnum) + perm(cellnum-1));
            bo   = ( Bo(cellnum - 1, t) + Bo(cellnum, t) ) / 2;
            bw   = ( Bw(cellnum - 1, t) + Bw(cellnum, t) ) / 2;
            vo   = ( Vo(cellnum - 1, t) + Vo(cellnum, t) ) / 2;
            vw   = ( Vw(cellnum - 1, t) + Vw(cellnum, t) ) / 2;
            an    = Bo(cellnum, t) * 0.00633 * dx * dz / dy * perm_ * kro / bo / vo ...
                + Bw(cellnum, t) * 0.00633 * dx * dz / dy * perm_ * krw / bw / vw;
        else
            an = 0;
        end
        temp(cellnum, cellnum + ny - 1) = -an;
        
        if mod(cellnum, ny) ~= 0  % South 여부
            I = P(cellnum,t) < P(cellnum + 1,t);
            kro = relperm_oil(cellnum+1*I);
            krw = relperm_wat(cellnum+1*I);
            perm_ = 2*perm(cellnum)*perm(cellnum+1)/(perm(cellnum) + perm(cellnum+1));
            bo   = ( Bo(cellnum + 1, t) + Bo(cellnum, t) ) / 2;
            bw   = ( Bw(cellnum + 1, t) + Bw(cellnum, t) ) / 2;
            vo   = ( Vo(cellnum + 1, t) + Vo(cellnum, t) ) / 2;
            vw   = ( Vw(cellnum + 1, t) + Vw(cellnum, t) ) / 2;
            as    = Bo(cellnum, t) * 0.00633 * dx * dz / dy * perm_ * kro / bo / vo ...
                + Bw(cellnum, t) * 0.00633 * dx * dz / dy * perm_ * krw / bw / vw;
        else
            as = 0;
        end
        temp(cellnum, cellnum + ny + 1) =  -as;
        
        if ( cellnum + ny ) < nx*ny + 1  % East 여부
            I = P(cellnum,t) < P(cellnum + ny,t);
            kro = relperm_oil(cellnum + ny*I);
            krw = relperm_wat(cellnum + ny*I);
            perm_ = 2*perm(cellnum)*perm(cellnum+ny)/(perm(cellnum) + perm(cellnum+ny));
            bo   = ( Bo(cellnum + ny, t) + Bo(cellnum, t) ) / 2;
            bw   = ( Bw(cellnum + ny, t) + Bw(cellnum, t) ) / 2;
            vo   = ( Vo(cellnum + ny, t) + Vo(cellnum, t) ) / 2;
            vw   = ( Vw(cellnum + ny, t) + Vw(cellnum, t) ) / 2;
            ae    = Bo(cellnum, t) * 0.00633 * dy * dz / dx * perm_ * kro / bo / vo ...
                + Bw(cellnum, t) * 0.00633 * dy * dz / dx * perm_ * krw / bw / vw;
        else
            ae = 0;
        end
        temp(cellnum, cellnum + ny + ny) = -ae;
        
        av = Vp(cellnum, t)*Ct(cellnum, t)/dt;
        
        if ismember(cellnum, cell_prod)
            lamda = (relperm_oil(cellnum)./Vo(cellnum, t) + relperm_wat(cellnum)./Vw(cellnum, t));
        elseif ismember(cellnum, cell_inj)
            lamda = (relperm_wat(cellnum)./Vw(cellnum, t));
%             lamda = (Bo(cellnum, t)/Bw(cellnum, t)*relperm_oil(cellnum)./Vo(cellnum, t) + relperm_wat(cellnum)./Vw(cellnum, t));
        else
            lamda = 0;
        end
        ac = (aw+as+ae+an+av + J_model(cellnum)*lamda);
        temp(cellnum, cellnum + ny) = ac;
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%% B matrix %%%%%%%%%
        if ismember(cellnum, cell_prod)
            temp2(cellnum, 1) = av*P(cellnum, t) + J_model(cellnum)*lamda*pwf_p;
        elseif ismember(cellnum, cell_inj)
            temp2(cellnum, 1) = av*P(cellnum, t) + J_model(cellnum)*lamda*pwf_i;
        else
            temp2(cellnum, 1) = av*P(cellnum, t);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    A = temp(:, ny+1:ny+ny*nx );
    B = temp2;
    
    P(:,t+1) = A\B;
    %     P(:,t+1) = TDMA(A,B,size(A,1));
    
    %%%%%%%%% parameter update %%%%%%%%%%
    Vp(:, t+1) = Vp(:,t) .* (1+cf*(P(:,t+1) - P(:,t)));
    Bo(:, t+1) = Bo(:,t) .* (1-co*(P(:,t+1) - P(:,t)));
    Bw(:, t+1) = Bw(:,t) .* (1-cw*(P(:,t+1) - P(:,t)));
    
    P_temp = [zeros(ny, 1); P(:,t+1); zeros(ny,1)];
    for cellnum = 1:nx*ny
        
        if ismember(cellnum, cell_inj)
            pwf = pwf_i; D = -1;
        elseif ismember(cellnum, cell_prod)
            pwf = pwf_p; D = 1;
        else
            pwf = P_temp(cellnum + ny); D = 0;
        end
        
        if ( cellnum - ny ) > 0   % West 여부
            I = P(cellnum,t) < P(cellnum - ny,t);
            kro = relperm_oil(cellnum - ny*I);
            perm_ = 2*perm(cellnum)*perm(cellnum-ny)/(perm(cellnum) + perm(cellnum-ny));
            bo  = ( Bo(cellnum - ny, t) + Bo(cellnum, t) ) / 2;
            vo  = ( Vo(cellnum - ny, t) + Vo(cellnum, t) ) / 2;
            aw  = Bo(cellnum, t) * 0.00633 * dy * dz / dx * perm_ * kro / bo / vo;
        else
            aw = 0;
        end
        
        if mod(cellnum - 1, ny) ~= 0  % North 여부
            I = P(cellnum,t) < P(cellnum - 1,t);
            kro = relperm_oil(cellnum-1*I);
            perm_ = 2*perm(cellnum)*perm(cellnum-1)/(perm(cellnum) + perm(cellnum-1));
            bo   = ( Bo(cellnum - 1, t) + Bo(cellnum, t) ) / 2;
            vo   = ( Vo(cellnum - 1, t) + Vo(cellnum, t) ) / 2;
            an   = Bo(cellnum, t) * 0.00633 * dx * dz / dy * perm_ * kro / bo / vo;
        else
            an = 0;
        end
        
        if mod(cellnum, ny) ~= 0  % South 여부
            I = P(cellnum,t) < P(cellnum + 1,t);
            kro = relperm_oil(cellnum+1*I);
            perm_ = 2*perm(cellnum)*perm(cellnum+1)/(perm(cellnum) + perm(cellnum+1));
            bo   = ( Bo(cellnum + 1, t) + Bo(cellnum, t) ) / 2;
            vo   = ( Vo(cellnum + 1, t) + Vo(cellnum, t) ) / 2;
            as   = Bo(cellnum, t) * 0.00633 * dx * dz / dy * perm_ * kro / bo / vo;
        else
            as = 0;
        end
        
        if ( cellnum + ny ) < nx*ny + 1  % East 여부
            I = P(cellnum,t) < P(cellnum + ny,t);
            kro = relperm_oil(cellnum + ny*I);
            perm_ = 2*perm(cellnum)*perm(cellnum+ny)/(perm(cellnum) + perm(cellnum+ny));
            bo   = ( Bo(cellnum + ny, t) + Bo(cellnum, t) ) / 2;
            vo   = ( Vo(cellnum + ny, t) + Vo(cellnum, t) ) / 2;
            ae   = Bo(cellnum, t) * 0.00633 * dy * dz / dx * perm_ * kro / bo / vo;
        else
            ae = 0;
        end
%         
        So(cellnum, t+1) = Bo(cellnum, t+1)/Vp(cellnum, t+1)*...
            ((aw * (P_temp(cellnum) - P_temp(cellnum + ny)) ...
            + an * (P_temp(cellnum + ny - 1) - P_temp(cellnum + ny)) ...
            + as * (P_temp(cellnum + ny + 1) - P_temp(cellnum + ny)) ...
            + ae * (P_temp(cellnum + ny + ny) - P_temp(cellnum + ny)) ...
            - D*J_model(cellnum) * (relperm_oil(cellnum) / Vo(cellnum, t) / Bo(cellnum, t)) * (P_temp(cellnum + ny) - pwf))*dt ...
            + Vp(cellnum, t)*So(cellnum, t)/Bo(cellnum, t));
% if cellnum == 19
% keyboard()
% end
%         Sw(cellnum, t+1) = Bw(cellnum, t+1)/Vp(cellnum, t+1)*...
%             ((aw * (P_temp(cellnum + ny - ny) - P_temp(cellnum + ny)) ...
%             + an * (P_temp(cellnum + ny - 1) - P_temp(cellnum + ny)) ...
%             + as * (P_temp(cellnum + ny + 1) - P_temp(cellnum + ny)) ...
%             + ae * (P_temp(cellnum + ny + ny) - P_temp(cellnum + ny)) ...
%             + J_model(cellnum) * (relperm_wat(cellnum)./Vw(cellnum, t) / Bw(cellnum, t)) * (P_temp(cellnum + ny) - pwf))*dt ...
%             + Vp(cellnum, t)*Sw(cellnum, t)/Bw(cellnum, t));
        
%                 So(cellnum, t+1) = Bo(cellnum, t+1)/Vp(cellnum, t+1)*...
%                                    ((aw * (cellnum, cellnum) * (P_temp(cellnum) - P_temp(cellnum + ny)) ...
%                                    + an * (cellnum, cellnum + ny - 1) * (P_temp(cellnum + ny - 1) - P_temp(cellnum + ny)) ...
%                                    + as * (cellnum, cellnum + ny + 1) * (P_temp(cellnum + ny + 1) - P_temp(cellnum + ny)) ...
%                                    + ae * (cellnum, cellnum + ny + ny) * (P_temp(cellnum + ny + ny) - P_temp(cellnum + ny)) ...
%                                    - D*J_model(cellnum) * (relperm_oil(cellnum) / Vo(cellnum, t) / Bo(cellnum, t)) * (P_temp(cellnum + ny) - pwf))*dt ...
%                                    + Vp(cellnum, t)*So(cellnum, t)/Bo(cellnum, t));
%         if abs(So(cellnum, t+1) - So(cellnum, t)) > 0.2
%             So(cellnum, t+1) = So(cellnum, t) + 0.1*(So(cellnum, t+1) - So(cellnum, t));
%         elseif 0.2 > abs(So(cellnum, t+1) - So(cellnum, t)) && abs(So(cellnum, t+1) - So(cellnum, t)) > 0.1
%             So(cellnum, t+1) = So(cellnum, t) + 0.2*(So(cellnum, t+1) - So(cellnum, t));
%         end
        if So(cellnum, t+1) < 0.25
            So(cellnum, t+1) = 0.25;
        elseif So(cellnum, t+1) > 0.75
            So(cellnum, t+1) = 0.75;
        end
        Sw(cellnum, t+1) = 1-So(cellnum, t+1);
        
    end
    
    Ct(:, t+1) = cf+co.*So(:,t+1)+cw.*Sw(:,t+1);
end

end


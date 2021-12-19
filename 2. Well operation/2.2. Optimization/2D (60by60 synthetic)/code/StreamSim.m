
function [TOF_beg, TOF_end] = StreamSim(pos, type, perm, swof, P, So, parameters)
names = fieldnames(parameters);
for k = 1:length(names)
    eval([names{k} '=parameters.' names{k} ';' ]);
end
% global nx ny nz dx dy dz dt T...
%     pi phi cf ps co ro soi boi voi ...
%     cw vwi bwi swi ...
%     s rw pwf_p pwf_i

perm = reshape(reshape(perm, nx, ny)', nx*ny, 1);

loc_p = reshape(pos(:,repelem(type,1,2) == 1), 2, length(find(type == 1)));
loc_i = reshape(pos(:,repelem(type,1,2) == 0), 2, length(find(type == 0)));

cell_prod = ny.*(loc_p(1,:) - 1) + loc_p(2,:);
cell_inj  = ny.*(loc_i(1,:) - 1) + loc_i(2,:);

nline = 10; max_tof = 10000;

nstart = length(find(type == 0));
X_w = cell(nstart,nline); Y_w = cell(nstart,nline); 
X_n = cell(nstart,nline); Y_n = cell(nstart,nline); 
X_s = cell(nstart,nline); Y_s = cell(nstart,nline); 
X_e = cell(nstart,nline); Y_e = cell(nstart,nline); 
Cell_w = cell(nstart,nline); Delt_w = cell(nstart,nline);
Cell_n = cell(nstart,nline); Delt_n = cell(nstart,nline);
Cell_s = cell(nstart,nline); Delt_s = cell(nstart,nline);
Cell_e = cell(nstart,nline); Delt_e = cell(nstart,nline);

relperm = [CalculateRelPerm(swof, So, 'o'), CalculateRelPerm(swof, 1-So, 'w')];
    
% figure();
% hold on;  grid on; axis ij
% xlim([0 nx*dx]); ylim([0 ny*dy]);

npass = zeros(nx*ny,1);
for k = 1:nstart
    startcell = cell_inj(k);
    
    if ( startcell - ny ) > 0 % West 여부
        
        coord0 = [ceil(startcell/ny)- 1, mod(startcell,ny)+ny*(mod(startcell,ny)==0)] .* [dx, dy]; % startcell의 왼쪽아래모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell - ny; % 들어가는 cell
            x = coord0(1); y = coord0(2) - (dy/nline*(numlin-1)+0.05);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];
                
                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'origin');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_prod) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_w{k, numlin} = temp1;
            Y_w{k, numlin} = temp2;
            Cell_w{k, numlin} = temp3;
            Delt_w{k, numlin} = temp4;
            
            npass(Cell_w{k, numlin}) = npass(Cell_w{k, numlin}) + 1;
%             plot(X_w{k, numlin}, Y_w{k, numlin}, 'b');
            
        end
    else
        
        
    end
    
    if mod(startcell - 1, ny) ~= 0 % North 여부
        
        coord0 = [ceil(startcell/ny)- 1, mod(startcell,ny)+ny*(mod(startcell,ny)==0)-1] .* [dx, dy]; % startcell의 왼쪽위모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell - 1; % 들어가는 cell
            x = coord0(1) + (dx/nline*(numlin-1)+0.05); y = coord0(2);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];

                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'origin');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_prod) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_n{k, numlin} = temp1;
            Y_n{k, numlin} = temp2;
            Cell_n{k, numlin} = temp3;
            Delt_n{k, numlin} = temp4;
            
            npass(Cell_w{k, numlin}) = npass(Cell_w{k, numlin}) + 1;
%             plot(X_n{k, numlin}, Y_n{k, numlin}, 'b');
            
        end
    else
        
        
    end
    
    if mod(startcell, ny) ~= 0  % South 여부
        
        coord0 = [ceil(startcell/ny), mod(startcell,ny)+ny*(mod(startcell,ny)==0)] .* [dx, dy]; % startcell의 왼쪽아래모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell + 1; % 들어가는 cell
            x = coord0(1) - (dx/nline*(numlin-1)+0.05); y = coord0(2);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];
%                 if numlin == 4
%                     keyboard()
%                 end
                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'origin');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_prod) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_s{k, numlin} = temp1;
            Y_s{k, numlin} = temp2;
            Cell_s{k, numlin} = temp3;
            Delt_s{k, numlin} = temp4;
            
            npass(Cell_w{k, numlin}) = npass(Cell_w{k, numlin}) + 1;
            
%             plot(X_s{k, numlin}, Y_s{k, numlin}, 'b');
            
        end
    else
        
        
    end
    
    if ( startcell + ny ) < nx*ny + 1  % East 여부
        
        coord0 = [ceil(startcell/ny), mod(startcell,ny)+ny*(mod(startcell,ny)==0)-1] .* [dx, dy]; % startcell의 왼쪽아래모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell + ny; % 들어가는 cell
            x = coord0(1); y = coord0(2) + (dy/nline*(numlin-1)+0.05);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];
                
                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'origin');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_prod) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_e{k, numlin} = temp1;
            Y_e{k, numlin} = temp2;
            Cell_e{k, numlin} = temp3;
            Delt_e{k, numlin} = temp4;
            
            npass(Cell_w{k, numlin}) = npass(Cell_w{k, numlin}) + 1;
            
%             plot(X_e{k, numlin}, Y_e{k, numlin}, 'b');            
            
        end
    else
        
        
    end

    
end
% hold off

%%%% TOF %%%%
TOF_beg = zeros(nx*ny,1);
Cell = [Cell_w;Cell_n;Cell_e;Cell_s];
Delt = [Delt_w;Delt_n;Delt_e;Delt_s];

tof_beg = zeros(nx*ny,1);
for k = 1:nstart*4
   for numlin = 1:nline
       tof_beg_temp = zeros(nx*ny,1);
       tof_beg_temp(Cell{k,numlin}) = cumsum(Delt{k,numlin});
       tof_beg = [tof_beg, tof_beg_temp]; 
   end
end
TOF_beg = sum(tof_beg(:,2:end), 2)./sum((tof_beg > 0), 2);
TOF_beg(isnan(TOF_beg)) = 0;
TOF_beg(TOF_beg == Inf) = max_tof; 
TOF_beg(TOF_beg == 0) = max_tof;
TOF_beg(TOF_beg > max_tof) = max_tof;
TOF_beg(cell_inj) = 0;
    

%%%%%%%%%%%%%%%%% source = prod %%%%%%%%%%%%%%%%%%%%%%
nstart = length(find(type == 1));
X_w = cell(nstart,nline); Y_w = cell(nstart,nline); 
X_n = cell(nstart,nline); Y_n = cell(nstart,nline); 
X_s = cell(nstart,nline); Y_s = cell(nstart,nline); 
X_e = cell(nstart,nline); Y_e = cell(nstart,nline); 
Cell_w = cell(nstart,nline); Delt_w = cell(nstart,nline);
Cell_n = cell(nstart,nline); Delt_n = cell(nstart,nline);
Cell_s = cell(nstart,nline); Delt_s = cell(nstart,nline);
Cell_e = cell(nstart,nline); Delt_e = cell(nstart,nline);

% figure();
% hold on;  grid on; axis ij
% xlim([0 nx*dx]); ylim([0 ny*dy]);

for k = 1:nstart
    startcell = cell_prod(k);
    
    if ( startcell - ny ) > 0 % West 여부
        
        coord0 = [ceil(startcell/ny)- 1, mod(startcell,ny)+ny*(mod(startcell,ny)==0)] .* [dx, dy]; % startcell의 왼쪽아래모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell - ny; % 들어가는 cell
            x = coord0(1); y = coord0(2) - (dy/nline*(numlin-1)+0.05);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];

                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'reverse');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_inj) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_w{k, numlin} = temp1;
            Y_w{k, numlin} = temp2;
            Cell_w{k, numlin} = temp3;
            Delt_w{k, numlin} = temp4;
            
%             plot(X_w{k, numlin}, Y_w{k, numlin}, 'b');
            
        end
    else
        
        
    end
    
    if mod(startcell - 1, ny) ~= 0 % North 여부
        
        coord0 = [ceil(startcell/ny)- 1, mod(startcell,ny)+ny*(mod(startcell,ny)==0)-1] .* [dx, dy]; % startcell의 왼쪽위모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell - 1; % 들어가는 cell
            x = coord0(1) + (dx/nline*(numlin-1)+0.05); y = coord0(2);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];
%                 if k == 1 && numlin == 10
%                     keyboard()
%                 end
                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'reverse');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_inj) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_n{k, numlin} = temp1;
            Y_n{k, numlin} = temp2;
            Cell_n{k, numlin} = temp3;
            Delt_n{k, numlin} = temp4;
            
%             plot(X_n{k, numlin}, Y_n{k, numlin}, 'b');
            
        end
    else
        
        
    end
    
    if mod(startcell, ny) ~= 0  % South 여부
        
        coord0 = [ceil(startcell/ny), mod(startcell,ny)+ny*(mod(startcell,ny)==0)] .* [dx, dy]; % startcell의 왼쪽아래모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell + 1; % 들어가는 cell
            x = coord0(1) - (dx/nline*(numlin-1)+0.05); y = coord0(2);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];
                
                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'reverse');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_inj) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_s{k, numlin} = temp1;
            Y_s{k, numlin} = temp2;
            Cell_s{k, numlin} = temp3;
            Delt_s{k, numlin} = temp4;
            
%             plot(X_s{k, numlin}, Y_s{k, numlin}, 'b');
            
        end
    else
        
        
    end
    
    if ( startcell + ny ) < nx*ny + 1  % East 여부
        
        coord0 = [ceil(startcell/ny), mod(startcell,ny)+ny*(mod(startcell,ny)==0)-1] .* [dx, dy]; % startcell의 왼쪽아래모서리 좌표  
        
        for numlin = 1:nline
            cellnum = startcell + ny; % 들어가는 cell
            x = coord0(1); y = coord0(2) + (dy/nline*(numlin-1)+0.05);
            temp1 = x; temp2 = y; temp3 = startcell; temp4 = 0;

            stopcondition = false;
            while ~stopcondition
                
                temp3 = [temp3, cellnum];
                
                [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, 'reverse');
                
                temp1 = [temp1, x]; temp2 = [temp2, y]; temp4 = [temp4, del_t];
                
                stopcondition = ismember(cellnum, cell_inj) || cellnum < 1 || cellnum > nx*ny || ismember(cellnum, temp3) || del_t > max_tof;
            end
            X_e{k, numlin} = temp1;
            Y_e{k, numlin} = temp2;
            Cell_e{k, numlin} = temp3;
            Delt_e{k, numlin} = temp4;
            
%             plot(X_e{k, numlin}, Y_e{k, numlin}, 'b');            
            
        end
    else
        
        
    end

    
end
% hold off

%%%% TOF %%%%
TOF_end = zeros(nx*ny,1);
Cell = [Cell_w;Cell_n;Cell_e;Cell_s];
Delt = [Delt_w;Delt_n;Delt_e;Delt_s];

tof_end = zeros(nx*ny,1);
for k = 1:nstart*4
   for numlin = 1:nline
       tof_end_temp = zeros(nx*ny,1);
       tof_end_temp(Cell{k,numlin}) = cumsum(Delt{k,numlin});
       tof_end = [tof_end, tof_end_temp]; 
   end
end
TOF_end = sum(tof_end(:,2:end), 2)./sum((tof_end > 0), 2);
TOF_end(isnan(TOF_end)) = 0;
TOF_end(TOF_end == Inf) = max_tof; 
TOF_end(TOF_end == 0) = max_tof;
TOF_end(TOF_end > max_tof) = max_tof;
TOF_end(cell_prod) = 0;

% figure; imagesc(reshape(TOF_beg, nx, ny)); colormap jet
% figure; imagesc(reshape(TOF_end, nx, ny)); colormap jet

end

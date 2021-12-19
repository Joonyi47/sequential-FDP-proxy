function [viot] = constViolation(pos, type)
  %%% pos ´Â ¿­º¤ÅÍ = {x1, y1, x2, y2, ... }
  
  global N dx dy area nx ny ACTIVE
  
  Nofconstraint = N * ( N - 1) / 2;
  loc = repmat([dx, dy],1,N) .* pos;
  radius = sqrt(43560*area/pi);
  active = reshape(ACTIVE, nx*ny, 2);
  active = active(:,1).*active(:,2); 
  active = reshape(active, nx, ny)';
  
  type(type >= 1/3) = 1; % prod.
  type(type <= -1/3) = -1; % injec. 
  type(- 1/3 < type & type < 1/3) = 0; % no well
  
  drWell.Index = repelem(type,1,2);
  drWell.loc   = loc(drWell.Index ~= 0);
  drWell.pos   = reshape(pos(drWell.Index ~= 0), 2, sum(type ~= 0));
  drWell.Dloc  = reshape(drWell.loc,2,sum(type ~= 0))';
  
  pWell.Index  = repelem(type,1,2);
  pWell.pos    = pos(pWell.Index == 1);
  pWell.loc    = loc(pWell.Index == 1);
  pWell.Dloc   = reshape(pWell.loc,2,sum(type == 1))';   % dimension form  {x1, y1; x2, y2; ,,, }
  pWell.Bloc   = [1,1,dx,dy] .* [3, 3, nx-5, ny-5];            % boundary
   
  if sum(type~=0) == 1
      all_perm = [];
  else
      all_perm = nchoosek(1:sum(type~=0), 2);
  end
  V = [];
  for i = 1:size(all_perm,1)
      V = [V; radius^2 - (pdist([drWell.Dloc(all_perm(i,1),:); ...
          drWell.Dloc(all_perm(i,2),:)]))^2];
  end
  for i = 1:sum(type~=0)   % check well on active cell
      V = [V; radius^2 * ~active(drWell.pos(2,i), drWell.pos(1,i))];
  end
      V = [V; radius^2 - Dbound(pWell.pos, active)];
%   for i = 1:sum(type~=0)
%       V = [V; radius^2*(active(pos(1,2*i), pos(1,2*(i-1)+1)) == 0)];
%   end
  V(V<0) = 0;        % satisfies 40acres rule
  viot = sum(V);

end
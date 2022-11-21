function [viot] = constViolation(pos, type)
  %%% pos ´Â ¿­º¤ÅÍ = {x1, y1, x2, y2, ... }
  
  global N dx area nx
  
  Nofconstraint = N * ( N - 1) / 2;
  loc = dx * pos;
  radius = sqrt(43560*area/pi);
  
  type(type >= 1/3) = 1; % prod.
  type(type <= -1/3) = -1; % injec. 
  type(- 1/3 < type & type < 1/3) = 0; % no well
  
  drWell.Index = repelem(type,1,2);
  drWell.loc   = loc(drWell.Index ~= 0);
  drWell.Dloc  = reshape(drWell.loc,2,sum(type ~= 0))';
  
  pWell.Index  = repelem(type,1,2);
  pWell.loc    = loc(pWell.Index == 1);
  pWell.Dloc   = reshape(pWell.loc,2,sum(type == 1))';   % dimension form  {x1, y1; x2, y2; ,,, }
  pWell.Bloc   = dx * [0, nx+1];            % boundary
  
  if sum(type~=0) > 1
      all_perm = nchoosek(1:sum(type~=0), 2);
      V = [];
      for i = 1:size(all_perm,1)
          V = [V; radius^2 - (pdist([drWell.Dloc(all_perm(i,1),:); ...
              drWell.Dloc(all_perm(i,2),:)]))^2];
      end
  else
      V = [];
  end

  for i = 1:sum(type == 1)
      V = [V; radius^2 - min(abs(diff([repelem(pWell.Bloc,1,2);repmat(pWell.Dloc(i,:),1,2)])))^2];
  end
  V(V<0) = 0;        % satisfies 40acres rule
  viot = sum(V);

end
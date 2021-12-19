function [viot] = activeViolation(pos, type)
  %%% pos ´Â ¿­º¤ÅÍ = {x1, y1, x2, y2, ... }
  
  global N area nx ny ACTIVE
  
  radius = sqrt(43560*area/pi);
  active = reshape(ACTIVE, nx*ny, 2);
  active = active(:,1).*active(:,2); 
  active = reshape(active, nx, ny)';
  
  drWell.pos   = reshape(pos, 2, N);
  
  V = [];

  for i = 1:size(type,2)
     V = [V; radius^2 * ~active(drWell.pos(2,i), drWell.pos(1,i))];
  end
 
  V(V<0) = 0;        % satisfies 40 acre's rule
  viot = max(V);

end
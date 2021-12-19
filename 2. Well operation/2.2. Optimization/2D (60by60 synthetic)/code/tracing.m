% cellnum 은 line이 들어가는 cell의 number, 출력은 다음으로 들어가는 cell number
% function [x, y, cellnum, del_t] = tracing(x, y, cellnum, P, So, perm, swof, direction)
function [x, y, cellnum, del_t] = tracing(parameters, x, y, cellnum, P, perm, relperm, direction)
  
  nx = parameters.nx;
  ny = parameters.ny;
  dx = parameters.dx;
  dy = parameters.dy;
  phi = parameters.phi;
  voi = parameters.voi;
  vwi = parameters.vwi;
  
     
  if strcmp(direction, 'reverse')
      D = -1; O = 1; W = 1;
  else
      D = 1; O = 1; W = 1;
  end
  
  P_ = nan*ones(ny+2, nx+2);
  P_(2:end-1, 2:end-1) = reshape(P,ny,nx);
  relperm_oil = relperm(:,1);
  relperm_wat = relperm(:,2);
  cellnum_ = (ny+2) * (ceil(cellnum/ny)) + 1 + mod(cellnum,ny)+ny*(mod(cellnum,ny)==0);
  
  base = [ceil(cellnum/ny)- 1, mod(cellnum,ny)+ny*(mod(cellnum,ny)==0)] .* [dx, dy];
  x0 = base(1); y0 = base(2);
  
  %} West velocity
  try
    isnan(P_(cellnum_ - (ny+2)));
  catch me
      keyboard()
  end
  if ~isnan(P_(cellnum_ - (ny+2)))
      I = P(cellnum) < P(cellnum - ny);
      kro = relperm_oil(cellnum - ny*I);
      krw = relperm_wat(cellnum - ny*I);
      vw = 0.00633/phi(cellnum)*2*perm(cellnum)*perm(cellnum-ny)/(perm(cellnum)+perm(cellnum-ny)) ...
          *(P_(cellnum_ - (ny+2)) - P_(cellnum_))/dx*(O*kro/voi+W*krw/vwi);
      vw = D*vw;
  else
      vw = 0;
  end
  
  %} East velocity
  if ~isnan(P_(cellnum_ + (ny+2)))
      I = P(cellnum) < P(cellnum + ny);
      kro = relperm_oil(cellnum + ny*I);
      krw = relperm_wat(cellnum + ny*I);
      ve = 0.00633/phi(cellnum)*2*perm(cellnum)*perm(cellnum+ny)/(perm(cellnum)+perm(cellnum+ny)) ...
          *(P_(cellnum_) - P_(cellnum_ + (ny+2)))/dx*(O*kro/voi+W*krw/vwi);
      ve = D*ve;
  else
      ve = 0;
  end
  
  %} North velocity
 if ~isnan(P_(cellnum_ - 1))
      I = P(cellnum) < P(cellnum - 1);
      kro = relperm_oil(cellnum - 1*I);
      krw = relperm_wat(cellnum - 1*I);
      vn = 0.00633/phi(cellnum)*2*perm(cellnum)*perm(cellnum-1)/(perm(cellnum)+perm(cellnum-1)) ...
          *(P_(cellnum_) - P_(cellnum_ - 1))/dx*(O*kro/voi+W*krw/vwi);
      vn = D*vn;
  else
      vn = 0;
  end
  
  %} South velocity
  if ~isnan(P_(cellnum_ + 1))
      I = P(cellnum) < P(cellnum + 1);
      kro = relperm_oil(cellnum + 1*I);
      krw = relperm_wat(cellnum + 1*I);
      vs = 0.00633/phi(cellnum)*2*perm(cellnum)*perm(cellnum+1)/(perm(cellnum)+perm(cellnum+1)) ...
          *(P_(cellnum_ + 1) - P_(cellnum_))/dx*(O*kro/voi+W*krw/vwi);
      vs = D*vs;
  else
      vs = 0;
  end
  
  fun_vx = @(a)vw+(ve-vw)/dx*(a-x0);
  fun_vy = @(b)vs+(vn-vs)/dx*(y0-b);
  if feval(fun_vx, x) >= 0
      vxe = ve;
  else
      vxe = vw;
  end
  if feval(fun_vx, x)*vxe < 0
      vxe = 0;
  end
  if all([vw, ve] == 0)
      del_t_x = inf;
  else
      del_t_x = abs(1/((ve - vw) / dx) * log(vxe/feval(fun_vx, x)));
  end
%   if isnan(del_t_x)
%       del_t_x = 0;
%   end
  
  if feval(fun_vy, y) >= 0
      vye = vn;
  else
      vye = vs;
  end
  if feval(fun_vy, y)*vye < 0
      vye = 0;
  end
%   if feval(fun_vy, y) == 0
%       keyboard()
%   end
  if all([vs, vn] == 0)
      del_t_y = inf;
  else
      del_t_y = abs(1/((vn - vs) / dy) * log(vye/feval(fun_vy, y)));
  end
%   if isnan(del_t_y)
%       del_t_y = 0;
%   end
  
if any([del_t_y, del_t_x] ~= inf)
    
    [del_t, I] = min([del_t_x, del_t_y]);
    
    if I == 2  % y쪽으로 나감
        
        x = 1/((ve - vw) / dx)*(feval(fun_vx, x)*exp((ve - vw) / dx * del_t) - vw)+x0;
        
        if feval(fun_vy, y) >= 0   % 북쪽으로 나감
            y = y0 - dy;
            cellnum = cellnum - 1;
        else                       % 남쪽으로 나감
            y = y0;
            cellnum = cellnum + 1;
        end
        
    else       % x쪽으로 나감
        
        y_ = ( 1/((vn - vs) / dy)*(feval(fun_vy, y)*exp((vn - vs) / dy * del_t) - vs) + (dy*ny-y0) );
        if y_ < 0.1 % del_t가 너무커서 y쪽으로 나가는 것처럼 보일때
            y = ny*dy - 0.1;
        else
            y = ny*dy - y_;
        end
        
        
        if feval(fun_vx, x) >= 0    % 동쪽
            x = x0 + dx;
            cellnum = cellnum + ny;
        else                        % 서쪽
            x = x0;
            cellnum = cellnum - ny;
        end
        
    end
else
    del_t = inf;
end
 


end
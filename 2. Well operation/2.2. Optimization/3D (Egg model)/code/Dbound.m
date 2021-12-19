function [mindist] = Dbound(pos, active)
global dx dy

mindist = [];
for i = 1:size(pos,2)/2
    n = find(active(pos(1,2*i):end,pos(1,2*i-1)) == 0, 1) - 1;
    s = find(flipud(active(1:pos(1,2*i),pos(1,2*i-1))) == 0, 1) - 1;
    e = find(active(pos(1,2*i),pos(1,2*i-1):end) == 0, 1) - 1;
    w = find(fliplr(active(pos(1,2*i),1:pos(1,2*i-1))) == 0, 1) - 1;
    mindist = [mindist; min([dy*n,dy*s,dx*e,dx*w].^2)];
end

end
function d = dominates(x,y)
    v.x = x(:,size(x,2));
    v.y = y(:,size(y,2));
    f.x = x(:,1:size(x,2)-1);
    f.y = y(:,1:size(y,2)-1);
    d   = zeros(size(x,1),1);
    for i = 1:size(x,1)
        if all([v.x(i,1),v.y(i,1)]>0,2)
            d(i,1) = (v.x(i,1) < v.y(i,1));
        elseif any([v.x(i,1),v.y(i,1)]==0,2) && ~all([v.x(i,1),v.y(i,1)]==0,2)
            d(i,1) = (v.x(i,1) == 0);
        else
            d(i,1) = all(f.x(i,:)<=f.y(i,:),2) & any(f.x(i,:)<f.y(i,:),2);
        end
    end
end
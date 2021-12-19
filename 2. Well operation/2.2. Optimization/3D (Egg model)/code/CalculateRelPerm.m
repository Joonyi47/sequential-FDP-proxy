function [RelPerm] = CalculateRelPerm(table, sat, phase)

if strcmp(phase, 'w')
    RelPerm = interp1(table(:,1), table(:,2), sat, 'linear', 'extrap');
    if any(sat < min(table(:,1)))
        RelPerm(sat < min(table(:,1))) = table(1,2);
    elseif any(sat > max(table(:,1)))
        RelPerm(sat > max(table(:,1))) = table(end-1,2);
    end
elseif strcmp(phase, 'o')
    RelPerm = interp1(table(:,1), table(:,3), 1-sat, 'linear', 'extrap');
    if any(sat < 1 - max(table(:,1)))
        RelPerm(sat < 1 - max(table(:,1))) = table(end-1,3);
    elseif any(sat > 1 - min(table(:,1)))
        RelPerm(sat > 1 - min(table(:,1))) = table(1,3);
    end
end

end
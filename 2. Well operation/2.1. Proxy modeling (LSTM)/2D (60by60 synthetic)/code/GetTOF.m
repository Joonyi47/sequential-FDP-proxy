function [so] = GetTOF(filename, Nt, keyword)

fclose all;
if fix(Nt/10) > 0
    so=get_griddata([filename '.F00' int2str(Nt)] ,keyword);
else
    so=get_griddata([filename '.F000' int2str(Nt)] ,keyword);
end
end
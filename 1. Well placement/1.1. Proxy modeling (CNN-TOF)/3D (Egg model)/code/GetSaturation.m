function [so] = GetSaturation(filename, Nt, keyword)

Restart_converter_writting(filename,Nt);
% copyfile('Restart_converter.log', ['PSO' int2str(iter) '/particle' int2str(best) '/' 'Restart_converter.log']);
% currentfolder=pwd;
% cd ([currentfolder '\' 'PSO' int2str(iter) '\' 'particle' int2str(best)]);
dos('$convert < Restart_converter.log > NUL');                               % 변환시키기.
fclose all;
if fix(Nt/10) > 0
    so=get_griddata([filename '.F00' int2str(Nt)] ,keyword);
else
    so=get_griddata([filename '.F000' int2str(Nt)] ,keyword);
end
end
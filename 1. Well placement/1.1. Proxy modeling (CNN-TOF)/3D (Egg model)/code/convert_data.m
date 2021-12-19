function [] = convert_data(filename, Nt)

Restart_converter_writting(filename,Nt);
% copyfile('Restart_converter.log', ['PSO' int2str(iter) '/particle' int2str(best) '/' 'Restart_converter.log']);
% currentfolder=pwd;
% cd ([currentfolder '\' 'PSO' int2str(iter) '\' 'particle' int2str(best)]);
[a,b] = dos('$convert < Restart_converter.log > NUL');                              % 변환시키기.
if size(b,2) > 20
    keyboard()
end
end
function [Graph_Mat, TOF, P, tcpu] = SLsimulate(pos, type, directory, datafile, permfile, ne)
    global Np  ...
           posfile constfile ...
           dtstep pmax ...
           slstep nsteps

    Np = size(pos,1);

    pos_fit = [];
    pos_vio = [];
    
    
    %% for parallel computing %%
    [home] =cd(directory);
    delete('*.S*', '*.F*', '*.X*');
    
    fid = fopen([datafile '.DATA'], 'r');
    total=100;
    if fid == -1
        disp('File open not successful')
    end
    Ecl=cell(300,300);
    i=1;
    while ~feof(fid)
        Ecl{i} = fgetl(fid);
        i=i+1;
    end
    row=find(strcmp(Ecl,'INCLUDE')==1,total);
    for j = 1:Np
        fid = fopen([datafile '_' int2str(j) '.DATA'], 'w');
        for i = 1:i
            if i == row(2,1)+1
                fprintf(fid, '%s\n', [char(39)  permfile '.DATA'  char(39) ' /']);
            elseif i == row(3,1)+1
                fprintf(fid, '%s\n', [char(39)  posfile '_' int2str(j) '.DATA'  char(39) ' /']);
            elseif i == row(4,1)+1
                fprintf(fid, '%s\n', [char(39)  constfile '_' int2str(j) '.DATA' char(39) ' /']);
            else
                fprintf(fid, '%s\n', Ecl{i});
            end
        end
        fclose(fid);
    end
    
    fclose('all');
    
   
  
    %%  Streamline simulation
    Graph_Mat   = {};
    TOF         = {};
    tcpu        = [];

        parfor j = 1:Np
            change_tstep([constfile '_' int2str(j)], nsteps, slstep/nsteps);
        end
        
        parfor j = 1:Np
           [~,~] = dos(['C:\ecl\2009.1\bin\pc\frontsim.exe ' datafile '_' int2str(j) ' > NUL']);
        end

        parfor j = 1:Np
            change_tstep([constfile '_' int2str(j)], pmax/dtstep, dtstep);
        end
        
        TOF = mat2cell(zeros(1,Np), 1,ones(1,Np));
        
        Restart_converter_writting_batch(datafile, nsteps, 1:Np);
        [~,~] = dos('$convert < Restart_converter.log > NUL');                              % 변환시키기.
        
        TOF_beg = [];
        TOF_end = [];
        PRESS   = [];
        
        parfor j= 1:Np
            
            TOF_beg = [TOF_beg, GetTOF([datafile '_' int2str(j)], nsteps, 'TIME_BEG')];
            TOF_end = [TOF_end, GetTOF([datafile '_' int2str(j)], nsteps, 'TIME_END')];
            PRESS   = [PRESS, GetTOF([datafile '_' int2str(j)], nsteps, 'PRESSURE')];

        end

    TOF = [{TOF_beg}, {TOF_end}];
    P   = {PRESS};
    
   %% delete all file except *.DATA, *.RSM

%     mkdir(['Restart files for field ' int2str(ne)]);
%     copyfile(fullfile('*.F*'), ['Restart files for field ' int2str(ne)]);
    
    delete('*.S*', '*.F*', '*.X*');
    
    fclose('all');
    cd(home);
    
    


end
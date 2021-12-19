function [Graph_Mat, TOF] = SLsimulate(pos, type, directory, datafile, permfile, ne)
    global Np  ...
           posfile constfile ...
           nx ny ...
           dstep tstep pmax ...
           parameters
    
    Np = size(pos,1);
    
    pos_fit = [];
    pos_vio = [];
    
    %% for parallel computing %%
    [home] =cd(directory);
    
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
        for i = 1:134
            if i == row(1,1)+1
                fprintf(fid, '%s\n', [char(39)  permfile '.DATA'  char(39) ' /']);
            elseif i == row(2,1)+1
                fprintf(fid, '%s\n', [char(39)  posfile '_' int2str(j) '.DATA'  char(39) ' /']);
            elseif i == row(3,1)+1
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

        %%%% frontsim %%%%
        parfor j = 1:Np
            change_tstep([constfile '_' int2str(j)], 1, tstep);
%             change_tstep([constfile '_' int2str(j)], 1, 10);
        end
        
        
        parfor j = 1:Np
           [~,~] = dos(['C:\ecl\2009.1\bin\pc\frontsim.exe ' datafile '_' int2str(j) ' > NUL']);
        end
                
        parfor j = 1:Np
            change_tstep([constfile '_' int2str(j)], pmax/dstep, dstep);
        end
    
        TOF = mat2cell(zeros(1,Np), 1,ones(1,Np));
        
        Restart_converter_writting_batch(datafile, 1, 1:Np);
        [~,~] = dos('$convert < Restart_converter.log > NUL');                              % 변환시키기.
        
        TOF_beg = [];
        TOF_end = [];
        
        parfor j= 1:Np
            
            TOF_beg = [TOF_beg, GetTOF([datafile '_' int2str(j)], 1, 'TIME_BEG')];
            TOF_end = [TOF_end, GetTOF([datafile '_' int2str(j)], 1, 'TIME_END')];

        end

%         %%%% impes & ss %%%%
%         params = parameters;
%         
%         perm = GetPermx([permfile '.DATA']);
%         
%         P = [];
%         So = [];
%          
%         parfor j = 1:Np
%             
%             [p, so] = IMPES(pos(j,:), type(j,:), perm, params.swof, params);
%             P = [P, p(:,end)];
%             So = [So, so(:, end)];
%             
%         end
%         
%         TOF_beg = [];
%         TOF_end = [];
%         
%         parfor j = 1:Np
%             
%             [tbeg, tend] = StreamSim(pos(j,:), type(j,:), perm, params.swof, P(:,j), So(:,j), params);
%             TOF_beg = [TOF_beg, tbeg];
%             TOF_end = [TOF_end, tend];
%             
%         end
%         
%         parfor j = 1:Np
%             
%             TOF_beg(:,j) = reshape(reshape(TOF_beg(:,j), nx, ny)', nx*ny, 1);
%             TOF_end(:,j) = reshape(reshape(TOF_end(:,j), nx, ny)', nx*ny, 1);
%             
%         end
%         parfor j= 1:Np
% %             TOF{j} = GetTOF([datafile '_' int2str(j)], 1, 'TIME_BEG') ...
% %                    + GetTOF([datafile '_' int2str(j)], 1, 'TIME_END'); 
% %             TOF{j} = reshape(log(TOF{j}), nx, ny)';
%             
%             source = mat2cell(pos(j,repelem(~type(j,:),1,2)), 1, 2*ones(1,length(find(type(j,:) == 0))));
%             sink   = mat2cell(pos(j,repelem(logical(type(j,:)),1,2)), 1, 2*ones(1,length(find(type(j,:) == 1))));
%                
% %             temp = reshape(log(TOF_beg(:,j)), nx, ny)';
% %             graph = transform_to_graph(source, sink, temp);
% %             Graph_Mat_beg = [Graph_Mat_beg, {graph}];
% %             
% %             temp = reshape(log(TOF_end(:,j)), nx, ny)';
% %             graph = transform_to_graph(sink, source, temp);
% %             Graph_Mat_end = [Graph_Mat_end, {graph}];
%             
%             graph = transform_to_graph(source, sink, reshape((TOF_beg(:,j) + TOF_end(:,j)), nx, ny)');
%             Graph_Mat = [Graph_Mat, {graph}];
%                                   
%             if mod(j, (Np/100) ) == 0
%                 disp([num2str(j/Np * 100) '%']);
%             end
% 
%         end
    % delete all file except *.DATA, *.RSM
    TOF = [{TOF_beg}, {TOF_end}];
%     mkdir(['Restart files for field ' int2str(ne)]);
%     copyfile(fullfile('*.F*'), ['Restart files for field ' int2str(ne)]);
    
    delete('*.S*', '*.F*', '*.X*');
    
    fclose('all');
    cd(home);
    
    


end
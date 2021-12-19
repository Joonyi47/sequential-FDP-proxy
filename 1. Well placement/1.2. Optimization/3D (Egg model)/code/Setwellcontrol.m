function [] = Setwellcontrol(POS, filename)
    global N ...
           dtstep pmax
    
    loc  = POS(1,1:2*N);
    type = POS(1,2*N+1:3*N);
    wset = POS(1,3*N+1:end);
    

    fid = fopen(filename, 'w');
    
        fprintf(fid, '\n%s\n', 'WCONPROD');
        for i = 1:N
           if type(i) >= 1/3 
%                fprintf(fid, '%s %s %s %s %s %s %d', [char(39) 'PROD' int2str(i) char(39)] ,'OPEN','BHP','3*','700','1*',wset(2) - (type(i) - 1/3) *3/2 * (wset(2) - wset(1)));
               fprintf(fid, '%s %s %s %s %s %s %d', [char(39) 'PROD' int2str(i) char(39)] ,'OPEN','BHP','3*','1*','1*',wset(2) - (type(i) - 1/3) *3/2 * (wset(2) - wset(1)));
               fprintf(fid, '/\n');
           end
        end
        fprintf(fid, '/\n');
        fprintf(fid, '\n%s\n', 'WCONINJE');
        for i = 1:N
            if type(i) <= -1/3
%                 fprintf(fid, '%s %s %s %s %d %s %s', [char(39) 'INJECT' int2str(i) char(39)] ,'WATER','OPEN','RATE',wset(3) + (- type(i) - 1/3) *3/2 * (wset(4) - wset(3)),'1*','800');
                fprintf(fid, '%s %s %s %s %s %s %d', [char(39) 'INJECT' int2str(i) char(39)] ,'WATER','OPEN','BHP', '1*', '1*', wset(3) + (- type(i) - 1/3) *3/2 * (wset(4) - wset(3)));
                fprintf(fid, '/\n');          
            end
        end
        fprintf(fid, '/\n');
        fprintf(fid, '\n%s\n', 'TSTEP');
        fprintf(fid, '%s%s /\n', [int2str(pmax/dtstep) '*'],int2str(dtstep));
        
    fclose(fid);
end
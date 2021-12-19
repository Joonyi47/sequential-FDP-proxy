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
               fprintf(fid, '%s %s %s %s %s %d %d', ['P' int2str(i)] ,'1*','BHP','5000','4*', wset(2) - (type(i) - 1/3) *3/2 * (wset(2) - wset(1)));
               fprintf(fid, '/\n');
           end
        end
        fprintf(fid, '/\n');
        fprintf(fid, '\n%s\n', 'WCONINJE');
        for i = 1:N
            if type(i) <= -1/3
                fprintf(fid, '%s %s %s %s %s %s %d', ['I' int2str(i)] ,'WATER','1*','BHP','5000','1*',wset(3) + (- type(i) - 1/3) *3/2 * (wset(4) - wset(3)));
%                 fprintf(fid, '%s %s %s %s %d %s %s', ['I' int2str(i)] ,'WATER','1*','RATE', wset(2),'1*', '1*');
                fprintf(fid, '/\n');          
            end
        end
        fprintf(fid, '/\n');
        fprintf(fid, '\n%s\n', 'TSTEP');
        fprintf(fid, '%s%s /\n', [int2str(pmax/dtstep) '*'],int2str(dtstep));
        
    fclose(fid);
end
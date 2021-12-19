function [] = Setwellpos(POS, filename)
    global N 
    
    loc  = POS(1,1:2*N);
    type = POS(1,2*N+1:end);
    
    fid = fopen(filename, 'w');
    fprintf(fid, '\n%s\n', 'WELSPECS');
    for i = 1:N
        if type(i) >= 1/3 % prod.
            fprintf(fid, '%s %s %d %d %s %s ', [char(39) 'PROD' int2str(i) char(39)] ,'1',loc(2*i-1),loc(2*i), '1*', 'OIL');
            fprintf(fid, '/\n');
        end
    end
    for i = 1:N
        if type(i) <= -1/3 % injec.
            fprintf(fid, '%s %s %d %d %s %s ', [char(39) 'INJECT' int2str(i) char(39)],'1',loc(2*i-1),loc(2*i), '1*', 'WATER');
            fprintf(fid, '/\n');
        end
    end
    fprintf(fid, '/\n');
    fprintf(fid, '\n%s\n', 'COMPDAT');
    for i = 1:N
        if type(i) >= 1/3
            fprintf(fid, '%s %s %d %d %s %s %s %s %s ', [char(39) 'PROD' int2str(i) char(39)],'2*', 1, 2 ,'OPEN','2*','0.2','1*','0');
            fprintf(fid, '/\n');
        end
    end
    for i = 1:N
        if type(i) <= -1/3
            fprintf(fid, '%s %s %d %d %s %s %s %s %s ', [char(39) 'INJECT' int2str(i) char(39)],'2*', 1, 2 ,'OPEN','2*','0.2','1*','0');
            fprintf(fid, '/\n');
        end
    end
    fprintf(fid, '/\n');
fclose(fid);
end
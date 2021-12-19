function [] = Setwellpos(POS, filename)
    global N Nstep
    
    loc  = POS(1,1:2*N);
    type = POS(1,2*N+1:end);
    
    fid = fopen(filename, 'w');
    fprintf(fid, '\n%s\n', 'WELSPECS');
    for i = 1:N
        if type(i) >= 1/3 % prod.
            fprintf(fid, '%s %s %d %d %s %s %s %s ', ['P' int2str(i)] ,'ALL',loc(2*i-1),loc(2*i), '1*', 'LIQ', '3*', 'NO');
            fprintf(fid, '/\n');
        end
    end
    for i = 1:N
        if type(i) <= -1/3 % injec.
            fprintf(fid, '%s %s %d %d %s %s %s %s ', ['I' int2str(i)],'ALL',loc(2*i-1),loc(2*i), '1*', 'WATER', '3*', 'NO');
            fprintf(fid, '/\n');
        end
    end
    fprintf(fid, '/\n');
    fprintf(fid, '\n%s\n', 'COMPDAT');
    for i = 1:N
        if type(i) >= 1/3
            fprintf(fid, '%s %d %d %d %d %s %s %s %s %s %s %s %s ', ['P' int2str(i)],loc(2*i-1),loc(2*i),1,1,'1*','1*','1*','1','1*','1*','1*','Z');
            fprintf(fid, '/\n');
        end
    end
    for i = 1:N
        if type(i) <= -1/3
            fprintf(fid, '%s %d %d %d %d %s %s %s %s %s %s %s %s ', ['I' int2str(i)],loc(2*i-1),loc(2*i),1,1,'1*','1*','1*','1','1*','1*','1*','Z');
            fprintf(fid, '/\n');
        end
    end
    fprintf(fid, '/\n');
fclose(fid);
end
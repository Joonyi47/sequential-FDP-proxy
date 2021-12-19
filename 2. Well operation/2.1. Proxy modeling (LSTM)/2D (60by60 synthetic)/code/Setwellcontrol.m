function [] = Setwellcontrol(POS, filename)
global N Nstep ...
       dstep tstep pmax ...
       opt_type

loc  = POS(1,1:2*N);
type = POS(1,2*N+1:2*N+4*N);
wset = POS(1,2*N+4*N+1:end);

type = reshape(type, 4, N);
schd = zeros(N, Nstep);
for i = 1:N
   schd(i, ceil(Nstep * type(1,i)):end) = 1;
end
fun = @(x,y,z,b,t)x*(1+cos(y*(1-t)+z))+b;

fid = fopen(filename, 'w');
for j = 1:Nstep
    fprintf(fid, '\n%s\n', 'WCONPROD');
    for i = 1:N
        wcont = feval(fun, (wset(2) - wset(1))/2*type(2,i), type(3,i), type(4,i), wset(1), [1/Nstep:1/Nstep:1]);
        if opt_type(i) >= 1/3            
            if schd(i,j) == 0
                fprintf(fid, '%s %s %s %s %s %d %d', ['P' int2str(i)] ,'SHUT','BHP','6000','4*', wset(2));
            else
%                 fprintf(fid, '%s %s %s %s %s %d %d', ['P' int2str(i)] ,'1*','BHP','5000','4*', wset(2) - (type(i,j) - 1/3) *3/2 * (wset(2) - wset(1)));
                fprintf(fid, '%s %s %s %s %s %d %d', ['P' int2str(i)] ,'1*','BHP','6000','4*', wcont(1,j));
                %                fprintf(fid, '%s %s %s %s %d %d', ['P' int2str(i)] ,'1*','LRAT', '3*', wset(2) - (type(i) - 1/3) *3/2 * (wset(2) - wset(1)));
            end
            fprintf(fid, '/\n');
        end
    end
    fprintf(fid, '/\n');
    fprintf(fid, '\n%s\n', 'WCONINJE');
    for i = 1:N
        wcont = feval(fun, (wset(4) - wset(3))/2*type(2,i), type(3,i), type(4,i), wset(3), [1/Nstep:1/Nstep:1]);
        if opt_type(i) <= -1/3
            if schd(i,j) == 0
                fprintf(fid, '%s %s %s %s %s %s %d', ['I' int2str(i)] ,'WATER','SHUT','BHP','6000','1*', wset(3));
            else
%                 fprintf(fid, '%s %s %s %s %s %s %d', ['I' int2str(i)] ,'WATER','1*','BHP','5000','1*',wset(3) + (- type(i,j) - 1/3) *3/2 * (wset(4) - wset(3)));
                fprintf(fid, '%s %s %s %s %s %s %d', ['I' int2str(i)] ,'WATER','1*','BHP','6000','1*', wcont(1,j));
                %                 fprintf(fid, '%s %s %s %s %d %s %s', ['I' int2str(i)] ,'WATER','1*','RATE', wset(3) + (- type(i) - 1/3) *3/2 * (wset(4) - wset(3)),'1*', '1*');
            end
            fprintf(fid, '/\n');
        end
    end
    fprintf(fid, '/\n');
    fprintf(fid, '\n%s\n', 'TSTEP');
    fprintf(fid, '%s%s /\n\n', [int2str(dstep/tstep) '*'],int2str(tstep));
end

fclose(fid);
end
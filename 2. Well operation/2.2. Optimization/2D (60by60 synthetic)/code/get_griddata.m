function [aa, bb]=get_griddata(filename,property)

fid=fopen(filename);
                             
keyword=property;
% keyword='10';

% 파일을 열고 연 파일을 스캐닝하기 시작.

i=0;
aa=[];
bb=0;

%{
while (~feof(fid))
    i=i+1;
    tline=fgets(fid);    
       
    
end
%}   

% 일단 이렇게 하면 i에는 전체의 줄이 다 잡히고
% tline에는 다음의 값이 기록된다.

while(~feof(fid))
    tline = fgetl(fid);
    index=strread(tline,'%s');         
    
    [location e]=regexp(index,keyword);
    pp=prod(cellfun(@isempty,location)*1);
    
    if pp==0 & location{1}==2                         %-------- yata! find!!
        
        a=str2num(index{2});
        if isempty(a)==1
            a=str2num(index{3});
        end
        
        while bb~=a
        
            tline = fgetl(fid);
            index=strread(tline,'%s');
        
            for ii=1:length(index)            
                aa=[aa str2num(index{ii})];
            end
            
            bb=length(aa);
            
        end
        
        break
    end


end
aa=aa';


%end

fclose('all');

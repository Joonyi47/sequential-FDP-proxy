function [aa, bb]=get_griddata(filename,property)

fid=fopen(filename);
                             
keyword=property;
% keyword='10';

% ������ ���� �� ������ ��ĳ���ϱ� ����.

i=0;
aa=[];
bb=0;

%{
while (~feof(fid))
    i=i+1;
    tline=fgets(fid);    
       
    
end
%}   

% �ϴ� �̷��� �ϸ� i���� ��ü�� ���� �� ������
% tline���� ������ ���� ��ϵȴ�.

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

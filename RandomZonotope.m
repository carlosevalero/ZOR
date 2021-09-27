function Z0=RandomZonotope(Kmax,NMD,NMO)

%% Creating Random zonotopes
Z0=cell(NMD,NMO,Kmax);
for k=1:Kmax
    for i=1:NMD       %number of dimension
        for j=1:NMO   %number of order
            o=(i+1)*(j+1);
%             T=zeros(i+1,o);
%             for l=1:o
%                 T(:,l)=rand(i+1,1);
%             end
            T=rand(i+1,o);
            alpha=randi(60,1,o);
            T=alpha.*T;
            c=zeros(i+1,1);
            Z0{i,j,k} = zonotope([c, T]);
        end
    end
end
path = pwd ;   % mention your path
myfolder = strcat('Test','d',string(NMD),'o',string(NMO),'s',string(Kmax));
myfolder=char(myfolder);
folder = mkdir([path,filesep,myfolder]);
Cadena=strcat(myfolder,'/','RandomZonotopes','d',string(NMD),'o',string(NMO),'s',string(Kmax));
save(Cadena,'Z0','-v7.3');
end
%% New Zonotope Order reduction. Root file for tests.
%Require CORA toolbox
clear all;
close all;
clc;
%% Creating Random Zonotopes
Kmax=10;%Test(z,1);  %Maximum number of zonotopes with NMD and NMO
NMD=15;%Test(z,2);   %Number Maximum of Dimension
NMO=15;%Test(z,3);   %Number Maximum of Orders
path = pwd ;   % mention your path
myfolder = strcat('Test','d',string(NMD),'o',string(NMO),'s',string(Kmax));
Cadena=strcat(myfolder,'/RandomZonotopes','d',string(NMD),'o',string(NMO),'s',string(Kmax),'.mat');
if isfile(Cadena)
    load(Cadena);
else
    Z0=RandomZonotope(Kmax,NMD,NMO);
end
%% Defining the method in CORA
Method={
    %     'cluster'
    'combastel'
    % %     'constOpt'
    'girard'
    % %     'methA'
    % %     'methB'
    % %     'methC'
    'pca'
    'scott'
    % %     'redistribute'
    };
Nm=length(Method); %number of methods
Cadena=strcat(myfolder,'/ReducedZonotopes','d',string(NMD),'o',string(NMO),'s',string(Kmax),'nm',string(Nm),'.mat');
if isfile(Cadena)
    load(Cadena);
else
    Z=cell(Nm+1,1);
    CPUT=cell(Nm+1,1);
    V=cell(Nm+1,1);
    for i=1:Nm+1
        Z{i}=cell(NMD-1,NMO,Kmax);
        CPUT{i}=cell(NMD-1,NMO,Kmax);
        V{i}=cell(NMD-1,NMO,Kmax);
    end
    zz=zeros(2,1);zz1=zeros(2,1);
    %% Main Code
    for i=1:NMD            %number of dimension
        for j=1:NMO        %number of generators
            for k=1:Kmax   %number of different zonotopes
                T=Z0{i,j,k}.generators;
                tic;
                Z{Nm+1}{i,j,k}=ZORMethod(T,i+1,2);
                %                     Z{Nm+1}{i,j,k}=optimalzonoReducProj(Z0{i,j,k}.generators,i+1);
                CPUT{Nm+1}{i,j,k}=toc;
                V{Nm+1}{i,j,k}=(2^i)*abs(det(Z{Nm+1}{i,j,k}));
                for l=1:Nm
                    tic
                    Z{l}{i,j,k}=reduce(Z0{i,j,k},Method{l},1);
                    CPUT{l}{i,j,k}=toc;
                    V{l}{i,j,k}=(2^i)*abs(det(Z{l}{i,j,k}.generators));
                end
                if V{3}{i,j,k}<=V{5}{i,j,k}
                    zz(1)=zz(1)+1;
                end
                if V{5}{i,j,k}<=V{3}{i,j,k}
                    zz1(1)=zz1(1)+1;
                end
                if V{4}{i,j,k}<=V{5}{i,j,k}
                    zz(2)=zz(2)+1;
                end
                if V{5}{i,j,k}<=V{4}{i,j,k}
                    zz1(2)=zz1(2)+1;
                end
            end
        end
    end
    Cadena=strcat(myfolder,'/ReducedZonotopes','d',string(NMD),'o',string(NMO),'s',string(Kmax),'nm',string(Nm),'.mat');
    save(Cadena,'CPUT','V','Z','-v7.3','zz1','Nm','NMD','NMO','Kmax','myfolder');
end
% end
disp(zz1);
%% Pre-Processing Data
% Normalizing Volume.The method with the major volume at each
VolumeProm=cell(Nm+1,1);
NormalVolume=cell(Nm+1,1);
VolumeScott=cell(Nm+1,1);
VolumePCA=cell(Nm+1,1);
VolumeApproach=cell(Nm+1,1);
for i=1:NMD
    for j=1:NMO
        Volmax=0;
        for l=1:5
            AV=0;
            for k=1:Kmax
                AV=AV+V{l}{i,j,k}/Kmax;
            end
            if AV>Volmax
                Volmax=AV;
            end
            VolumeProm{l}(i,j)=AV;
        end
        for l=1:5
            NormalVolume{l}(i,j)=VolumeProm{l}(i,j)/Volmax;
            VolumeScott{l}(i,j)=min(5,VolumeProm{l}(i,j)/VolumeProm{4}(i,j));
            VolumePCA{l}(i,j)=min(5,VolumeProm{l}(i,j)/VolumeProm{3}(i,j));
            VolumeApproach{l}(i,j)=min(5,VolumeProm{l}(i,j)/VolumeProm{5}(i,j));
        end
    end
end

%% Graphs normalized volume
Method={
    'Combastel'
    'Girard'
    'PCA'
    'Scott'
    'Author Approach'
    };

%Regarding of the number of possible graphs we will show only five.
index=[1 2 round(2*NMD/5)-1  round(3*NMD/5)-1 round(4*NMD/5)-1 NMD-1];
% index=1:15;
symbol={'x-','s-','d-','+-','o-'};
wideline=[0.5,0.5,0.5,1.5,2];
t=2:NMD+1;
for i=1:length(index)
    cadena=strcat(myfolder,'/rVAppR',string(index(i)+1),'d',string(NMD),'o',string(NMO),'s',string(Kmax),'.fig');
    if ~isfile(cadena)
        
        figure
        for l=1:Nm+1
            plot(t,VolumeApproach{l}(index(i),:),symbol{l},'LineWidth',wideline(l));
            hold on;
        end
        grid on;
        legend(Method);
        xlabel('Zonotope Order','interpreter','latex');
        ylabel('$ min (5,\frac{V_i}{V_{Ap}})$','interpreter','latex');
        cadena=strcat('Zonotope Order Reduction in $R^{',string(index(i)+1),'}$');
        title(cadena,'interpreter','latex');
        cadena=strcat(myfolder,'/rVAppR',string(index(i)+1),'d',string(NMD),'o',string(NMO),'s',string(Kmax));
        saveas(gcf,cadena,'fig');
        saveas(gcf,cadena,'eps');
    else
        openfig(cadena);
    end
end
%% Computational Time
% % %Checking CPUT
% % % CPUT{M}{D,G,Z} M-Method D-Dimension G-Generators Z-Zonotopes
ACPUT=cell(Nm+1,1);
AV=cell(Nm+1,1);
CPU=zeros(NMD,NMO);
AVm=zeros(NMD,NMO);
for i=1:Nm+1
    for j=1:NMD
        for k=1:NMO
            c=0;
            d=0;
            for l=1:Kmax
                c=c+CPUT{i}{j,k,l}/Kmax;
                d=d+V{i}{j,k,l}/Kmax;
            end
            CPU(j,k)=c;
            AVm(j,k)=d;
        end
    end
    ACPUT{i}=CPU;
    AV{i}=AVm;
end
for i=1:length(index)
    cadena=strcat(myfolder,'/CPUR',string(index(i)+1),'d',string(NMD),'o',string(NMO),'s',string(Kmax),'.fig');
    if ~isfile(cadena)
        
        figure
        for j=1:Nm+1
            plot(t,ACPUT{j}(index(i),:),symbol{j},'LineWidth',wideline(j));
            hold on;
        end
        grid on;
        legend(Method);
        xlabel('Zonotope Order','interpreter','latex');
        ylabel('CPU Time','interpreter','latex');
        cadena=strcat('Zonotope Order Reduction in $R^{',string(index(i)+1),'}$');
        title(cadena,'interpreter','latex');
        cadena=strcat(myfolder,'/CPUR',string(index(i)+1),'d',string(NMD),'o',string(NMO),'s',string(Kmax));
        saveas(gcf,cadena,'fig');
        saveas(gcf,cadena,'eps');
    else
        openfig(cadena);
    end
end
% 

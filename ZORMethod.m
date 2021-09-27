function Tn=ZORMethod(T,g)
% Tn=ZORMethod(T,g)
%    This function computes an matrix Tn with g columns, taking into
%    account the columms of the matrix T using the Gram-Schmidt Algorithm
%    T Matrix to reduce (n x p)
%    g Number of generators wished (1x1).
[n,m]=size(T);
if g<n || g > m
    msgbox("The number of generators wished is not feasible","ERROR",'error');
end
%% Projection into the maximum angle
    a=vecnorm(T);
    [~,index]=sort(a,'descend');
    T=T(:,index);
    T3=T;
    volF=0;
    for i=1:1
        T=T3;
        Tr=zeros(n);
        Tr(:,1)=T(:,i);
        T(:,i)=[];
        taux=Tr./vecnorm(Tr);
        %% Gram-Schmidt Criteria
        In=eye(n);
        p=1;
        while p<n
            %Orthogonal projection of all vectors over Taux
            c1=taux(:,1:p)*taux(:,1:p)';
            b1=(In-c1);
            aux=b1*T;
            %Reorder
            Cr=vecnorm(aux);
            [~,index]=sort(Cr,'descend');
            T=T(:,index);
            aux=aux(:,index);
            Tr(:,p+1)=T(:,1);
            taux(:,p+1)=aux(:,1)/norm(aux(:,1));
            T(:,1)=[];
            p=p+1;
        end
        volTn=abs(det(Tr));
        if volTn>volF
            Trn=Tr;
            volF=volTn;
        end
    end
    Tr=Trn;
    a=vecnorm(T);
    [~,index]=sort(a,'ascend');
    T=T(:,index);
    %% Linear Combination
    alpha= ones(n,1)+sum(abs(Tr\T(:,1:end-g+n)),2);
    %% Selection of the order
    T(:,1:end-g+n)=[];
    Tn=[alpha'.*Tr T];
end
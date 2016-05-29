clearvars; clc; load data.mat; global J T; data=log(D); 
for i=1:size(data,1)
    for j=1:size(data,2)
        for k=1:size(data,3)
            ndata(k,j,i)=data(i,j,k);
        end
    end
end
data=ndata;
clear i j k D;
%%
I=size(data,1); J=size(data,2)-1; K=size(data,3);
f_plus=data(:,2:end,:); f=data(:,1:end-1,:); D=f_plus-f; T=1:J;
Mean_vector=zeros(K,J); Mean_m=zeros(I,J,K);
for k=1:K
    Mean_vector(k,:)=mean(D(:,:,k));
    Mean_m(:,:,k)=repmat(mean(D(:,:,k)),I,1);
end
CUT=zeros(1,K);
for k=1:K
    CUT(k)=detector(D(:,:,k),4);
    CUT1(k)=detector1(D(:,:,k));
end
%%
f_plus=D(:,2:end,:); f=D(:,1:end-1,:); SE_D=f_plus-f;
for II=1:K
    for l=1:11
        SCP=detector(SE_D(:,:,II),l);
        RE(l)=SCP;
    end
    [Test beta_0 V_B]=wald_drevetive(D(:,:,II),min(RE)-1);
    UCL=beta_0+1.96*V_B;
    LCL=beta_0-1.96*V_B;
    disp([min(RE)-2 beta_0 Test LCL UCL]);
    Result(II,:)=[min(RE)-2 Test LCL UCL];
    %%
    figure;
    subplot(2,2,1)
    J=size(ndata,2);
    T=1:J;
    for i=1:I
        p=plot(T,ndata(i,:,II));
        set(p,'Color',rand(3,1));
        hold on;
    end
    DA=ndata(:,:,II);
    p=plot(T,mean(DA));
    set(p,'LineWidth',2,'color','black','linestyle','--');
    plot([Result(II,1)+1,Result(II,1)+1],[min(DA(:)),max(DA(:))],'--','LineWidth',2);
    title(sprintf('original gene %d',II));
    %%
    T=1:J-1;
    subplot(2,2,2)
    for i=1:I
        p=plot(T,D(i,:,II));
        set(p,'Color',rand(3,1));
        hold on;
    end
    DA=D(:,:,II);
    p=plot(T,mean(DA));
    set(p,'LineWidth',2,'color','black','linestyle','--');
    plot([Result(II,1)+1,Result(II,1)+1],[min(DA(:)),max(DA(:))],'--','LineWidth',2);
    p=plot([Result(II,1)+1 J-1],[Result(II,4) Result(II,4)]);
    p=plot([Result(II,1)+1 J-1],[Result(II,3) Result(II,3)]);
    title(sprintf('First derevative of gene number %d',II));
    %%
    subplot(2,2,3)
    hold on;
    p=plot(T,mean(DA));
    set(p,'LineWidth',2,'color','black','linestyle','--');
    plot([Result(II,1)+1,Result(II,1)+1],[min(DA(:)),max(DA(:))],'--','LineWidth',2);
    p=plot([Result(II,1)+1 J-1],[Result(II,4) Result(II,4)]);
    p=plot([Result(II,1)+1 J-1],[Result(II,3) Result(II,3)]);
    title(sprintf('mean First derevative of gene number %d',II));
    %%
    T=1:J-2;
    subplot(2,2,4)
    for i=1:I
        p=plot(T,SE_D(i,:,II));
        set(p,'Color',rand(3,1));
        hold on;
    end
    DA=SE_D(:,:,II);
    p=plot(T,mean(DA));
    set(p,'LineWidth',2,'color','black','linestyle','--');
    plot([Result(II,1)+1,Result(II,1)+1],[min(DA(:)),max(DA(:))],'--','LineWidth',2);
    title(sprintf('Second derevative of gene number %d',II));    
end
figure;
cnames = {'Gene Number','Change Point','Test Value','Upper Limit','Lower Limit'};
DA=[(1:II)' Result];
t = uitable('Data',DA,'ColumnName',cnames,'Position',[0 300 600 600]);

function F_T=detector(data,l)
I=size(data,1); J=size(data,2);
X_bar=mean(data);
X_bar=repmat(X_bar,I,1); Z=data-X_bar; S=(Z'*Z)/(I);
control=0;
control1=0;
t=J-5;
F0=sum((data<=0))/I;
while control==0
    JJ=min(t+l,J);
    %JJ=J;
    F_0=sum(F0(t:JJ));
    Exp=(JJ-t+1)*0.5;
    FF=0;
    for i=t:JJ
        for j=i:JJ
            if i~=j
                SIGMA=[S(i,i) S(i,j); S(j,i) S(j,j)];
                mu=[X_bar(1,i); X_bar(1,j)];
                FF=mvncdf([0;0],mu,SIGMA)+FF;
            end
        end
    end
    STD=sqrt((0.25*(JJ-t+1)*(t-JJ+1)+2*FF)/I);
    Z=(F_0-Exp)/STD;
    if abs(Z)>norminv(1-0.2/2) && control1==0
        t_warning=t;
        control1=1;
    end
    if abs(Z)>norminv(1-0.05/2)
        control=1;
    elseif t==1
        control=1;
    else
        t=t-1;
    end
end
control2=0;
Z=(F0-0.5)/sqrt(0.25/I);
if  t~=t_warning
    for i=t_warning:-1:t+1
        if control2==0
            if abs(Z(i))>1.96
                F_T=i;
                control2=1;
            end
        end
        if control2==0 && i==t+1
            F_T=t;
        end
    end
else
    F_T=t;
end
end

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
    set(p,'LineWidth',2,'color','black','line','--');
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
    set(p,'LineWidth',2,'color','black','line','--');
    plot([Result(II,1)+1,Result(II,1)+1],[min(DA(:)),max(DA(:))],'--','LineWidth',2);
    p=plot([Result(II,1)+1 J-1],[Result(II,4) Result(II,4)]);
    p=plot([Result(II,1)+1 J-1],[Result(II,3) Result(II,3)]);
    title(sprintf('First derevative of gene number %d',II));
    %%
    subplot(2,2,3)
    hold on;
    p=plot(T,mean(DA));
    set(p,'LineWidth',2,'color','black','line','--');
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
    set(p,'LineWidth',2,'color','black','line','--');
    plot([Result(II,1)+1,Result(II,1)+1],[min(DA(:)),max(DA(:))],'--','LineWidth',2);
    title(sprintf('Second derevative of gene number %d',II));    
end
figure;
cnames = {'Gene Number','Change Point','Test Value','Upper Limit','Lower Limit'};
DA=[(1:II)' Result];
t = uitable('Data',DA,'ColumnName',cnames,'Position',[0 300 600 600]);



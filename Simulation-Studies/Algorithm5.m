% Derevative Function
clearvars -global; clearvars; 
clc;
disp(['  ' 'condition' '    ' 'alpha' '     ' 'mean' '        ' 'mse']);
global I J K f fun data; 
for alpha1=[0.1 0.05 .025];
I=30; J=50; K=1; R=zeros(J,J); MU=zeros(1,J); D=zeros(1,J); C=25;
T=linspace(0,1,J); ro=0.5; alpha=2; ro=ro^(1/J); sigma2=.2;
t=linspace(0,1,J); 
%f=@(x)sin(pi*x)-x.^2-.5;
%f=@(x)sin(pi*x)-x.^2+1;
f=@(x)4*pi*cos(4*pi*x)-6*(x-3/4).^2-12.1936;
fun=f(t);
fun(C+1:end)=0;
plot(t,fun);
title('Derevative function');
xlabel('Time');
ylabel('Derevative function');
%%
T1=T;
f1=f;
index=1;
accept=0;
for rep=1:5000
    i=1;
    data=RAND(I, J, alpha, f1, sigma2, ro, C, T1);
    f_plus=data(:,2:end,:); f=data(:,1:end-1,:); D=f_plus-f; T=1:J-1; data=D;
    %s_mean_0=SRD(I, J, T, data);
    X_bar=mean(data);
    X_bar=repmat(X_bar,I,1); Z=data-X_bar; S=(Z'*Z)/(I);
    control=0;
    control1=0;
    t=J-4;
    F0=sum((data<=0))/I;
    while control==0
        JJ=min(t+1,J);
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
        if abs(Z)>norminv(1-0.1/2) && control1==0
            t_warning=t;
            control1=1;
        end
        if abs(Z)>norminv(1-alpha1/2)
            control=1;
        elseif t==1
            control=1;
        else
            t=t-1;
        end
    end
    control2=0;
    Z=(F0-0.5)./sqrt((F0.*(1-F0))/I);
    if  t~=t_warning
        for i=t_warning:-1:t+1
            if control2==0
                if abs(Z(i))>norminv(1-alpha1/2)
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
    RES(rep)=F_T;
end
mm=mean(RES(RES<30))+1;
mse=mean((RES(RES<30)-C));
RES=RES(RES<30);
disp([I alpha1 mm mse mean(RES==(C))]);
end







% Derevative Function
clearvars;
global I J K T X NB B; 
% I is the number of conditions 
% J number of time points
% C is the cut point
I=80; J=25; C=19;
K=1; 
%---------------------------------------------------------%
R=zeros(J,J); MU=zeros(1,J); D=zeros(1,J);
T=linspace(0,1,J); ro=0.5; alpha=2; ro=ro^(1/J); sigma2=.2;
t=0:0.04:1; %t is the time vector is should define based on the value of J
%f=@(x)sin(pi*x)-x.^2-.5;
%f=@(x)sin(pi*x)-x.^2+1;
f=@(x)sin(2*pi*x)-4*(x-3/4).^2+1;
% f is the mean function which is defied as a handel function
fun=f(t);
fun(C+1:end)=0;
plot(t,fun);
title('Derevative function');
xlabel('Time');
ylabel('Derevative function');
%%
for rep=1:100
    for L=1:1
        for i=1:J
            for j=1:J
                R(i,j)=ro^abs(i-j);
            end
            if i<=C
                fun=f(T(i));
            else
                fun=0;
            end
            D(i)=exp(0.5*alpha*fun);
            MU(i)=fun;
        end
        RR=R;
        D1=diag(D); R=sigma2*D1*R*D1;
        data=mvnrnd(MU',R,I);
        if L==1
            DATA=data;
        else
            DATA=DATA+data;
        end
    end
    X=DATA/L;
    data=X;
    Mean_m(:,:)=repmat(mean(data),I,1);
    
    %%
    C=0; back=2;
    S=cov(X);
    while C==0
        SS=S(J-2+1:J,J-2+1:J);
        SS=sum(sum(SS));
        T_sq=mean(sum(X(:,J-back+1:J)'))/sqrt(sum(sum(S))/I);
        disp(T_sq);
        Z=norminv(0.975,0,1);
        if abs(T_sq)>=Z||back==J
            C=1;
        else
            back=back+1;
        end
    end
    CUT(rep)=back+1;
    disp(rep);
end
mean(CUT)

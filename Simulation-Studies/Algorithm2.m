clearvars; clc;
global I J K T X;
I=60; J=25; K=1; R=zeros(J,J); MU=zeros(1,J); D=zeros(1,J); C=19;
T=linspace(0,1,J); ro=0.5; alpha=2; ro=ro^(1/J); sigma2=.2;
t=T;
f=@(x)sin(2*pi*x)-4*(x-3/4).^2+1;
fun=f(t);
fun(C+1:end)=0;
plot(t,fun);
title('Derevative function');
xlabel('Time');
ylabel('Derevative function');
%%
rep=10000;
for M=1:rep
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
    %plot(X(:,:)');
    control=0;
    R=1;
    while control==0
        XX=X(:,24-R:25);%6-R);
        A=XX-repmat(mean(XX),size(XX,1),1);
        S=A'*A/(I-1);
        T_sq=I*(mean(XX))*cov(XX)^-1*(mean(XX))';
        Test=(((I-1)*3)/(J-3))*finv(0.95,3,I-3);
        Result=T_sq<Test;
        if Result==0
            Point(M)=24-R;
            control=1;
        end
        R=R+1;
        if R>20
            control=1;
        end
    end
    R=1;
    control=0;
    while control==0
        XX=X(:,24-R:25);
        [COEFF,SCORE,latent] = princomp(XX);
        TE=XX*COEFF;
        n_alpha=1-(1-.01)^(1/size(XX,2));
        n_alpha=.05/size(TE,2);
        Z=tinv(1-n_alpha/2,I-1);
        TEST=mean(TE)./sqrt(var(TE)/size(TE,1));
        Result=floor(sum(abs(TEST)<Z)/size(TE,2));
        if Result==0
            Point1(M)=24-R;
            control=1;
        end
        R=R+1;
        if R>20
            control=1;
        end
    end
end
disp(mean(Point)+1);
disp(mean(Point1)+1);
hold on;
plot([(mean(Point)+1)/J (mean(Point)+1)/J],[-1 1],'color','red')
plot([(mean(Point1)+1)/J (mean(Point1)+1)/J],[-1 1])
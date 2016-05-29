% Derevative Function
clearvars;
global I J K T X; 
I=15; J=50; K=1; R=zeros(J,J); MU=zeros(1,J); D=zeros(1,J); C=24;
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
for i=1:J
    for j=1:J
        R(i,j)=ro^abs(i-j);
    end
    if i<=C
        fun1=f(T(i));
    else
        fun1=0;
    end
    D(i)=exp(0.5*alpha*fun1);
    MU(i)=fun1;
end
RR=R;
D1=diag(D); R=sigma2*D1*R*D1;
h1=@(z) z;
h2=@(z) z.^2-1;
h3=@(z) z.^3-3*z;
h4=@(z) z.^4-6*z.^2+3;
h5=@(z) z.^5-10*z.^3+15*z;
h6=@(z) z.^6-15*z.^4+45*z.^2-15;
for rep=1:5000
    DATA=[];
    for L=1:1
        data=mvnrnd(MU',R,I);
        if L==1
            DATA=data;
        else
            DATA=DATA+data;
        end
    end
    X=DATA/L;
    data=X;
     f_plus=data(:,2:end,:); f=data(:,1:end-1,:); D=f_plus-f; data=D; X=D;
    SS=cov(X);
    %%
    C=0; back=3;
    S=cov(X);
    Tresh=.6;
    n1=I;
    J=size(D,2);
    Tresh=.6;
    while C==0
        WEI=zeros(1,back);
        WEI(1)=Tresh;
        WEI(2:end)=(1-Tresh)/(back-1);
        WEI=ones(1,back);
        WEI=repmat(WEI,I,1);
        data1=WEI.*data(:,J-back+1:J);
        S=cov(data1);
        y=data(:,J-back+1:J);
        y_b=mean(y)';
        TD=I*(y_b'*y_b)-trace(S);
        p=back;
        TD1=n1*(y_b'*y_b)/trace(S);
        R2=sqrt(p)*(TD1-1);    
        b0=n1*(n1^3+6*n1^2+21*n1+18);
        b1=2*n1*(2*n1^2+6*n1+9);
        b2=2*n1*(3*n1+2);
        b3=n1*(2*n1^2+5*n1+7);
        a1=trace(S)/p;
        a2=n1^2/(p*(n1-1)*(n1+2))*(trace(S^2)-(trace(S))^2/n1);
        a3=n1/((n1-1)*(n1-2)*(n1+2)*(n1+4))*(trace(S^3)/p-3*(n1+2)*(n1-1)*a1*a2-n1*p^2*a1^3);
        a4=1/b0*(trace(S^4)/p-p*b1*a1-p^2*b2*a1^2*a2-p*b3*a2^2-n1*p^3*a1^4);
        c2=.5;
        c3=sqrt(2)*a3/(3*sqrt(a2^3));
        c4=a4/(2*a2^2);
        c6=a3^2/(9*a2^3);
        z1=R2/sqrt(2*a2/a1^2);
        z=TD/sqrt(2*p*a2);
        data(i)=z;
        
        Re(i)=normcdf(z)-normpdf(z).*(1/sqrt(p)*c3*h2(z)+1/p*(c4*h3(z)+c6*h5(z))+1/n1*c2*h1(z));
        Z=norminv(.955);
        Z_a=Z+(1/sqrt(p))*(sqrt(2)*a3/(3*sqrt(a2^3)))*(Z^2-1)+1/p*(a4/(2*a2^2)*Z*(Z^2-3)-2*a3^2/(9*a2^3)*Z*(2*Z^2-5))+1/(2*n1)*Z;
        PP(i)=z<Z_a;  
        %disp([z Z_a Re(i)])
        if abs(z)>=Z_a||back==J
            C=1;
        else
            back=back+1;
        end
    end
    CUT(rep)=back-1;
    disp(rep);
end
disp([mean(CUT) mean(abs(CUT-26))])

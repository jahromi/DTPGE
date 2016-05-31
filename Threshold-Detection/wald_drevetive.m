function [T, beta_0 V_B]=wald_drevetive(data,cut)
X=data(:,cut:end); I=size(X,1); J=size(X,2);
B=ones(J,1); NB=size(B,2);
yy=reshape(X',I*J,1); BB=repmat(B,I,1); X_bar=mean(X);
X_bar=repmat(X_bar,I,1); Z=X-X_bar; S=(Z'*Z)/(I);
for i=1:I 
    block{i}=S;
end
S_b=blkdiag(block{:}); beta_0=(BB'*S_b*BB)^-1*BB'*S_b*yy;
LF=(beta_0-0*abs(beta_0))'; UF=(beta_0+0*abs(beta_0))';
s=diag(S); CC=B*beta_0; N_E=1;
for i=1:J
    for j=1:J
        rho_m(i,j)=S(i,j)/sqrt(s(i)*s(j));
        if i~=j
            Rho(N_E)=rho_m(i,j)^(1/abs(i-j));  
            N_E=N_E+1;
        end
    end
end
rho=mean(real(Rho));
N_E=1;
for i=1:J
    for j=1:J
        Alpha_p(N_E)=S(i,j)/(rho^abs(i-j));
    end
end
alpha_p=mean(Alpha_p);
L=[beta_0 rho alpha_p];
for i=1:J
    for j=1:J
        R(i,j)=rho^(abs(i-j));     
    end
end
V_B=sqrt(alpha_p*(B'*R^-1*B)^-1);
T=abs(beta_0/V_B);
end
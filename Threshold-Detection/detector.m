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

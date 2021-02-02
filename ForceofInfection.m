function [Lambda]=ForceofInfection(x,beta,epsi,Nage,Ncomp,VEI,r,nGD)

eps=1e-10;

%Force of infection
rhoN=zeros(Nage,1);

% rhoN=sum(x(:,1:8),2)+sum(x(:,11:10+nGD+7+(nGD-1)),2);
rhoN=sum(x(:,1:16+3*(nGD-1)),2);

rhoNtot=sum(rhoN(1:Nage));

%Mixing Matrix
%%%Age group
II=eye(Nage);
H=zeros(Nage,Nage);
for a=1:Nage
    H(:,a)=epsi.*II(:,a)+(1-epsi).*(rhoN(a)/(rhoNtot+eps)).*ones(Nage,1);
end

LambdaR=zeros(Nage,1);

for a=1:Nage %Age
    LambdaR(:)=LambdaR(:)+beta(:).*H(:,a).*((x(a,3)+x(a,4)+x(a,6)+(1-VEI)*(1+r)*(x(a,11+2*(nGD-1))+x(a,11+2*(nGD-1)+1)+x(a,11+2*(nGD-1)+3)))/(rhoN(a)+eps));
end
Lambda=LambdaR;

end

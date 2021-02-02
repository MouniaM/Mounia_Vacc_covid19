function [dxL]=risk_Corona19_ChinaVA2(t,x,tf,a0,a1,b1,c1,ML,Beta,tEas,Slope,sigma,alpha,delta,NAge,mu,epsi,nst,eta,fM,fS,fC,fM1,fS1,fC1,nuM,nuS,nuC,nuSID,nuCID,tVacc,DurVacc,AV,gamma,gamma2,VES,VEI,VEP1,VEP2,omega,r,nGD,RateGD,nuRGD)
%%

xt=zeros(NAge,nst);
iiii=1:nst;
for j=1:NAge
    isj= iiii+(j-1)*nst;
    xt(j,iiii)=x(isj);
end 
%%
beta=zeros(NAge,1);

% for j=1:NAge
% %      beta(j)=a0*(1+(a1/(1+exp((t-b1)/c1)))).*sigma(j); %(ZA*exp(-((j-b2)/a2)^2)); %*(1/(1+exp(-a2*(j-b2))));
% beta(j)=ML.*Beta.*sigma(j);
% end

for a=1:NAge
    if (t>=tVacc && t<=tVacc+tEas)
        beta(a)=(ML.*Beta+Slope*((t-1)-tVacc)).*sigma(a);
    elseif (t>tVacc+tEas && t<=tf)
        beta(a)=(ML.*Beta+Slope*((tVacc+tEas)-tVacc)).*sigma(a);
    else
        beta(a)=ML.*Beta.*sigma(a);
    end
end
% beta
lambda=ForceofInfection(xt,beta,epsi,NAge,8,VEI,r,nGD);

%%
% if nGD>1
dx=zeros(NAge,nst);
dx(1,1)=0-(mu+eta(1)+lambda(1)+VaccRate(t,tVacc,DurVacc,NAge,1,AV,gamma,gamma2))*xt(1,1)+RateGD*xt(1,8+nGD+(nGD-1))+nuRGD*xt(1,8+(nGD-1))+RateGD*xt(1,8+2*nGD+6+(nGD-1)); %(S)0
dx(1,2)=lambda(1)*xt(1,1) - (mu+delta+eta(1))*xt(1,2) ;       %(E)

dx(1,3)=fM(1)*delta*xt(1,2) - (mu+nuM+eta(1))*xt(1,3);      %(IM)

dx(1,4)=fS(1)*delta*xt(1,2) - (mu+nuSID+eta(1))*xt(1,4);      %(IS) 
dx(1,5)=nuSID*xt(1,4) - (mu+nuS+eta(1))*xt(1,5);      %(DS)

dx(1,6)=fC(1)*delta*xt(1,2) - (mu+nuCID+eta(1))*xt(1,6);      %(IC)
dx(1,7)=nuCID*xt(1,6) - (mu+nuC+eta(1)+alpha(1))*xt(1,7);      %(DC)
dx(1,8)=nuM*xt(1,3)+nuS*xt(1,5)+nuC*xt(1,7)-(mu+eta(1)+VaccRate(t,tVacc,DurVacc,NAge,1,AV,gamma,gamma2)+nuRGD)*xt(1,8);      %(R_1)
for ii=2:nGD
    dx(1,8+(ii-1))=nuRGD*xt(1,8+(ii-1)-1)-(mu+eta(1)+VaccRate(t,tVacc,DurVacc,NAge,1,AV,gamma,gamma2)+nuRGD)*xt(1,8+(ii-1)); %(R_2:R_nGD)
end
dx(1,8+nGD)=VaccRate(t,tVacc,DurVacc,NAge,1,AV,gamma,gamma2)*xt(1,1)-(mu+eta(1)+(1-VES)*(1+r)*lambda(1)+RateGD)*xt(1,8+nGD); %(SV_1)
for ii=2:nGD
    dx(1,8+nGD+(ii-1))=RateGD*xt(1,8+nGD+(ii-1)-1)-(mu+eta(1)+(1-VES)*(1+r)*lambda(1)+RateGD)*xt(1,8+nGD+(ii-1)); % (SV_2:SV_nGD)
end
dx(1,8+2*nGD)=(1-VES)*(1+r)*lambda(1)*(sum(xt(1,8+nGD:8+2*nGD-1)))-(mu+delta+eta(1))*xt(1,8+2*nGD) ;       %(EV)

dx(1,8+2*nGD+1)=fM1(1)*delta*xt(1,8+2*nGD) - (mu+nuM+eta(1))*xt(1,8+2*nGD+1);      %(IMV)

dx(1,8+2*nGD+2)=fS1(1)*delta*xt(1,8+2*nGD) - (mu+nuSID+eta(1))*xt(1,8+2*nGD+2);      %(ISV) 
dx(1,8+2*nGD+3)=nuSID*xt(1,8+2*nGD+2) - (mu+nuS+eta(1))*xt(1,8+2*nGD+3);      %(DSV)

dx(1,8+2*nGD+4)=fC1(1)*delta*xt(1,8+2*nGD) - (mu+nuCID+eta(1))*xt(1,8+2*nGD+4);      %(ICV)
dx(1,8+2*nGD+5)=nuCID*xt(1,8+2*nGD+4) - (mu+nuC+eta(1)+alpha(1))*xt(1,8+2*nGD+5);      %(DCV)

dx(1,8+2*nGD+6)=VaccRate(t,tVacc,DurVacc,NAge,1,AV,gamma,gamma2)*xt(1,8)+(nuM)*xt(1,8+2*nGD+1)+nuS*xt(1,8+2*nGD+3)+nuC*xt(1,8+2*nGD+5)-(mu+eta(1)+RateGD)*xt(1,8+2*nGD+6);          %(RV_1)
for ii=2:nGD
   dx(1,8+2*nGD+6+(ii-1))=VaccRate(t,tVacc,DurVacc,NAge,1,AV,gamma,gamma2)*xt(1,8+(ii-1))+RateGD*xt(1,8+2*nGD+6+(ii-1)-1)-(mu+eta(1)+RateGD)*xt(1,8+2*nGD+6+(ii-1)); %% (RV_2:RV_nGD)
end
dx(1,16+3*(nGD-1)+1)=delta*xt(1,2);    %%Cumulative incidence non vaccinated
dx(1,16+3*(nGD-1)+2)=alpha(1)*xt(1,7); %%Cumulative deaths  non vaccinated

dx(1,16+3*(nGD-1)+3)=delta*xt(1,10+2*(nGD-1));    %%Cumulative incidence vaccinated
dx(1,16+3*(nGD-1)+4)=alpha(1)*xt(1,15+2*(nGD-1)); %%Cumulative deaths vaccinated

for rst=2:NAge
dx(rst,1)=eta(rst-1)*xt(rst-1,1)-(mu+eta(rst)+lambda(rst)+VaccRate(t,tVacc,DurVacc,NAge,rst,AV,gamma,gamma2))*xt(rst,1)+RateGD*xt(rst,8+2*(nGD-1)+1)+nuRGD*xt(rst,8+(nGD-1))+RateGD*xt(rst,8+2*nGD+6+(nGD-1));             %(S)
dx(rst,2)=eta(rst-1)*xt(rst-1,2)+lambda(rst)*xt(rst,1) - (mu+delta+eta(rst))*xt(rst,2) ;       %(E)

dx(rst,3)=eta(rst-1)*xt(rst-1,3)+fM(rst)*delta*xt(rst,2) - (mu+nuM+eta(rst))*xt(rst,3);      %(IM)

dx(rst,4)=eta(rst-1)*xt(rst-1,4)+fS(rst)*delta*xt(rst,2) - (mu+nuSID+eta(rst))*xt(rst,4);      %(IS) 
dx(rst,5)=eta(rst-1)*xt(rst-1,5)+nuSID*xt(rst,4) - (mu+nuS+eta(rst))*xt(rst,5);      %(DS)

dx(rst,6)=eta(rst-1)*xt(rst-1,6)+fC(rst)*delta*xt(rst,2) - (mu+nuCID+eta(rst))*xt(rst,6);      %(IC)
dx(rst,7)=eta(rst-1)*xt(rst-1,7)+nuCID*xt(rst,6) - (mu+nuC+eta(rst)+alpha(rst))*xt(rst,7);      %(DC)
dx(rst,8)=eta(rst-1)*xt(rst-1,8)+nuM*xt(rst,3)+nuS*xt(rst,5)+nuC*xt(rst,7)-(mu+eta(rst)+VaccRate(t,tVacc,DurVacc,NAge,rst,AV,gamma,gamma2)+nuRGD)*xt(rst,8);      %(R_1)
for ii=2:nGD
    dx(rst,8+(ii-1))=eta(rst-1)*xt(rst-1,8+(ii-1))+nuRGD*xt(rst,8+(ii-1)-1)-(mu+eta(rst)+VaccRate(t,tVacc,DurVacc,NAge,rst,AV,gamma,gamma2)+nuRGD)*xt(rst,8+(ii-1)); %(R_2:R_nGD)
end
dx(rst,8+nGD)=eta(rst-1)*xt(rst-1,8+nGD)+VaccRate(t,tVacc,DurVacc,NAge,rst,AV,gamma,gamma2)*xt(rst,1)-(mu+eta(rst)+(1-VES)*(1+r)*lambda(rst)+RateGD)*xt(rst,8+nGD); %(SV_1)
for ii=2:nGD
    dx(rst,8+nGD+(ii-1))=eta(rst-1)*xt(rst-1,8+nGD+(ii-1))+RateGD*xt(rst,8+nGD+(ii-1)-1)-(mu+eta(rst)+(1-VES)*(1+r)*lambda(rst)+RateGD)*xt(rst,8+nGD+(ii-1)); %%(SV_2:SV_nGD)
end
dx(rst,8+2*nGD)=eta(rst-1)*xt(rst-1,8+2*nGD)+(1-VES)*(1+r)*lambda(rst)*(sum(xt(rst,8+nGD:8+2*nGD-1)))-(mu+delta+eta(rst))*xt(rst,8+2*nGD) ;       %(EV)

dx(rst,8+2*nGD+1)=eta(rst-1)*xt(rst-1,8+2*nGD+1)+fM1(rst)*delta*xt(rst,8+2*nGD) - (mu+nuM/(1-VEP1)+eta(rst))*xt(rst,8+2*nGD+1);      %(IMV)

dx(rst,8+2*nGD+2)=eta(rst-1)*xt(rst-1,8+2*nGD+2)+fS1(rst)*delta*xt(rst,8+2*nGD) - (mu+nuSID+eta(rst))*xt(rst,8+2*nGD+2);      %(ISV) 
dx(rst,8+2*nGD+3)=eta(rst-1)*xt(rst-1,8+2*nGD+3)+nuSID*xt(rst,8+2*nGD+2) - (mu+nuS+eta(rst))*xt(rst,8+2*nGD+3);      %(DSV)

dx(rst,8+2*nGD+4)=eta(rst-1)*xt(rst-1,8+2*nGD+4)+fC1(rst)*delta*xt(rst,8+2*nGD) - (mu+nuCID+eta(rst))*xt(rst,8+2*nGD+4);      %(ICV)
dx(rst,8+2*nGD+5)=eta(rst-1)*xt(rst-1,8+2*nGD+5)+nuCID*xt(rst,8+2*nGD+4) - (mu+nuC+eta(rst)+alpha(rst))*xt(rst,8+2*nGD+5);      %(DCV)

dx(rst,8+2*nGD+6)=eta(rst-1)*xt(rst-1,8+2*nGD+6)+VaccRate(t,tVacc,DurVacc,NAge,rst,AV,gamma,gamma2)*xt(rst,8)+(nuM)*xt(rst,8+2*nGD+1)+nuS*xt(rst,8+2*nGD+3)+nuC*xt(rst,8+2*nGD+5)-(mu+eta(rst)+RateGD)*xt(rst,8+2*nGD+6);          %(RV_1)
for ii=2:nGD
   dx(rst,8+2*nGD+6+(ii-1))=eta(rst-1)*xt(rst-1,8+2*nGD+6+(ii-1))+VaccRate(t,tVacc,DurVacc,NAge,rst,AV,gamma,gamma2)*xt(rst,8+(ii-1))+RateGD*xt(rst,8+2*nGD+6+(ii-1)-1)-(mu+eta(rst)+RateGD)*xt(rst,8+2*nGD+6+(ii-1)); %(RV_2:RV_nGD)
end
dx(rst,16+3*(nGD-1)+1)=delta*xt(rst,2);    %%Cumulative incidence non vaccinated
dx(rst,16+3*(nGD-1)+2)=alpha(rst)*xt(rst,7); %%Cumulative deaths  non vaccinated

dx(rst,16+3*(nGD-1)+3)=delta*xt(rst,10+2*(nGD-1));    %%Cumulative incidence vaccinated
dx(rst,16+3*(nGD-1)+4)=alpha(rst)*xt(rst,15+2*(nGD-1)); %%Cumulative deaths vaccinated
end

dxL=zeros(NAge*nst,1);
iii=1:nst;
for j=1:NAge
    isj=iii+(j-1)*nst;
dxL(isj)= dx(j,iii);
end

end




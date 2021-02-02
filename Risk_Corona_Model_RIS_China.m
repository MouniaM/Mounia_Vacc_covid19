clc
clear all
%% GAmma Distribution
DurImm=1*365;
DurNImm=1*365;
nGD=15;
DurGD=DurImm/nGD;
DurNGD=DurNImm/nGD;
RateGD=(1/DurGD);
nuRGD=1/DurNGD;
%% define the number of stages and INFECTION groups.
nst=16+3*(nGD-1)+4; % NUMBER OF STATES AND STAGES S, E, I, R
NAge=9; %number of age groupss
%% US Population
DemoChina='CopyofUSADemo.xlsx';
PopSize= xlsread(DemoChina);

InciChina='DataCOVIDChina.xlsx';
IncidenceChina=xlsread(InciChina);
N0=zeros(NAge,1);
% % % %China Population
 for a=2:2:16
     N0(round(a/2))=PopSize(a,3)+PopSize(a-1,3);  
 end
N0(9)=PopSize(17,3)+PopSize(18,3)+PopSize(19,3)+PopSize(20,3);
mu=1/(79.10*365);   %US     % 1/life expectancy

%% China Population
% DemoChina='ChineseDemo.xlsx';
% PopSize= xlsread(DemoChina);
% 
% InciChina='DataCOVIDChina.xlsx';
% IncidenceChina= xlsread(InciChina);
% N0=zeros(NAge,1);
% % % % %China Population
%  for a=2:2:16
%      N0(round(a/2))=PopSize(a,3)+PopSize(a-1,3);  
%  end
%  N0(9)=PopSize(17,3)+PopSize(18,3)+PopSize(19,3)+PopSize(20,3);
% mu=1/(76.47*365);  % China     % 1/life expectancy
%% Qatar Population
% % % 
% DemoQatar='QatarDemoH.xlsx';
% PopSize= xlsread(DemoQatar);
% InciChina='DataCOVIDChina.xlsx';
% IncidenceChina= xlsread(InciChina);
% % N0=zeros(NAge,1);
% % for a=2:2:16
% %     N0(round(a/2))=PopSize(a,1)+PopSize(a-1,1);  
% % end
% % N0(9)=PopSize(17,1);
% N0(1:9)=PopSize(1:9,3);
% 
% mu=1/(80.7*365);  % Qatar     % 1/life expectancy
% 

%% Vaccinaation parameters
tVacc=626;% US
DurVacc=365;
tVacc80=tVacc+DurVacc %US
% tVacc=75.5; % China
% tVacc80=tVacc+91; %China
gamma=0.0038;
gamma2=0;
VES=0;
VEI=0;
VEP1=0;
VEP2=0;0.95;
omega=0;1/365;
r=0
AV=1;
dir3=['US']

% filename=['R01p2NoeffAG' num2str(AV)]
% filename=['R01p2VES07AG' num2str(AV)]
   filename=['R01p2Noeff1']

%% Distribution of cases according to severity in the non-vaccinated population
% Identical to the excel sheet
fSAG4=0.011;
RRS=[0.14 0.14 0.49 1 1.21 2.32 4.49 7.76 27.63];
fS=fSAG4.*RRS;

fCAG4=0.002;
RRC=[0.21 0.21 0.33 1 1.83 4.67 10.58 13.61 8.67];
fC=fCAG4.*RRC;

fM=1-fS-fC;

Alpha4=0.0002;
RRAlpha=[0.1 0.1 0.4 1 3 10 45 120 505];
alpha=Alpha4.*RRAlpha;
%%%%%

%% Distribution of cases according to severity in the vaccinated population
for a=1:NAge
fM1(a)=((1-VEP2)+VEP2./fM(a)).*fM(a);
fS1(a)=(1-VEP2).*fS(a);
fC1(a)=(1-VEP2).*fC(a);
end
% sumfractions=fM1+fS1+fC1;

delta=1/(3.69);  % infection rate (from exposed)
nuM=1/(3.48); %
nuS=1/(3.5*7+7-3.48);
nuC=1/(5.5*7+7-3.48);
nuSID=1/(3.48);
nuCID=1/(3.48);
nuRT=0;%1/(0.1*7);
fRT=0;

[M,I]=max(IncidenceChina(:,1));

%Good fitts 
load('EstimatedParam.mat','EstimatedParam')
zz=EstimatedParam;
tinit=zz(1);
a0=zz(2);
a1=zz(3);
b1=I+tinit-zz(17); 
c1=zz(4);
Beta=1.06;
%%ML=1.1956;%% R0=2.5
%  ML=0.4075;%% R0=1.5 
%  ML=0.815;%% R0=3
% ML=1.2;%% R0=4
% ML=0.353;  %% R0=1.3
 ML=0.3277; %%R0=1.2 US
Slope=0.00442;
%  Slope=0.0024;
% ML=0.277; %R0=1;
% Slope=0.0047;% R0 4 china
 tEas=182.5;

sigma=ones(NAge,1);
epsi=zz(15);
Anch=zz(16);



IncidenceDelay=zz(18);
MortalityDelay=zz(19);
%%
%%%% first Initial conditions
initprev=1/NAge; %initialnumber of people with corona
x0L=zeros(NAge,nst);
for a=1:NAge
x0L(a,1)=N0(a)-initprev;
x0L(a,2)=0+initprev;
end
x0=zeros(nst*NAge,1);
ii=1:nst;
for j=1:NAge
    isj= ii+(j-1)*nst;
    x0(isj)=x0L(j,ii);
end

%%%%%%%%%%%Vaccination parameters%%%%%%%%%%%%%%%%%
%% Time scale (days)
t0=0;         % Start time
tf=4*365;         % Stop time (Based on Data Points)
dt=0.5;          % Time interval
tspan=t0:dt:tf;  % Timespan to use in the ode45
eta=(ones(NAge,1)).*(1/(10*365));
eta(NAge)=0;


%% Deterministic model
options = odeset('NonNegative',1:NAge*nst-4);
[T,x]=ode45(@(t,x)risk_Corona19_ChinaVA2(t,x,tf,a0,a1,b1,c1,ML,Beta,tEas,Slope,sigma,alpha,delta,NAge,mu,epsi,nst,eta,fM,fS,fC,fM1,fS1,fC1,nuM,nuS,nuC,nuSID,nuCID,tVacc,DurVacc,AV,gamma,gamma2,VES,VEI,VEP1,VEP2,omega,r,nGD,RateGD,nuRGD),tspan,x0,options);
TT=length(T);
size(x)
%%%% Compartments by Age%%%%%%%%%%%%
%%%%Non vaccinated%%%%%%%%%

SusA=x(:,1:nst:nst*NAge); 
LatentA=x(:,2:nst:nst*NAge);
% 
InfectedIMA=x(:,3:nst:nst*NAge);

InfectedISA=x(:,4:nst:nst*NAge);
InfectedDSA=x(:,5:nst:nst*NAge);

InfectedICA=x(:,6:nst:nst*NAge);
InfectedDCA=x(:,7:nst:nst*NAge);
for jj=1:nGD
RecoveredAGD(:,(jj-1)*NAge+1:jj*NAge)=x(:,8+jj-1:nst:nst*NAge); 
end
for tt=1:length(tspan)
    for a=1:NAge
     RecoveredA(tt,a)=sum(RecoveredAGD(tt,a:NAge:nGD*NAge)); 
    end
end
%%%%Vaccinated%%%%%%%%%
for jj=1:nGD
SusVaccAGD(:,(jj-1)*NAge+1:jj*NAge)=x(:,8+nGD+jj-1:nst:nst*NAge); 
end
for tt=1:length(tspan)
    for a=1:NAge
      SusVaccA(tt,a)=sum(SusVaccAGD(tt,a:NAge:nGD*NAge)); 
    end
end

LatentVaccA=x(:,8+2*nGD:nst:nst*NAge);

InfectedIMVaccA=x(:,8+2*nGD+1:nst:nst*NAge);

InfectedISVaccA=x(:,8+2*nGD+2:nst:nst*NAge);
InfectedDSVaccA=x(:,8+2*nGD+3:nst:nst*NAge);

InfectedICVaccA=x(:,8+2*nGD+4:nst:nst*NAge);
InfectedDCVaccA=x(:,8+2*nGD+5:nst:nst*NAge);

for jj=1:nGD
RecoveredVaccAGD(:,(jj-1)*NAge+1:jj*NAge)=x(:,16+2*(nGD-1)+jj-1:nst:nst*NAge); 
end
for tt=1:length(tspan)
    for a=1:NAge
      RecoveredVaccA(tt,a)=sum(RecoveredVaccAGD(tt,a:NAge:nGD*NAge)); 
    end
end
% 
CumulativeIncidenceA=x(:,16+3*(nGD-1)+1:nst:nst*NAge);

CumulativeDeathsA=x(:,16+3*(nGD-1)+2:nst:nst*NAge);

CumulativeIncidenceVaccA=x(:,16+3*(nGD-1)+3:nst:nst*NAge);

CumulativeDeathsVaccA=x(:,16+3*(nGD-1)+4:nst:nst*NAge);

TotalA=SusA+LatentA+InfectedIMA+(InfectedISA+InfectedDSA)+(InfectedICA+InfectedDCA)+RecoveredA;
TotalVaccA=SusVaccA+LatentVaccA+InfectedIMVaccA+(InfectedISVaccA+InfectedDSVaccA)+(InfectedICVaccA+InfectedDCVaccA)+RecoveredVaccA;
PropVA=TotalVaccA./(TotalA+TotalVaccA);
PropA=(TotalA+TotalVaccA)./(sum(TotalA,2)+sum(TotalVaccA,2));

% figure;
% hold on
% for a=1:NAge
%     plot(tspan,PropVA(:,a).*100);
% end
% PropVA(2*(tVacc+365)+1,AV).*100;

Sus=sum(x(:,1:nst:nst*NAge),2)'; 
Latent=sum(x(:,2:nst:nst*NAge),2)';
InfectedIM=sum(x(:,3:nst:nst*NAge),2)';
InfectedIS=sum(x(:,4:nst:nst*NAge),2)';
InfectedDS=sum(x(:,5:nst:nst*NAge),2)';
InfectedIC=sum(x(:,6:nst:nst*NAge),2)';
InfectedDC=sum(x(:,7:nst:nst*NAge),2)';
Recovered=sum(RecoveredA,2).';
%%%%Vaccinated%%%%%%%%%
SusVacc=sum(SusVaccA,2)'; 
LatentVacc=sum(LatentVaccA,2)';
InfectedIMVacc=sum(InfectedIMVaccA,2)';
InfectedISVacc=sum(InfectedISVaccA,2)';
InfectedDSVacc=sum(InfectedDSVaccA,2)';
InfectedICVacc=sum(InfectedICVaccA,2)';
InfectedDCVacc=sum(InfectedDCVaccA,2)';
RecoveredVacc=sum(RecoveredVaccA,2)';

Total=Sus+Latent+InfectedIM+(InfectedIS+InfectedDS)+(InfectedIC+InfectedDC)+Recovered;
TotalVacc=SusVacc+LatentVacc+InfectedIMVacc+(InfectedISVacc+InfectedDSVacc)+(InfectedICVacc+InfectedDCVacc)+RecoveredVacc;
PropV=TotalVacc./(Total+TotalVacc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% second ODEs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PropV(2*tVacc80+1).*100





%% without vaccination
for t=1:TT
    for a=1:NAge
        IncidenceInfectedAge(t,a)=delta*(LatentA(t,a)+LatentVaccA(t,a));    %%%% Number of new infected in each age group   
        IncMortalityCasesA(t,a)=alpha(a).*(InfectedDCA(t,a)+InfectedDCVaccA(t,a));
        IncRateA(t,a)=(delta*(LatentA(t,a)+LatentVaccA(t,a)))./(SusA(t,a)+SusVaccA(t,a)); %%% incidence rate in each age group
        AttackRateA(t,a)=10000*(CumulativeIncidenceA(t,a)+CumulativeIncidenceVaccA(t,a))./(TotalA(t,a)+TotalVaccA(t,a)); %%% attack rate in each age group
      
        PrevA(t,a)=(InfectedIMA(t,a)+InfectedISA(t,a)+InfectedICA(t,a)+InfectedIMVaccA(t,a)+InfectedISVaccA(t,a)+InfectedICVaccA(t,a))./(TotalA(t,a)+TotalVaccA(t,a));
        InfA(t,a)=InfectedIMA(t,a)+InfectedISA(t,a)+InfectedICA(t,a)+InfectedIMVaccA(t,a)+InfectedISVaccA(t,a)+InfectedICVaccA(t,a);
        
                  
        CumMildIncidenceA(t,a)=dt*(fM(a)*delta*sum(Latent(1:t))+fM1(a)*delta*sum(LatentVacc(1:t))); %%% cumulative incident mild infections
        CumSevereIncidenceA(t,a)=dt*(fS(a)*delta*sum(Latent(1:t))+fS1(a)*delta*sum(LatentVacc(1:t))); %%% cumulative incident severe infections
        CumCriticalIncidenceA(t,a)=dt*(fC(a)*delta*sum(Latent(1:t))+fC1(a)*delta*sum(Latent(1:t)));%%% cumulative incident critical infcetions

    end
    AttackRate(t)=10000*(sum(CumulativeIncidenceA(t,:))+sum(CumulativeIncidenceVaccA(t,:)))./(Total(t)+TotalVacc(t));
    IncRate(t)=sum((delta*(LatentA(t,:)+LatentVaccA(t,:))))./sum(SusA(t,:)+SusVaccA(t,:));
    IncidenceInfectedDS(t)=nuSID*(InfectedIS(t)+InfectedISVacc(t));%%%%  incident severe diseases
    IncidenceInfectedDC(t)=nuCID*(InfectedIC(t)+InfectedICVacc(t));%%%%  incident critical diseases

    CumInfectedDS(t)=dt*nuSID*sum(InfectedISVacc(1:t));%%% cumulative incident severe diseases 
    CumInfectedDC(t)=dt*nuCID*sum(InfectedICVacc(1:t));%%%cumulative incident critical diseases
    
    
    IncidenceInfected(t)=sum(IncidenceInfectedAge(t,:)); %% number of new infceted
end
%%
%%% MainPaper and appendix results%%%%%%%
for t=1:TT
    for a=1:NAge
              
          IncA(t,a)=(delta*(LatentA(t,a)+LatentVaccA(t,a)));
          IncSevereDiseasesA(t,a)=nuSID*(InfectedISA(t,a)+InfectedISVaccA(t,a));
          IncCriticalDiseasesA(t,a)=nuCID*(InfectedICA(t,a)+InfectedICVaccA(t,a));
          IncMortalityCasesA(t,a)=alpha(a).*(InfectedDCA(t,a)+InfectedDCVaccA(t,a));

          CumIncInfA(t,a)=CumulativeIncidenceA(t,a)+CumulativeIncidenceVaccA(t,a); %%% cumulative numer of incident infections in each age group
          CumSevereDiseaseIncA(t,a)=dt*sum(IncSevereDiseasesA(1:t,a)); %%% cumulative incident severe Diseases
          CumCriticalDiseaseIncA(t,a)=dt*sum(IncCriticalDiseasesA(1:t,a));%%% cumulative incident critical Diseases
          CumIncMortalityCasesTA(t,a)=CumulativeDeathsA(t,a)+CumulativeDeathsVaccA(t,a); %%%%Cumulative deaths in each age group


    end
          Inc(t)=sum(IncA(t,:));
         Prev(t)=(sum(InfectedIMA(t,:))+sum(InfectedISA(t,:))+sum(InfectedICA(t,:))+sum(InfectedIMVaccA(t,:))+sum(InfectedISVaccA(t,:))+sum(InfectedICVaccA(t,:)))./(sum(TotalA(t,:))+sum(TotalVaccA(t,:)));

          IncSevereDiseases(t)=sum(IncSevereDiseasesA(t,:));
          IncCriticalDiseases(t)=sum(IncCriticalDiseasesA(t,:));
          IncMortalityCases(t)=sum(IncMortalityCasesA(t,:));
          
          CumInc(t)=sum(CumIncInfA(t,:));
          CumSevereDisInc(t)=sum(CumSevereDiseaseIncA(t,:));
          CumCriticalDisInc(t)=sum(CumCriticalDiseaseIncA(t,:));
          CumIncMortality(t)=sum(CumIncMortalityCasesTA(t,:)); %%% cumulative deaths in the total population
end


%% Repoductive number with mixing
% rhoN=zeros(NAge,1);
% 
% rhoN=sum(x(:,1:8),2)+sum(x(:,11:18),2);
% rhoNtot=sum(rhoN(1:NAge));

%Mixing Matrix
%%%Age group
t00=0
II=eye(NAge);
for ttt=1:length(tspan)
H=zeros(NAge,NAge);
rhoN=TotalA(ttt,:);
rhoNtot=Total(ttt);
for a=1:NAge
    if (tspan(ttt)>=tVacc && tspan(ttt)<=tVacc+tEas)
        beta(a)=(ML.*Beta+Slope*((tspan(ttt)-1)-tVacc)).*sigma(a);
    elseif (tspan(ttt)>=tVacc+tEas)
        beta(a)=(ML.*Beta+Slope*((tVacc+tEas)-tVacc)).*sigma(a);
    else
        beta(a)=ML.*Beta.*sigma(a);
    end
end
for a=1:NAge
    H(:,a)=epsi.*II(:,a)+(1-epsi).*(rhoN(a)/(rhoNtot+eps)).*ones(NAge,1);
end

for a=1:NAge
    for b=1:NAge
        HH(a,b)=H(a,b).*((fM(b).*(delta)./((delta+mu+eta(b)).*(nuM+mu+eta(b))))...
        +(fS(b).*((delta)./((delta+mu+eta(b)).*(nuSID+mu+eta(b)))))...
        +(fC(b).*((delta)./((delta+mu+eta(b)).*(nuCID+mu+eta(b))))));
%     end
%     HHH(a)=((fM(a).*(beta(a).*delta)./((delta+mu+eta(a)).*(nuM+mu+eta(a))))...
%         +(fS(a).*((beta(a).*delta)./((delta+mu+eta(a)).*(nuSID+mu+eta(a)))))...
%         +(fC(a).*((beta(a).*delta)./((delta+mu+eta(a)).*(nuCID+mu+eta(a))))));
    HHH(a)=sum(HH(a,:));
    VA(a)=beta(a).*HHH(a);
%     VA(a)=HHH(a);
    R0Amixing(a)=PropA(2*t00+1,a).*VA(a);
    end
        R0At(ttt,a)=PropA(ttt,a).*VA(a);

end
R0mixing(ttt)=sum(R0Amixing,2);
R0mixing1(ttt)=sum(R0At(ttt,:));

end
 R0mixing(1)
 R0mixing(2*tVacc+1)
 tVacc+tEas
 R0mixing(2*(tVacc+tEas)+1)
figure;
plot(tspan,R0mixing,'LineWidth',2);

 %% Cumulative number of vaccinations
 for t=1:length(tspan)
    if (t<=2*tVacc+1)
        for a=1:NAge
             
         VaccCovVES05A(t,a)=0.*(SusA(t,a)+RecoveredA(t,a));
        end
         elseif (t>=2*(tVacc)+1 && t<=2*(tVacc+DurVacc)+1)
             for a=1:NAge
             if (1<=a<=AV)
         VaccCovVES05A(t,a)=gamma.*(SusA(t,a)+RecoveredA(t,a));
         else
           
         VaccCovVES05A(t,a)=0.*(SusA(t,a)+RecoveredA(t,a));
             end
             end
    else
        for a=1:NAge
        
         VaccCovVES05A(t,a)=0.*(SusA(t,a)+RecoveredA(t,a));
        end
    end

    
  VaccCovVES05(t)=sum(VaccCovVES05A(t,:));
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Cumulative numbner of vaccinations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CumNumVaccVES05=zeros(1,length(tspan));
CumNumVaccVES05=cumtrapz(tspan,VaccCovVES05);
% NbVacc=CumNumVaccVES05(2*(tVacc+14)+1)
NbVacc=CumNumVaccVES05(2*(tVacc80)+1);
PropVCov=CumNumVaccVES05./(Total+TotalVacc);
PV80=PropVCov(2*(tVacc80)+1).*100
% for jj=1:12
%     PVM(jj)=PropVCov(2*(tVacc+jj*(30))+1).*100;
% end
% figure;
% plot(1:12,PVM);

figure(4);
hold on;
plot(tspan,PropV.*100,'LineWidth',2);
plot(tspan,PropVCov.*100,'LineWidth',2);
% plot(tspan,CumNumVaccVES05,'LineWidth',2);

%%
  cd ResultsVES95
  cd (dir3)
  save(filename);
  
txlim=tf;
figure(6);
subplot(3,2,1)
hold on
plot(tspan,AttackRate./100,'LineWidth',2);
xlabel('Days')
ylabel('Attack rate (%)')
xlim([0 txlim])
subplot(3,2,2)
hold on
plot(tspan,IncRate,'LineWidth',2);
xlabel('Days')
ylabel('Incidence rate')
xlim([0 txlim])
subplot(3,2,3)
hold on
plot(tspan,IncSevereDiseases,'LineWidth',2);
xlabel('Days')
ylabel('Incidence of severe diseases')
xlim([0 txlim])
subplot(3,2,4)
hold on
plot(tspan,IncCriticalDiseases,'LineWidth',2);
xlabel('Days')
ylabel('Incidence of critical diseases')
legend('Susceptible','Recovered')
xlim([0 txlim])
subplot(3,2,5)
hold on
plot(tspan,IncMortalityCases,'LineWidth',2);
xlabel('Days')
ylabel('Incidence of deaths')
xlim([0 txlim])
subplot(3,2,6)
hold on
plot(tspan,Inc,'LineWidth',2);
xlabel('Days')
ylabel('Incidence')
xlim([0 txlim])

function rates = HBV2021_Eqs_report(time,state,params)
%rates = HBV2018_Eqs(time, state, params) returns the rate of change of
%each model state
%time
%params.nA = params.params.nA; %number age groups.
%time
% z = find(time<params.tvac_cover,1)-1;
%time
if time==0
t = 1;
%tt=0;
else
t = ceil(time/12);
%tt=t;
end


%% Setup treatment uptake over time - functions defined at end of this script:
% t = the year, eg t=71 refers to 2021. 
if (t<50) % no treatment before year 2000
params.tau1= params.T1(:,t); 
params.tau2= params.T2(:,t); 
params.tau3= params.T3(:,t); 
params.tau4= params.T4(:,t); 
params.tau5= params.T5(:,t); 
params.tau6= params.T6(:,t); 
params.tau7= params.T7(:,t); 
params.tau8= params.T8(:,t);

elseif (t>49 && t<=62) 
params.tau1= treat(t).*params.T1(:,t);
params.tau2= treat(t).*params.T2(:,t);
params.tau3= treat(t).*params.T3(:,t);
params.tau4= treat(t).*params.T4(:,t);
params.tau5= treat(t).*params.T5(:,t);
params.tau6= treat(t).*params.T6(:,t);
params.tau7= treat(t).*params.T7(:,t); 
params.tau8= treat(t).*params.T8(:,t);

elseif (t>62 && t<=71) 

params.tau1= treat2021(t).*params.T1(:,t); 
params.tau2= treat2021(t).*params.T2(:,t); 
params.tau3= treat2021(t).*params.T3(:,t); 
params.tau4= treat2021(t).*params.T4(:,t); 
params.tau5= treat2021(t).*params.T5(:,t); 
params.tau6= treat2021(t).*params.T6(:,t); 
params.tau7= treat2021(t).*params.T7(:,t); 
params.tau8= treat2021(t).*params.T8(:,t);

else 

params.tau1= treat_proj(t).*params.T1(:,t); 
params.tau2= treat_proj(t).*params.T2(:,t); 
params.tau3= treat_proj(t).*params.T3(:,t); 
params.tau4= treat_proj(t).*params.T4(:,t); 
params.tau5= treat_proj(t).*params.T5(:,t); 
params.tau6= treat_proj(t).*params.T6(:,t); 
params.tau7= treat_proj(t).*params.T7(:,t); 
params.tau8= treat_proj(t).*params.T8(:,t);
    
end

%%%%


% Rate of going off treatment
params.Eta1= params.Eta1(:,t);
params.Eta2= params.Eta2(:,t);
params.Eta3= params.Eta3(:,t);
params.Eta4= params.Eta4(:,t);
params.Eta5= params.Eta5(:,t);
params.Eta6= params.Eta6(:,t);
params.Eta7= params.Eta7(:,t);
params.Eta8= params.Eta8(:,t);

% Births:
births = (params.births(t)./12);

% Mortality
Mu = params.mu(:,t); %age by 1 vector

% Vacc coverage
Nu = params.nu(:,t); %age by 1 vector


% This allows us to vary the estimated prevalence of source countries when
% we run Latin hypercube sampling. Only for migrants from source countries
% prior to 1990:
% Countries categorised into low, intermediate and high prevalence
params.prev_COB_1951to1974 =[params.low_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.int_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.int_prev,	params.low_prev,	params.high_prev,	params.int_prev,	params.low_prev,	params.low_prev,	params.int_prev, params.int_prev,	0]';
params.prev_COB_1975to1990 = [params.int_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.high_prev,	params.high_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.int_prev,    params.int_prev,	params.int_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.high_prev,	params.int_prev,	params.high_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.high_prev,	params.low_prev,	params.low_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.high_prev,	params.int_prev,	params.high_prev,	params.int_prev,	params.high_prev,	params.int_prev,	params.int_prev,	params.int_prev,	params.high_prev,	params.int_prev,	params.int_prev,	params.high_prev,	params.int_prev,	params.int_prev,	params.low_prev,	params.low_prev,	params.low_prev,	params.high_prev]';

%% Run the migration input scripts to estimate number of migrants entering the population
% in each phase of CHB:
[params.MIG_S, params.MIG_I, params.MIG_I1, params.MIG_I2, params.MIG_I3, params.MIG_R] = MigrationEstimates_1951to1990(params); % 18 x 40 (1951 - 1990)
[params.mig_S, params.mig_I, params.mig_I1, params.mig_I2, params.mig_I3, params.mig_R] = ABS_MigrationEstimates_1991to2003(params); %18 x 13 (1991 - 2003)
[params.mig_Sa, params.mig_Ia, params.mig_I1a, params.mig_I2a, params.mig_I3a, params.mig_Ra] = ABS_MigrationEstimates_2004onwards(params); %18 x 47 (2004 - 2050) 

%% Change inputs to be per month (instead of per year):
params.MIG_S = [params.MIG_S, params.mig_S, params.mig_Sa]./12; %age by 100 matrix (phay o day thi sao nhi?)'
params.MIG_I = [params.MIG_I, params.mig_I, params.mig_Ia]./12; 
params.MIG_I1 = [params.MIG_I1, params.mig_I1, params.mig_I1a]./12;
params.MIG_I2 = [params.MIG_I2, params.mig_I2, params.mig_I2a]./12;
params.MIG_I3 = [params.MIG_I3, params.mig_I3, params.mig_I3a]./12;
params.MIG_R = [params.MIG_R, params.mig_R, params.mig_Ra]./12;
%params.Mig_all_chron = params.MIG_I + params.MIG_I1 + params.MIG_I2 + params.MIG_I3 + params.MIG_C1 + params.MIG_C2 + params.MIG_C3;


params.Mig_all_chron = params.MIG_I + params.MIG_I1 + params.MIG_I2 + params.MIG_I3;
params.All_mig = params.Mig_all_chron + params.MIG_S + params.MIG_R;

%% Estimate number of migrants in each phase with cirrhosis:
for i=1:18
    %i
    if i<5
    params.MIG_C1(i,:) = params.MIG_I1(i,:)*params.dis_phase_cirr1;
    params.MIG_C2(i,:) = params.MIG_I2(i,:)*params.dis_phase_cirr1;
    params.MIG_C3(i,:) = params.MIG_I3(i,:)*params.dis_phase_cirr1;
    %params.dis_phase_cirr(1)
    elseif  i>4 && i<9
    params.MIG_C1(i,:) = params.MIG_I1(i,:)*params.dis_phase_cirr2;
    params.MIG_C2(i,:) = params.MIG_I2(i,:)*params.dis_phase_cirr2;
    params.MIG_C3(i,:) = params.MIG_I3(i,:)*params.dis_phase_cirr2;
    %params.dis_phase_cirr(2)
    else 
    params.MIG_C1(i,:) = params.MIG_I1(i,:)*params.dis_phase_cirr3;
    params.MIG_C2(i,:) = params.MIG_I2(i,:)*params.dis_phase_cirr3;
    params.MIG_C3(i,:) = params.MIG_I3(i,:)*params.dis_phase_cirr3; 
    %params.dis_phase_cirr(3)
    end
end

% Then adjust numbers entering I1, I2 & I3 (= 1-prop entering C1,C2, C3 etc)
for i=1:18
    %i
    if i<5
    params.MIG_I1(i,:) = params.MIG_I1(i,:)*(1-params.dis_phase_cirr1);
    params.MIG_I2(i,:) = params.MIG_I2(i,:)*(1-params.dis_phase_cirr1);
    params.MIG_I3(i,:) = params.MIG_I3(i,:)*(1-params.dis_phase_cirr1);
    %params.dis_phase_cirr(1)
    elseif  i>4 && i<9
    params.MIG_I1(i,:) = params.MIG_I1(i,:)*(1-params.dis_phase_cirr2);
    params.MIG_I2(i,:) = params.MIG_I2(i,:)*(1-params.dis_phase_cirr2);
    params.MIG_I3(i,:) = params.MIG_I3(i,:)*(1-params.dis_phase_cirr2);
    %params.dis_phase_cirr(2)
    else 
    params.MIG_I1(i,:) = params.MIG_I1(i,:)*(1-params.dis_phase_cirr3);
    params.MIG_I2(i,:) = params.MIG_I2(i,:)*(1-params.dis_phase_cirr3);
    params.MIG_I3(i,:) = params.MIG_I3(i,:)*(1-params.dis_phase_cirr3); 
    %params.dis_phase_cirr(3)
    end
end
%%%%
% Collate total migration numbers by phase:
M_S = params.MIG_S(:,t); %age by 100
M_I = params.MIG_I(:,t); 
M_I1 = params.MIG_I1(:,t);
M_I2 = params.MIG_I2(:,t);
M_I3 = params.MIG_I3(:,t);
M_C1 = params.MIG_C1(:,t);
M_C2 = params.MIG_C2(:,t);
M_C3 = params.MIG_C3(:,t);
M_R = params.MIG_R(:,t);

%Define Health States (variables that are modelled):
S=state(1:params.nA); % Susceptible
V=state(params.nA+1:2*params.nA); % Vaccinated/Immune
A=state(2*params.nA+1:3*params.nA); % Acute Infection
I=state(3*params.nA+1:4*params.nA); % Immune Tolerant 
I1=state(4*params.nA+1:5*params.nA); % Immune Clearance
I2=state(5*params.nA+1:6*params.nA); % Immune Control
I3=state(6*params.nA+1:7*params.nA); % Immune Escape
C1=state(7*params.nA+1:8*params.nA); % Cirrhotic Immune Clearance
C2=state(8*params.nA+1:9*params.nA); % Cirrhotic Immune Control
C3=state(9*params.nA+1:10*params.nA); % Cirrhotic Immune Escape
I1T=state(10*params.nA+1:11*params.nA); % Immune Clearance on Treatment
I2T=state(11*params.nA+1:12*params.nA); % Immune Control on Treatment
I3T=state(12*params.nA+1:13*params.nA); % Immune Escape on Treatment
C1T=state(13*params.nA+1:14*params.nA); % Cirrhotic Immune Clearance on Treatment
C2T=state(14*params.nA+1:15*params.nA); % Cirrhotic Immune Control on Treatment
C3T=state(15*params.nA+1:16*params.nA); % Cirrhotic Immune Escape on Treatment
H=state(16*params.nA+1:17*params.nA);  % Hepatocellular Carcinoma
HT=state(17*params.nA+1:18*params.nA); % Hepatocellular Carcinoma on Treatment
D=state(18*params.nA+1:19*params.nA);  % Decompensated Cirrhosis
DT=state(19*params.nA+1:20*params.nA); % Decompensated Cirrhosis on Treatment
R=state(20*params.nA+1:21*params.nA);  % Recovered/Resolved HBV infection

%TotalPop = S+V+A+I+I1+I2+I3+C1+C2+C3+I1T+I2T+I3T+C1T+C2T+C3T+H+HT+D+DT+R;

%% Force of Infection %%       
% FoI point est based on number of incident notifications around the year 2000 by age
%  this will change for state/territory. Baseline = NNDSS multiplied by 2.
% BC foi pt = [6 6 80 20]; ages 0-4/5-14/15-44/45+
% KM foi pt (from Vic Surv report 2004 - same source as BC):
foi_pt_allages = (params.foi_scale*[3 2 1 2 11 13 7 5 5 1 2 2 1 1 1 1 0.5 0.5]*5)/(1000000*12); %infs/million/month

% Assuming chronic infection stages are 0.16 as infectious as acute stage
eps = 0.16*ones(10,1); eps(1) = 1;
%Total infectious individuals in each age group:
inf = [A I I1 I2 I3 C1 C2 C3 H D]*eps; %gives infectious indivs in ea age group/ assumes those on treatment have negligible trans.

% Run model with FoI = foi_pt_allages' (static foi) and extract number of
% infectious units in the year 2000 (peak of notifs) x10.
%FoI = foi_pt_allages';
% inf_ages_2000= round([A(50,:)' I(50,:)' I1(50,:)' I2(50,:)' I3(50,:)' C1(50,:)' C2(50,:)' C3(50,:)' H(50,:)' D(50,:)']*eps)

pt_inf_2000 = params.inf_units_agegroups_2000; %derived from previous modelled output
% FoI after run using static to 2000:
%FoI = (foi_pt_allages'.*(pt_inf_2000./sum(pt_inf_2000)))./pt_inf_2000;
rel_contrib = pt_inf_2000./sum(pt_inf_2000); %ages by 1 vector
pt_est = (foi_pt_allages'.*rel_contrib)./pt_inf_2000; 

pij_base = ones(18)*diag(pt_est); 

% pij_1 accounts for vert trans and IDU/Sexual contact better
% 0-4s only acquire inf from child bearing age
pij_base_half = pij_base*0.5;
%pij_1 = pij_base_half;
pij_1 = zeros(18);

for ag =5:18
pij_1(ag,ag) = 8*pij_base(ag,ag);
end
for ag1=2:4
    pij_1(ag1,ag1) = 8*pij_base(ag1,ag1);
end
pij_1 = 6*pij_1;
%pij_1(1,[1:4,9:18]) = 0; %only vert trans to 0-4s possible.
pij_1(1,5:8) = pij_base_half(1,5:8); %*2 or mult by 0.5
pij_1(:,1)=0; %0-4s don't transmit infection to anyone else.

% Uncomment:
FoI= pij_1*inf;%; %homogenous mixing: ones(18)*diag(pt_est) gives pij = prob sus indiv i gets inf from infect indiv j x number inf indivs in j (inf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Differential Equations:
dS=  [births;S(1:params.nA-1)]./params.Austep - FoI.*S - (Nu + Mu).*S - S.*params.AgeDur + M_S; % Susceptible nho sua giau * thanh giau /
dV=  [0;V(1:params.nA-1)]./params.Austep + Nu.*S - Mu.*V - V.*params.AgeDur ; % Vaccinated/Immune
dA=  [0;A(1:params.nA-1)]./params.Austep + FoI.*S - (params.gamma(:,1) + params.lambda(:,1) + Mu + params.alpha(:,1)).*A - A.*params.AgeDur; % Acute Infection
dI=  [0;I(1:params.nA-1)]./params.Austep  + params.alpha(:,1).*A +  params.psi(:,25).*H - ( params.gamma(:,2) + Mu + params.alpha(:,2) + params.kappa(:,1)).*I - I.*params.AgeDur + M_I; % Immune Tolerant 
dI1= [0;I1(1:params.nA-1)]./params.Austep + params.alpha(:,2).*I + params.psi(:,16).*H + (params.psi(:,1)).*C1 + (params.Eta1).*I1T - (params.gamma(:,3) + Mu + params.alpha(:,3) + params.kappa(:,2) + params.tau1 + params.omega(:,1) ).*I1  - I1.*params.AgeDur + M_I1; % Immune Clearance
dI2= [0;I2(1:params.nA-1)]./params.Austep + params.alpha(:,3).*I1 + params.psi(:,17).*H + (params.psi(:,3)).*C2 + (params.Eta2).*I2T - (params.gamma(:,4) + Mu + params.alpha(:,4) + params.kappa(:,3) + params.tau2 + params.omega(:,2)).*I2  - I2.*params.AgeDur + M_I2; % Immune Control
dI3= [0;I3(1:params.nA-1)]./params.Austep + params.alpha(:,4).*I2 + params.psi(:,18).*H + params.psi(:,5).*C3 + (params.Eta3).*I3T - (params.gamma(:,5) + Mu + params.kappa(:,4) + params.tau3 + params.omega(:,3) ).*I3  - I3.*params.AgeDur + M_I3 ; % Immune Escape
dC1= [0;C1(1:params.nA-1)]./params.Austep + params.omega(:,1).*I1 + params.psi(:,13).*H + (params.Eta4).*C1T + params.psi(:,7).*D  - (params.gamma(:,6) + Mu + params.alpha(:,5) + params.kappa(:,5) + (params.psi(:,1)) + params.tau4 + params.theta(:,1)).*C1  - C1.*params.AgeDur + M_C1; % Cirrhotic Immune Clearance
dC2= [0;C2(1:params.nA-1)]./params.Austep + params.alpha(:,5).*C1 + params.omega(:,2).*I2 + params.psi(:,14).*H + (params.Eta5).*C2T + params.psi(:,8).*D  - (params.gamma(:,7) + Mu + params.alpha(:,6) + params.kappa(:,6) + params.psi(:,3)  + params.tau5 + params.theta(:,2)).*C2  - C2.*params.AgeDur + M_C2; % Cirrhotic Immune Control
dC3= [0;C3(1:params.nA-1)]./params.Austep + params.alpha(:,6).*C2 + params.omega(:,3).*I3 + params.psi(:,15).*H + (params.Eta6).*C3T + params.psi(:,9).*D  - (params.gamma(:,8) + Mu + params.kappa(:,7) + params.psi(:,5)  + params.tau6 + params.theta(:,3)).*C3  - C3.*params.AgeDur + M_C3; % Cirrhotic Immune Escape
dI1T= [0;I1T(1:params.nA-1)]./params.Austep + params.tau1.*I1 + params.psi(:,22).*HT + (params.psi(:,2)).*C1T - (params.gamma(:,9) + Mu + params.alpha(:,7) + params.kappa(:,8) + params.omega(:,4) + (params.Eta1)).*I1T - I1T.*params.AgeDur ; % Immune Clearance on Treatment
dI2T= [0;I2T(1:params.nA-1)]./params.Austep + params.alpha(:,7).*I1T + params.tau2.*I2 + params.psi(:,23).*HT + params.psi(:,4).*C2T - (params.gamma(:,10) + Mu + params.alpha(:,8) + params.kappa(:,9) + params.omega(:,5) + (params.Eta2)).*I2T - I2T.*params.AgeDur ; % Immune Control on Treatment
dI3T= [0;I3T(1:params.nA-1)]./params.Austep + params.alpha(:,8).*I2T + params.tau3.*I3 + params.psi(:,24).*HT + params.psi(:,6).*C3T -(params.gamma(:,11) + Mu + params.kappa(:,10) + params.omega(:,6) + (params.Eta3)).*I3T - I3T.*params.AgeDur ; % Immune Escape on Treatment
dC1T= [0;C1T(1:params.nA-1)]./params.Austep + params.omega(:,4).*I1T +  params.psi(:,19).*HT + params.tau4.*C1 + params.psi(:,10).*DT -(params.gamma(:,12) + Mu + params.alpha(:,9) + params.kappa(:,11) + (params.psi(:,2)) + (params.Eta4) + params.theta(:,4)).*C1T - C1T.*params.AgeDur ; % Cirrhotic Immune Clearance on Treatment
dC2T= [0;C2T(1:params.nA-1)]./params.Austep + params.alpha(:,9).*C1T+ params.omega(:,5).*I2T + params.psi(:,20).*HT + params.tau5.*C2 + params.psi(:,11).*DT -(params.gamma(:,13) + Mu + params.alpha(:,10) + params.kappa(:,12) + params.psi(:,4) + (params.Eta5) + params.theta(:,5)).*C2T - C2T.*params.AgeDur ; % Cirrhotic Immune Control on Treatment
dC3T= [0;C3T(1:params.nA-1)]./params.Austep + params.alpha(:,10).*C2T + params.omega(:,6).*I3T + params.psi(:,21).*HT + params.tau6.*C3 + params.psi(:,12).*DT -(params.gamma(:,14) + Mu + params.kappa(:,13) + params.psi(:,6) + (params.Eta6) + params.theta(:,6)).*C3T - C3T.*params.AgeDur ; % Cirrhotic Immune Escape on Treatment

dH=   [0;H(1:params.nA-1)]./params.Austep   + params.kappa(:,1).*I + params.kappa(:,2).*I1 + params.kappa(:,3).*I2 + params.kappa(:,4).*I3 + params.kappa(:,5).*C1 + params.kappa(:,6).*C2 + params.kappa(:,7).*C3  + (params.Eta7).*HT + params.kappa(:,14).*D - ( params.gamma(:,15) + Mu + params.theta(:,7) + params.lambda(:,2) +  params.psi(:,13) +  params.psi(:,14) +  params.psi(:,15) + params.psi(:,16) +  params.psi(:,17) + params.psi(:,18) + params.tau7 + params.psi(:,25) ).*H  - H.*params.AgeDur ;  % Hepatocellular Carcinoma %+ (1./params.gamma(:,15)).*R
dHT=  [0;HT(1:params.nA-1)]./params.Austep  + params.kappa(:,8).*I1T + params.kappa(:,9).*I2T + params.kappa(:,10).*I3T + params.kappa(:,11).*C1T + params.kappa(:,12).*C2T + params.kappa(:,13).*C3T  + params.tau7.*H + params.kappa(:,15).*DT - (params.gamma(:,16) + Mu + (params.Eta7) + params.theta(:,8) + params.lambda(:,3) +  params.psi(:,19) +  params.psi(:,20) +  params.psi(:,21) +  params.psi(:,22) +  params.psi(:,23) +  params.psi(:,24) ).*HT  - HT.*params.AgeDur ; % Hepatocellular Carcinoma on Treatment  %+ (1./params.gamma(:,16)).*R %+ params.tau7.*H
dD=   [0;D(1:params.nA-1)]./params.Austep  + params.theta(:,1).*C1 + params.theta(:,2).*C2 + params.theta(:,3).*C3 + params.theta(:,7).*H + (params.Eta8).*DT  - ( params.tau8 + Mu + params.lambda(:,4) + params.psi(:,7) + params.psi(:,8) + params.psi(:,9) + params.kappa(:,14)).*D  - D.*params.AgeDur ;  % Decompensated Cirrhosis
dDT=  [0;DT(1:params.nA-1)]./params.Austep  + params.theta(:,4).*C1T + params.theta(:,5).*C2T + params.theta(:,6).*C3T + params.theta(:,8).*HT + params.tau8.*D - (params.kappa(:,15) + (params.Eta8) + Mu + params.lambda(:,5) + params.psi(:,10) + params.psi(:,11) + params.psi(:,12) + params.gamma(:,17)).*DT  - DT.*params.AgeDur ; % Decompensated Cirrhosis on Treatment %+ params.tau8.*D
dR=   [0;R(1:params.nA-1)]./params.Austep   + params.gamma(:,1).*A + params.gamma(:,2).*I + params.gamma(:,3).*I1 + params.gamma(:,4).*I2 + params.gamma(:,5).*I3 + params.gamma(:,6).*C1 + params.gamma(:,7).*C2 + params.gamma(:,8).*C3 + params.gamma(:,9).*I1T + params.gamma(:,10).*I2T + params.gamma(:,11).*I3T + params.gamma(:,12).*C1T + params.gamma(:,13).*C2T + params.gamma(:,14).*C3T + params.gamma(:,15).*H + params.gamma(:,16).*HT + params.gamma(:,17).*DT - Mu.*R   - R.*params.AgeDur + M_R;  % Recovered/Resolved HBV infection


% Cumulative infections and mortality:  
cum_inc_ac = FoI.*S; %Cumulative Incidence of Acute Infections
cum_inc_ImTol = params.alpha(:,2).*I; %Cumulative Incidence of Chronic Infections %Escape!
% cum_inc_HCC = 0 ; %Cumulative Incidence of HCC cases
% cum_inc_DC = 0 ; %Cumulative Incidence of DC cases
HCC_deaths = params.lambda(:,2).*H; %cumulative HCC deaths
HCC_deaths_T = params.lambda(:,3).*HT; %cumulative HCC deaths on T
DC_deaths = params.lambda(:,4).*D; %cumulative DC deaths
DC_deaths_T =  params.lambda(:,5).*DT; %cumulative DC deaths on T
bg_deaths = Mu.*(I+I1+I2+I3+C1+C2+C3+I1T+I2T+I3T+C1T+C2T+C3T+D+DT+H+HT); %bg moret in chronic to other causes
acute_deaths = params.lambda(:,1).*A;
HCC_incidence = params.kappa(:,1).*I + params.kappa(:,2).*I1 + params.kappa(:,3).*I2 + params.kappa(:,4).*I3 + params.kappa(:,5).*C1 + params.kappa(:,6).*C2 + params.kappa(:,7).*C3 + params.kappa(:,14).*D; %+ (params.Eta7).*HT
HCC_T_incidence =  params.kappa(:,8).*I1T + params.kappa(:,9).*I2T + params.kappa(:,10).*I3T + params.kappa(:,11).*C1T + params.kappa(:,12).*C2T + params.kappa(:,13).*C3T + params.kappa(:,15).*DT; % + params.tau7.*H 
DC_incidence = params.theta(:,1).*C1 + params.theta(:,2).*C2 + params.theta(:,3).*C3 + params.theta(:,7).*H ; %+ (params.Eta8).*DT
DC_T_incidence = params.theta(:,4).*C1T + params.theta(:,5).*C2T + params.theta(:,6).*C3T + params.theta(:,8).*HT ; %+ params.tau8.*D
Treat_uptake = params.tau4.*C1; %params.tau1.*I1 + params.tau2.*I2 + params.tau3.*I3 + params.tau4.*C1 + params.tau5.*C2 + params.tau6.*C3 + params.tau7.*H + params.tau8.*D;

rates = [dS;dV;dA;dI;dI1;dI2;dI3;dC1;dC2;dC3;dI1T;dI2T;dI3T;dC1T;dC2T;dC3T;dH;dHT;dD;dDT;dR ;cum_inc_ac; cum_inc_ImTol; HCC_deaths; HCC_deaths_T; DC_deaths; DC_deaths_T; bg_deaths; acute_deaths; HCC_incidence; HCC_T_incidence; DC_incidence; DC_T_incidence; Treat_uptake];%HCC_I; HCC_I1; HCC_I2; HCC_I3; HCC_C1; HCC_C2; HCC_C3; HCC_DC]; %cum_inc_HCC; cum_inc_DC];

%Output HCC Incidence by Phase
%rates = [dS;dV;dA;dI;dI1;dI2;dI3;dC1;dC2;dC3;dI1T;dI2T;dI3T;dC1T;dC2T;dC3T;dH;dHT;dD;dDT;dR ;cum_inc_ac; cum_inc_ImTol; HCC_deaths; HCC_deaths_T; DC_deaths; DC_deaths_T; bg_deaths; acute_deaths; HCC_incidence; HCC_T_incidence; DC_incidence; DC_T_incidence; Treat_uptake;HCC_I; HCC_I1; HCC_I2; HCC_I3; HCC_C1; HCC_C2; HCC_C3; HCC_DC];

    function [tau_treat] = treat(x)
      x=x-50;  
      a=90; %20
      b=2; %high b results in slower increase %2.1 increase too slow %orig = 2 %1.8 on lhs...good for pbs public
      c=0.1; %a =20, b=2, c=0.16/0.17 plus surv_params.pt = [0.09, 0, 0.1, 0; 0.09, 0.06, 0.36, 0.3]; slighty under for Nat Est but good shape!
      tau_treat = (a*c*exp(-c*x))./(exp(-c*x + b).^2);
    end

 % S curve like growth for treatment uptake in future:
function [tau_treat] = treat2021(x)  %current best x-45, a=18, b=1.6,c=0.17 (max 20%)/ a=24, b=1.8 good (max 23%)/ a=30, b=1.8 good (max (25%)/a=50,b=2.5 reaches 28%
        x = x-45;      
        a= 12;%12 %14; %17; %23; %27; %30; % a shifts plateau up or down
        b= 0.75; %0.7;%0.9 %1.2; %1.7; %1.9; %2.2;  %high b results in slower increase. %1.6 perfect %was 1.55;
        c= 0.18; %0.17;%a =20, b=2, c=0.16/0.17 plus surv_params.pt = [0.09, 0, 0.1, 0; 0.09, 0.06, 0.36, 0.3]; slighty under for Nat Est but good shape!
        tau_treat = (a*c)./(1 + exp(-c.*x + 2*b));
end


function [tau_treat] = treat_proj(x)  % future tx inc similar to end 2015 to 2020
        x = x-45;      
        a= 180; %39; %33.5; %37; %36; % a shifts plateau up or down
        b= 2.3; %2.95; %2.85; %3.1; %3.1; %2.9; %high b results in slower increase. %1.6 perfect %was 1.55;
        c= 0.1; %0.19; %0.18; 
        tau_treat = (a*c)./(1 + exp(-c.*x + 2*b));
end



end

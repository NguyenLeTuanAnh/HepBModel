function [ params ] = Nat_Params()

model_labels = {'ACT', 'NSW','NT','QLD','SA','TAS','VIC','WA','National'}; 
idx = 9; % 1 - ACT; 2- NSW; 3-NT; 4-QLD; 5-SA; 6-TAS; 7-VIC; 8-WA; 9-AUS
params.idx = 9;

%% Treatment included or excluded: t=0 excluded; t=1 included.
tt = 1; %1 %0.48; %0.48; %0.42; %0.48; 
vt = 1; %vacc on/off

% Births: per year
% Historical for 1951-2016 inclusive & ABS series B projections from 2017 onwards
pop_params.births = xlsread('Inputs/population.xlsx',model_labels{idx},'A2:A101');


% Total Population (1951 - 2050)
pop_params.population = xlsread('Inputs/population.xlsx',model_labels{idx},'B2:B101'); %dlmread('Inputs/population2018_National.csv',',',[1 1 100 1]);

% Historical age distribution of Australia 1951-2017 
% ABS Population Age Groups:
% 0-4;5-9;10-14;15-19;20-24;25-29;30-34;35-39;40-44;45-49;50-54;55-59;60-64
% 65-69;70-74;75-79;80-84;85+
pop_params.AgeDist = xlsread('Inputs/PopAgeDist.xlsx',model_labels{idx},'B3:S69'); % dlmread('Inputs/PopAgeDist_National.csv',',',[2 1 68 18]); % row2 = 1951, row3=1952 etc.

% Background Mortality Rates: mu(a,t) = rate for ages group a in year t.
pop_params.mu = xlsread('Inputs/MortalityRates.xlsx',model_labels{idx},'B3:S102'); %dlmread('Inputs/MortalityRates_Nat_adj.csv',',',[2 1 101 18]);

% % Assume forward projections of probability of dying in each age group
% % are the same as in 2016. 
pop_params.mu(1:68,:) = 1*pop_params.mu(1:68,:);
pop_params.mu(69:end,:) = ones(32,1)*pop_params.mu(68,:); %[2:(20-2)/32:20]'.
pop_params.mu(50:end,:) = 0.85*pop_params.mu(50:end,:);

pop_params.mu=pop_params.mu'./12;

% Run model with FoI = foi_pt_allages' (static foi) and extract number of
% infectious units in the year 2000 (peak of notifs).
% inf_ages_2000= round([A(50,:)' I(50,:)' I1(50,:)' I2(50,:)' I3(50,:)' C1(50,:)' C2(50,:)' C3(50,:)' H(50,:)' D(50,:)']*eps)
params.inf_units_agegroups_2000 = xlsread('Inputs/inf_units_age_2000.xlsx',model_labels{idx},'A1:A18'); %dlmread('Inputs/inf_units_agegroups_2000.csv',',',[0 0 17 0]);


% Migration: by age and infection state per year
params.prev_COB_1951to1974 = dlmread('Inputs/Prev_ArrCOB_1951to1974_BC.csv',',',[1 1 35 1]);
params.prev_COB_1975to1990 = dlmread('Inputs/Prev_MigCOB_1975to1990_BC.csv',',',[1 1 65 1]);

params.low_prev = 0.005;
params.int_prev = 0.05;
params.high_prev = 0.1;

params.DSS_MigPrev = 6; %mp(params.x); %(1 - Jenn's, 4-Lower CI or
% % 5-Upper CI)
params.DSS_MigPrev2019 = 7; %2019 syst review update

params.dis_phase_cirr1 = 0; %age group prop cirrhotic  <20yrs
params.dis_phase_cirr2 = 0.03; %age group prop cirrhotic 20-40yo
params.dis_phase_cirr3 = 0.05; %age group prop cirrhotic >40yo %prev 0.06


% Notifications: per year 1971 to 2020 incl.
pop_params.notifs_time = zeros(100,1);
pop_params.notifs1971to2021 = dlmread('Inputs/notifications.csv',',',[0 1 50 1]); % change var to 1971to2021
pop_params.cum_notifs1971to2021 = dlmread('Inputs/notifications-cum.csv',',',[0 1 50 1]); % 49-> 50 to include cum notif of 2021
pop_params.notifs_time(21:71) = pop_params.notifs1971to2021;


% Treatment: per year
% create one csv file for Nat + STT
pop_params.treated = dlmread('Inputs/PBS_Nat_Treat2000to2021.csv',',',[0 1 99 1]);
pop_params.treated2000to2021 = pop_params.treated(50:71);


%% Simulation Parameters (i.e. time related) %%
%(yearly_rate/365)*(days_month) = new param rate per month
sim_params.days_month = [31 28 31 30 31 30 31 31 30 31 30 31];
sim_params.time_start=0*12; %1951 %timestep in months
sim_params.time_end=100*12; %2050
params.time_entries= (sim_params.time_end - sim_params.time_start)/12; %reported timesteps in years


%% Demographic Parameters  %%
pop_params.ages = 0:5:86; % by 2months up to 6months; 6 months to 1 year; every 5 yrs.
params.nA = length(pop_params.ages); 
pop_params.nA = length(pop_params.ages); %Number of age groups
pop_params.births_months = ((pop_params.births./365)*sim_params.days_month)'; 
pop_params.births_m = pop_params.births_months(:);
pop_params.pop_months = ((pop_params.population./365)*sim_params.days_month)'; %12 by 100 matrix. 
pop_params.AgeDur = [(1./([diff(0:5:81),5]*12))';0]; %duration in months of age groups
pop_params.Austep = 1./[1;pop_params.AgeDur(1:pop_params.nA-1)];
pop_params.birth_rate = pop_params.births_months./pop_params.pop_months;


%% Initialize proportion of population in each age group who are in
%  Susceptible, Immune Tolerant and Resolved states.
%  init0(i,g) = prop of population in infection state i and in age group g.
%  i = 1-Susceptible, 2-Immune Tolerant, 3-Resolved
%  See Table 5.8, page 155 in B.Cowie thesis.
pop_params.init0 = dlmread('Inputs/Init0_Nat.csv',',',[1 1 6 18]);


%% Disease Parameters %%
dis_params.foi_scale = 1;
dis_params.sig = sig_av(pop_params.ages);
%Rate of resolution of HBV infection in each age group
dis_params.ac_duration = 3;
dis_params.ac_res = 1./dis_params.ac_duration; % 1/dis_params.ac_res = 3months = av duration of acute infection. =4/12;
dis_params.gamma= zeros(pop_params.nA,16)./12; 
dis_params.gamma(:,1)= dis_params.ac_res*(1- dis_params.sig./12); %1./sig_av(pop_params.ages)./12; %A -> R % gamma(nA,i) = gamma_i etc.
dis_params.gamma(:,2)=0./12; %I -> R *

dis_params.gamma_a = [0.0016*ones(6,1);   0.0022*ones(2,1); 0.0034*ones(2,1); 0.0038*ones(8,1)]./12;
dis_params.gamma_b = [0.0077*ones(6,1);  0.0107*ones(2,1); 0.0165*ones(2,1); 0.0183*ones(8,1)]./12;

dis_params.gamma(:,3) = dis_params.gamma_a; %I1 -> R
dis_params.gamma(:,4)= dis_params.gamma_b; %I2 -> R
dis_params.gamma(:,5)= dis_params.gamma_b; %I3 -> R
dis_params.gamma(:,6)= dis_params.gamma_a; %C1 -> R
dis_params.gamma(:,7)= dis_params.gamma_b;%C2 -> R
dis_params.gamma(:,8)= dis_params.gamma_b; %C3 -> R
dis_params.gamma(:,9)= dis_params.gamma_a;%I1T -> R
dis_params.gamma(:,10)= dis_params.gamma_b; %I2T -> R
dis_params.gamma(:,11)= dis_params.gamma_b; %I3T -> R
dis_params.gamma(:,12)= dis_params.gamma_a; %C1T -> R
dis_params.gamma(:,13)= dis_params.gamma_b; %C2T -> R
dis_params.gamma(:,14)= dis_params.gamma_b; %C3T -> R
dis_params.gamma(:,15)= 0; %H -> R
dis_params.gamma(:,16)= 0.035./12; %HT -> R
dis_params.gamma(:,17)= 0.035./12; %DT -> R 


%proportion of HCC cases that go back to chronic state of origin.
dis_params.hr=0.1; 
%proportion of cirrhotic cases that revert back to non-cirrhotic
dis_params.cr=0.1; 
dis_params.psi_fact1 = 0.15;
dis_params.psi_fact2 = 0.003;
dis_params.psi(:,1) = zeros(pop_params.nA,1); %C1 -> I1 
dis_params.psi(:,2) = (dis_params.psi_fact1*ones(pop_params.nA,1))./12; %C1T -> I1T 
dis_params.psi(:,3) = (dis_params.psi_fact2*ones(pop_params.nA,1))./12; %C2 -> I2 
dis_params.psi(:,4) = (dis_params.psi_fact2*ones(pop_params.nA,1))./12; %C2T -> I2T **new was _fact1
dis_params.psi(:,5) = zeros(pop_params.nA,1); %C3 ->I3 
dis_params.psi(:,6) = (dis_params.psi_fact1*ones(pop_params.nA,1))./12; %C3T ->I3T 
dis_params.psi(:,7) = zeros(pop_params.nA,1); %D -> C1 
dis_params.psi(:,8) = zeros(pop_params.nA,1); %D -> C2 
dis_params.psi(:,9) = zeros(pop_params.nA,1); %D -> C3 
dis_params.psi(:,10) = zeros(pop_params.nA,1); %DT -> C1T 
dis_params.psi(:,11) = (dis_params.psi_fact1*ones(pop_params.nA,1))./12; %DT -> C2T 
dis_params.psi(:,12) = zeros(pop_params.nA,1); %DT -> C3T

% No longer incorporating reversion of HCC
dis_params.psi(:,13) = zeros(pop_params.nA,1); %H -> C1 
dis_params.psi(:,14) = zeros(pop_params.nA,1); %H -> C2 
dis_params.psi(:,15) = zeros(pop_params.nA,1); %H -> C3 
dis_params.psi(:,16) = zeros(pop_params.nA,1); %H -> I1 
dis_params.psi(:,17) = zeros(pop_params.nA,1); %H -> I2 
dis_params.psi(:,18) = zeros(pop_params.nA,1); %H -> I3 
dis_params.psi(:,19) = zeros(pop_params.nA,1); %HT -> C1T
dis_params.psi(:,20) = zeros(pop_params.nA,1); %HT -> C2T 
dis_params.psi(:,21) = zeros(pop_params.nA,1); %HT -> C3T 
dis_params.psi(:,22) = zeros(pop_params.nA,1); %HT -> I1T 
dis_params.psi(:,23) = zeros(pop_params.nA,1); %HT -> I2T 
dis_params.psi(:,24) = zeros(pop_params.nA,1); %HT -> I3T 
dis_params.psi(:,25) = zeros(pop_params.nA,1); %H -> I 

%proportion of DC cases that revert back to chronic state of origin.
dis_params.dc=0.1; 
%dis_params.p_alpha = 1; 

%Progression rates through HBV phases
dis_params.alpha= zeros(pop_params.nA,10)./12; 
dis_params.alpha(:,1)= dis_params.ac_res*(dis_params.sig./12); % A -> I
%dis_params.alpha(:,2)= [ones(2,1)*0.03; ones(4,1)*0.02; ones(2,1)*0.03; ones(10,1)*0.06]./12;  % I -> I1 **
dis_params.alpha2 = ( [ones(2,1)*0.03; ones(4,1)*0.02; ones(2,1)*0.03; ones(2,1)*0.06; 0.07*ones(2,1); ones(6,1)*0.08]./12);
dis_params.alpha(:,2)= dis_params.alpha2;  % I -> I1 **
dis_params.alpha3 = (  [0.05; 0.12*ones(3,1); 0.184*ones(3,1); 0.136*ones(3,1); 0.15*ones(8,1) ]./12); %1.1***
dis_params.alpha(:,3)= dis_params.alpha3; % I1 -> I2  
dis_params.alpha4 = (  [0.0087*ones(6,1); 0.0149*ones(2,1); 0.0278*ones(2,1); 0.0202*ones(8,1)]./12); % I2 -> I3  0.8***
dis_params.alpha(:,4) = dis_params.alpha4;
dis_params.alpha5 =  ( [0.05; 0.12*ones(3,1); 0.184*ones(3,1); 0.136*ones(3,1); 0.15*ones(8,1) ]./12); % C1 -> C2 1.1*
dis_params.alpha(:,5) = dis_params.alpha5;
dis_params.alpha6 = ( [0.0087*ones(6,1); 0.0149*ones(2,1); 0.0278*ones(2,1); 0.0202*ones(8,1)]./12);  % C2 -> C3 
dis_params.alpha(:,6)= dis_params.alpha6;
dis_params.alpha7 = (  [0.11; 0.2664*ones(3,1); 0.4085*ones(3,1); 0.3019*ones(3,1); ones(8,1)*0.33]./12); % I1T -> I2T 1.1***
dis_params.alpha(:,7) = dis_params.alpha7;
dis_params.alpha(:,8)= zeros(18,1)./12; % I2T -> I3T
dis_params.alpha9 = (  [0.11; 0.2664*ones(3,1); 0.4085*ones(3,1); 0.3019*ones(3,1); ones(8,1)*0.33]./12 );% C1T -> C2T 1.1*
dis_params.alpha(:,9)= dis_params.alpha9;
dis_params.alpha(:,10)=(zeros(18,1)./12); % C2T -> C3T


% %% June 14th HCC updated rates 1:
% Progression rates to HCC
z1 = 1; 
z2 = 1; 
dis_params.kappa= zeros(pop_params.nA,13)./12;
dis_params.kappa1 = [0.0001*ones(6,1); 0.0005*ones(2,1); 0.0008*ones(2,1); 0.001*ones(8,1)]./12;
dis_params.kappa(:,1) = dis_params.kappa1; %I -> H  %[0.0005*ones(4,1); 0.001*ones(4,1); 0.0031*ones(10,1)]./12; %I -> H
                       %[0.0005*ones(4,1); 0.001*ones(4,1); 0.0031*ones(10,1)]./12; %I -> H
dis_params.kappa2 = [z1*0.0005*ones(6,1); z1*0.001*ones(2,1); z2*0.004*ones(2,1); 0.005*ones(2,1); 0.007*ones(6,1)]./12;
dis_params.kappa(:,2) = dis_params.kappa2 ; %I1 -> H
dis_params.kappa3 = [0.0001*ones(6,1); 0.0005*ones(2,1); 0.0008*ones(2,1); 0.001*ones(2,1); 0.002*ones(6,1)]./12;
dis_params.kappa(:,3)= dis_params.kappa3; %I2 -> H %0.002
dis_params.kappa4 = [z1*0.0005*ones(6,1); z1*0.001*ones(2,1); z2*0.003*ones(2,1); 0.005*ones(2,1); 0.007*ones(6,1)]./12; %0.9*
dis_params.kappa(:,4) = dis_params.kappa4;  %0.01*ones(18,1)./12; %I3 -> H

dis_params.kappa5 = ([z1*0.0035*ones(6,1); z1*0.007*ones(2,1);  z2*0.028*ones(2,1); 0.035*ones(2,1); 0.05*ones(6,1) ]./12);
dis_params.kappa(:,5)= dis_params.kappa5; %0.037*ones(18,1)./12; %C1 -> H
dis_params.kappa6 = ([z1*0.0035*ones(6,1); z1*0.0056*ones(2,1); z2*0.007*ones(2,1); 0.014*ones(2,1); 0.018*ones(6,1)]./12);
dis_params.kappa(:,6) = dis_params.kappa6; %0.02*ones(18,1)./12; %C2 -> H

%dis_params.kappa7 = ([z1*0.0035*ones(6,1); z1*0.007*ones(2,1);  z2*0.023*ones(2,1); 0.03*ones(2,1); 0.045*ones(6,1)]./12);
dis_params.kappa7 = ([z1*0.003*ones(6,1); z1*0.006*ones(2,1);  z2*0.018*ones(2,1); 0.025*ones(2,1); 0.04*ones(6,1)]./12);
dis_params.kappa(:,7) = dis_params.kappa7; %0.03*ones(18,1)./12; %C3 -> H

tx_adj=0.75; %was 0.5 then 0.75 after BC comments
%dis_params.kappa8 = (0.5*[z1*0.0005*ones(6,1); z1*0.001*ones(2,1); z2*0.004*ones(2,1); 0.005*ones(2,1); 0.007*ones(6,1)])./12;
dis_params.kappa8 = (tx_adj*[z1*0.0005*ones(6,1); z1*0.001*ones(2,1); z2*0.004*ones(2,1); 0.005*ones(2,1); 0.007*ones(6,1)])./12;
dis_params.kappa(:,8) = dis_params.kappa8 ; %I1T -> HT
%dis_params.kappa9 = (0.5*[0.0001*ones(6,1); 0.0005*ones(2,1); 0.0008*ones(2,1); 0.001*ones(2,1); 0.002*ones(6,1)])./12;
dis_params.kappa9 = (tx_adj*[0.0001*ones(6,1); 0.0005*ones(2,1); 0.0008*ones(2,1); 0.001*ones(2,1); 0.002*ones(6,1)])./12;
dis_params.kappa(:,9) = dis_params.kappa9; %0.001*ones(18,1)./12; %I2T -> HT

%dis_params.kappa10 = (0.5*[z1*0.0005*ones(6,1); z1*0.001*ones(2,1); z2*0.003*ones(2,1); 0.005*ones(2,1); 0.007*ones(6,1)])./12;
dis_params.kappa10 = (tx_adj*[z1*0.0005*ones(6,1); z1*0.001*ones(2,1); z2*0.003*ones(2,1); 0.005*ones(2,1); 0.007*ones(6,1)])./12;
dis_params.kappa(:,10)= dis_params.kappa10; %0.005*ones(18,1)./12; %I3T -> HT

%dis_params.kappa11 = (0.5*([z1*0.0035*ones(6,1); z1*0.007*ones(2,1);  z2*0.028*ones(2,1);  0.035*ones(2,1); 0.05*ones(6,1) ])./12);
dis_params.kappa11 = (tx_adj*([z1*0.0035*ones(6,1); z1*0.007*ones(2,1);  z2*0.028*ones(2,1);  0.035*ones(2,1); 0.05*ones(6,1) ])./12);
dis_params.kappa(:,11)= dis_params.kappa11; %0.017*ones(18,1)./12; %C1T -> HT

%dis_params.kappa12 = (0.5*[z1*0.0035*ones(6,1); z1*0.0056*ones(2,1); z2*0.007*ones(2,1); 0.014*ones(2,1); 0.018*ones(6,1)]./12); 
dis_params.kappa12 = (tx_adj*[z1*0.0035*ones(6,1); z1*0.0056*ones(2,1); z2*0.007*ones(2,1); 0.014*ones(2,1); 0.018*ones(6,1)]./12); 
dis_params.kappa(:,12)= dis_params.kappa12; %0.01*ones(18,1)./12; %C2T -> HT

%dis_params.kappa13 = (0.5*[z1*0.0035*ones(6,1); z1*0.007*ones(2,1);  z2*0.023*ones(2,1); 0.03*ones(2,1); 0.045*ones(6,1)])./12;
%dis_params.kappa13 = (0.5*[z1*0.003*ones(6,1); z1*0.006*ones(2,1);  z2*0.018*ones(2,1); 0.025*ones(2,1); 0.04*ones(6,1)])./12;
% Adj after discussion with BC:
dis_params.kappa13 = (tx_adj*[z1*0.003*ones(6,1); z1*0.006*ones(2,1);  z2*0.018*ones(2,1); 0.025*ones(2,1); 0.04*ones(6,1)])./12;
dis_params.kappa(:,13) = dis_params.kappa13;  % 0.015*ones(18,1)./12; %C3T -> HT

dis_params.kappa14 = 0.048*ones(18,1)./12;
dis_params.kappa(:,14) = dis_params.kappa14;  %D -> H
%dis_params.kappa15 = 0.024*ones(18,1)./12;
dis_params.kappa15 = tx_adj*0.048*ones(18,1)./12;
dis_params.kappa(:,15) = dis_params.kappa15; %DT -> HT
%dis_params.kappa=0.25*dis_params.kappa;


% Progression rates to Decompensated Cirrhosis
dis_params.theta1 = 0.033;
dis_params.theta2 = 0.01; 
dis_params.theta3 = 0.0149;
dis_params.theta4 = 0.0045;
dis_params.theta = zeros(pop_params.nA,8)./12;
dis_params.theta(:,1) = dis_params.theta1*ones(18,1)./12; %C1 -> D 
dis_params.theta(:,2) = dis_params.theta2*ones(18,1)./12; %C2 -> D
dis_params.theta(:,3) = dis_params.theta1*ones(18,1)./12; %C3 -> D
dis_params.theta(:,4) = dis_params.theta3*ones(18,1)./12; %C1T -> DT
dis_params.theta(:,5) = dis_params.theta4*ones(18,1)./12; %C2T -> DT
dis_params.theta(:,6) = dis_params.theta3*ones(18,1)./12; %C3T -> DT
dis_params.theta(:,7) = zeros(18,1)./12; %H -> D
dis_params.theta(:,8) = zeros(18,1)./12; %HT -> DT


%Rates of developing cirrhosis
dis_params.omega = zeros(pop_params.nA,6)./12; 
dis_params.omega1 = 0.016;
dis_params.omega(:,1) = (ones(18,1)*dis_params.omega1)./12; %I1 -> C1
%dis_params.omega(:,2) = 0.9*[0.00038*ones(6,1); 0.00049*ones(2,1); 0.00068*ones(2,1); 0.00150*ones(8,1)]./12; %I2 -> C2
dis_params.omega2 = [0.00038*ones(6,1); 0.00049*ones(2,1); 0.00068*ones(2,1); 0.00150*ones(8,1)]./12;
dis_params.omega(:,2) = dis_params.omega2; %I2 -> C2
dis_params.omega3 = [0.0026*ones(4,1); 0.0052*ones(2,1); 0.0078*ones(2,1); 0.0261*ones(2,1); 0.0414*ones(8,1)]./12;
dis_params.omega(:,3) = dis_params.omega3; %I3 -> C3
dis_params.omega4 = 0.0088;
dis_params.omega(:,4) = (ones(18,1)*dis_params.omega4)./12; %I1T -> C1T
dis_params.omega(:,5) = zeros(18,1)./12; %I2T -> C2T
dis_params.omega6 = [0.00143*ones(4,1); 0.00286*ones(2,1); 0.00429*ones(2,1); 0.0144*ones(2,1); 0.0228*ones(8,1)]./12;
dis_params.omega(:,6) = dis_params.omega6; %I3T -> C3T


%Mortality due to HBV infection
dis_params.lambda = zeros(pop_params.nA,5);
%dis_params.l1 = ([0.14; 0.14; 0.14; 0.350*ones(15,1)])./12; 
dis_params.l1 = ([0.001; 0.0014; 0.0014; 0.00350*ones(15,1)])./12; % BC model: [0.001, 0.0014, 0.0035, 0.0035] per year
dis_params.lambda(:,1) = dis_params.l1; %Ac HBV mortality
dis_params.l2 = 0.33;
dis_params.lambda(:,2) = dis_params.l2*ones(18,1)./12; %H mortality put back to 0.33 from 0.3 Aug 13
dis_params.l3 = 0.2475; %0.33*0.75 %or half non-treat estimate = 0.1650? % orig = 0.1213 
dis_params.lambda(:,3) = dis_params.l3*ones(18,1)./12; %HT mortality put back to 0.1213 from 0.1 Aug 13
dis_params.l4 = 0.2;
dis_params.lambda(:,4) = dis_params.l4*ones(18,1)./12; %D mortality
dis_params.lambda(:,5) = (0.5*dis_params.l4)*ones(18,1)./12; %DT mortality 


% Vaccination: different proportion of coverage from 1999 to 2018
%Proportion of vaccinated individuals that attain immunity 
surv_params.vacc_prop= [0.9, 0.6, 0.05, 0.03, 0];
surv_params.coverage = dlmread('Inputs/AustraliaSTT1999to2021InfantVacc.csv',',',[1 idx 23 idx])'./100; 
surv_params.efficacy = [0.95*ones(3,1); 0.9*ones(6,1); 0.75*ones(5,1); 0.7*ones(4,1)];

surv_params.nu = zeros(pop_params.nA,params.time_entries); %vacc coverage by age group

% Historical Infant Vacc Coverage & Assuming it will stay around 95% moving
% forward from 2017 tp 2050:
surv_params.nu(1,49:71) = (surv_params.coverage*surv_params.efficacy(1)); % 49:71 = 2000 - 2021 (vaccine efficacy: 0.95 for 0-14 from 2000-2021)
surv_params.nu(1,72:end) = (surv_params.coverage(end)*ones(1,29)*surv_params.efficacy(1)); 
% Adolescant Catch Up Vaccination:
surv_params.nu(3,49:end) = surv_params.efficacy(3)'*ones(1,52)*surv_params.vacc_prop(2);
surv_params.nu(4,49:end) = surv_params.efficacy(4)'*ones(1,52)*surv_params.vacc_prop(2);
% Adult & Maternal Vaccination:
surv_params.nu(5:7,49:end) = (surv_params.efficacy(5)*surv_params.vacc_prop(3))*ones(3,52);
surv_params.nu(8:9,49:end) = (surv_params.efficacy(5)*surv_params.vacc_prop(4))*ones(2,52);

surv_params.nu = vt*surv_params.nu./12;


%% Rates of treatment uptake - time varying. 
% Proportion of total people on treatment in phase (I1, I2, I3):
% Row1 = non-cirrhotic; %Row2 = cirrotic
%surv_params.pt = [20, 0, 25, 0; 20 , 30*0.5 , 1.5*60, 75]./100; %Aug 13 changed
surv_params.pt = [0.08, 0, 0.04, 0; 0.13, 0.15, 0.22, 0.45];  %Im Es was 0.05. just incr from 0.03 (tx to was 18.5k)


%% Fits PBS well and less bumpy than above. 
params.treat_age = ones(length(pop_params.ages),1)*1/12; %t_age(pop_params.ages)./12; 
surv_params.T1 = tt*( [zeros(18,49), (params.treat_age.*surv_params.pt(1,1))*ones(1,51) ]);   %I1 -> I1T **
surv_params.T2 = tt*( [zeros(18,49), (params.treat_age.*surv_params.pt(1,2))*ones(1,51) ]);  %I2 -> I2T (1/0.15)
surv_params.T3 = tt*( [zeros(18,49), (params.treat_age.*surv_params.pt(1,3))*ones(1,51) ]);    %I3 -> I3T
surv_params.T4 = tt*( [zeros(18,49), (params.treat_age.*surv_params.pt(2,1))*ones(1,51) ]);   %C1 -> C1T
surv_params.T5 = tt*( [zeros(18,49), (params.treat_age.*surv_params.pt(2,2))*ones(1,51) ]);     %C2 -> C2T  0.5
surv_params.T6 = tt*( [zeros(18,49), (params.treat_age.*surv_params.pt(2,3))*ones(1,51) ]);      %C3 -> C3T  1.5
%surv_params.tau7 = tt*( [zeros(18,49), (surv_params.pt(2,4)*surv_params.rho(:,50:end)) ./(12*HCC_den(:,50:end))]); %H -> HT
%surv_params.tau8 = ( [zeros(18,49), (surv_params.pt(1,4)*ones(pop_params.nA,51)) ]); %D -> DT

surv_params.T8 = tt*( [zeros(18,49), (params.treat_age.*surv_params.pt(2,4))*ones(1,51)]); %D -> DT
%surv_params.T7 = tt*[zeros(pop_params.nA,49), ones(pop_params.nA,51)*diag([0:(0.1-0)/5:0.1, 0.12:(0.38-0.12)/5:0.38, 0.4*ones(1,39)])]./12; %H -> HT
%surv_params.T7 = tt*[zeros(pop_params.nA,49), ones(pop_params.nA,51)*diag([0:(0.1-0)/5:0.1, 0.12:(0.35-0.12)/5:0.35, 0.36*ones(1,39)])]./12; %H -> HT
surv_params.T7 = tt*[zeros(pop_params.nA,49), ones(pop_params.nA,51)*diag([0:(0.1-0)/5:0.1, 0.12:(0.35-0.12)/8:0.35, 0.35*ones(1,36)])]./12; %H -> HT


zz =1;
% Rates of going off treatment:
surv_params.eta_scale = 0.15;
surv_params.eta_scale2 = 0.1;
surv_params.Eta1= zz*[zeros(pop_params.nA,49), surv_params.eta_scale*ones(pop_params.nA,51)]./12; %IT1 -> I1
surv_params.Eta2= zz*[zeros(pop_params.nA,49), surv_params.eta_scale*ones(pop_params.nA,51)]./12; %I2T -> I2
surv_params.Eta3= zz*[zeros(pop_params.nA,49), surv_params.eta_scale*ones(pop_params.nA,51)]./12; %I3T -> I3
surv_params.Eta4= zz*[zeros(pop_params.nA,49), surv_params.eta_scale*ones(pop_params.nA,51)]./12; %CT1 -> C1
surv_params.Eta5= zz*[zeros(pop_params.nA,49), surv_params.eta_scale*ones(pop_params.nA,51)]./12; %CT2 -> C2
surv_params.Eta6= zz*[zeros(pop_params.nA,49), surv_params.eta_scale*ones(pop_params.nA,51)]./12; %CT3 -> C3
surv_params.Eta7= zz*[zeros(pop_params.nA,49), zeros(pop_params.nA,51)]./12; %HT -> H
surv_params.Eta8= zz*[zeros(pop_params.nA,49), surv_params.eta_scale2*ones(pop_params.nA,51)]./12; %DT -> D


%%% Read in data for

%-----------------------------------------------------------------%
% Make ABS_MigrationEstimates_2004onwards.m more efficient:
%-----------------------------------------------------------------%
years = 2004:2050;
age_index_l = [4:5:85, 85]; %[9:10:59, 60];
age_index_u = (0:5:85); 
params.yj_a =zeros(47,18);
params.yj_b =zeros(47,18);

for y=1:47
    for jy=1:18
    params.yj_a(y,jy) = length(1919:years(y)-age_index_l(jy));  
    params.yj_b(y,jy) = length(1919:years(y)-age_index_u(jy));
    end
end
%-----------------------------------------------------------------%

%%% DSS_MigrationEstimates_1991to2016_future_projections.m:
idx=9; %line 17
params.NOM_Series = 0; %0 = Moderate, 1 = Short return to pre covid, 2 = Longer return to pre covid.

params.NOM_tot = dlmread('Inputs/NOM_all_v2.csv',',',[1 (idx+params.NOM_Series*10) 100 (idx+params.NOM_Series*10)]); %line 20

params.NOM_age = xlsread('Inputs/NOM_Age_STT.xlsx',model_labels{idx},'A2:CV8');
%lines 81 to 84:


%% Read in Top 30 COB by age group (numbers)
% DSS 1991 - 2003
params.DSS_Top30Age = xlsread('Inputs/Top30COB_AgeDist_1991_to_2019_DSS_ABS.xlsx',model_labels{idx},'C2:I391'); %1991 to 2003 inclusive %10yr age group

params.Top30_Prev_high = dlmread('Inputs/Top30_COB_Prev.csv',',',[1 2 30 2]); %higher prev estimates for COBs (Cowie PhD)
params.regions_top30_DSS = dlmread('Inputs/Top30_COB_Prev.csv',',',[1 3 30 3]); %current prev estimates for COBs (Jenns sys rev, 2017) 

% ABS 2004 - 2019/2020
params.ABS_Top30Age2004to2021 = xlsread('Inputs/ABS_Top30_COB_5yrAgeDist_2004to2021.xlsx',model_labels{idx},'C3:T543'); % đã update, chuyển tên biến theo tên file, và khoảng datta input từ T512 -> T543 %2004 - 2021 %5yr age group(ABS2022 update)

params.Top30_Prev_high_from2004 = dlmread('Inputs/Top30_COB_Prev_from2004.csv',',',[1 2 30 2]); %higher prev estimates for COBs (Cowie PhD)
params.regions_top30_DSS_from2004 = dlmread('Inputs/Top30_COB_Prev_from2004.csv',',',[1 3 30 3]); %current prev estimates for COBs (Jenns sys rev, 2017) 

%% Other COB - read in total numbers by year 
% DSS 1991 - 2003
params.DSSbyCOB_other_1991to2003 = xlsread('Inputs/DSS_other_COB_numbers_1991to2003.xlsx',model_labels{idx},'B2:N222');
params.PrevDSSbyCOB_other_1991to2003 = dlmread('Inputs/DSS_other_COB_prev.csv',',',[1 1 221 1])/100; %MacLachlan 2017 Ests
params.PrevDSSbyCOB_other2019 = dlmread('Inputs/DSS_other_COB_prev.csv',',',[1 3 221 3])/100; %MacLachlan 2019 Ests
params.regions_other_DSS = dlmread('Inputs/DSS_other_COB_prev.csv',',',[1 2 221 2]);

% ABS 2004 - 2019/2020
params.ABSbyCOB_other_from2004 = xlsread('Inputs/ABS_other_COB_numbers_from2004.xlsx',model_labels{idx},'B2:S253'); % đã update file này từ R253 -> S253
params.ABS_Other_agedist = xlsread('Inputs/ABS_other_COB_AvAgeDist_from2004.xlsx',model_labels{idx},'B2:AV19'); % đã update input nhưng k cần update code %up to and incl 2050 - applies overall average age distribution for 2015 -2019 for all 'other' countries forward
params.PrevABSbyCOB_other_from2004 = dlmread('Inputs/ABS_other_COB_prev_from2004.csv',',',[1 1 252 1])/100; %MacLachlan 2017 Ests
params.PrevABSbyCOB_other2019 = dlmread('Inputs/ABS_other_COB_prev_from2004.csv',',',[1 3 252 3])/100; %MacLachlan 2019 Ests
params.regions_other_ABS = dlmread('Inputs/ABS_other_COB_prev_from2004.csv',',',[1 2 252 2]);


% lines 483 to 488:
params.phase_dist_Africa = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[2 1 8 4]);
params.phase_dist_Asia = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[12 1 18 4]);
params.phase_dist_combined = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[45 1 51 4]); % orig = Phase_a_gram_BC_estimates2.csv ? for combined only it differs. 
params.phase_dist_Africa_2004onwards = dlmread('Inputs/Phase_a_gram_BC_estimates2_18age_groups.csv',',',[2 1 19 4]);
params.phase_dist_Asia_2004onwards = dlmread('Inputs/Phase_a_gram_BC_estimates2_18age_groups.csv',',',[23 1 40 4]);
% Age group 60+ for Indigenous not known:
%phase_dist_Indigenous = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[32 1 37 4]);
params.phase_dist_combined_2004onwards = dlmread('Inputs/Phase_a_gram_BC_estimates2_18age_groups.csv',',',[65 1 82 4]); % orig = Phase_a_gram_BC_estimates2.csv ? for combined only it differs. 


%% AMALGAMATE PARAMETERS TO BE PASSED TO MODEL

%params = dis_params;

% Population parameters
names = fieldnames(pop_params);
for i =1:length(names)
    params.(names{i}) = pop_params.(names{i});
end

names = fieldnames(dis_params);
for i =1:length(names)
    params.(names{i}) = dis_params.(names{i});
end

% Other parameters
names = fieldnames(sim_params);
for i =1:length(names)
    params.(names{i}) = sim_params.(names{i});
end

% Surveillance/Intervention parameters
names = fieldnames(surv_params);
for i =1:length(names)
    params.(names{i}) = surv_params.(names{i});
end




%% Derive average prob of prog to Im Tol for age groups:    
    function [sigma_av] = sig_av(ages)
      % Age-dependent probability of progressing to Immune Tol from Acute:
        a=0:100; %1 by 101 vect. 
        sigma=zeros(1,101);
        for z=1:length(a)
            if a(z)< 1 %(a(i) == 0) || (a(i)==0.5)
                sigma(z) =0.885;
            else
                sigma(z) = exp(-0.645*a(z).^(0.455));
            end
        end  
        av_prob = zeros(1,length(ages));
      %sigma(i) = prob of prog to Im tol for ages = 0:0.5:100.
        for k = 1:length(ages)  
           if k==length(ages) %Last age group = 85+ = 85 to 100. 
             j=ages(k):100;
             av_prob(k) = sum( sigma(j+1) )/length(j);
           else
               j=ages(k):ages(k+1)-1;
               av_prob(k) = sum( sigma(j+1) )/length(j);
           end
        end
      sigma_av = av_prob;
    end



end



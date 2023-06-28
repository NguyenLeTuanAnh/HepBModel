%% Run a single simulation of the HBV2018 model:
model = HBV2021_struct_Nat();
tic
% define the parameter distributions
params = model.params();

[time, params, state, events] = single_simulation(model, params);
toc
%% 

%% 

% DEFINE STATES
% Disease states
S=state(:,1:params.nA);
V=state(:,params.nA+1:2*params.nA);
A=state(:,2*params.nA+1:3*params.nA);
I=state(:,3*params.nA+1:4*params.nA);
I1=state(:,4*params.nA+1:5*params.nA);
I2=state(:,5*params.nA+1:6*params.nA);
I3=state(:,6*params.nA+1:7*params.nA);
C1=state(:,7*params.nA+1:8*params.nA);
C2=state(:,8*params.nA+1:9*params.nA);
C3=state(:,9*params.nA+1:10*params.nA);
I1T=state(:,10*params.nA+1:11*params.nA);
I2T=state(:,11*params.nA+1:12*params.nA);
I3T=state(:,12*params.nA+1:13*params.nA);
C1T=state(:,13*params.nA+1:14*params.nA);
C2T=state(:,14*params.nA+1:15*params.nA);
C3T=state(:,15*params.nA+1:16*params.nA);
H=state(:,16*params.nA+1:17*params.nA);
HT=state(:,17*params.nA+1:18*params.nA);
D=state(:,18*params.nA+1:19*params.nA);
DT=state(:,19*params.nA+1:20*params.nA);
R=state(:,20*params.nA+1:21*params.nA);

Ac_incidence =state(:,21*params.nA+1:22*params.nA); %ac incidence
ImTol_incidence =state(:,22*params.nA+1:23*params.nA); % # new chronic inf
HCC_deaths = state(:,23*params.nA+1:24*params.nA); %deaths HCC
HCC_deaths_T= state(:,24*params.nA+1:25*params.nA); % deaths HCC on T
DC_deaths = state(:,25*params.nA+1:26*params.nA); % deaths DC
DC_deaths_T = state(:,26*params.nA+1:27*params.nA); % deaths DC T
bg_deaths = state(:,27*params.nA+1:28*params.nA); %bg deaths in people living with CHB 
acute_deaths = state(:,28*params.nA+1:29*params.nA); %acute hbv deaths
HCC_incidence_cum = state(:,29*params.nA+1:30*params.nA);
HCC_T_incidence_cum = state(:,30*params.nA+1:31*params.nA);
DC_incidence_cum = state(:,31*params.nA+1:32*params.nA);
DC_T_incidence_cum = state(:,32*params.nA+1:33*params.nA);
Treat_uptake_cum = state(:,33*params.nA+1:34*params.nA);
 

%% Total in infection state by age: 
Total_Chronic_no_T = I + I1 + I2 + I3 + C1 + C2 + C3 + H + D;
Total_Chronic_T = I1T + I2T + I3T + C1T + C2T + C3T + HT + DT;
Total_Chronic = Total_Chronic_no_T + Total_Chronic_T;
%% 

SumTotal_Chronic = sum(Total_Chronic,2);
SumTotal_Chronic_onT = sum(Total_Chronic_T');
SumTotal_Chronic_no_T = sum(Total_Chronic_no_T');
%SumTotal_Chronic_report = round(SumTotal_Chronic)


%% Cirrhotic vs Non-Cirrhotic by Age
Cirrhotic_age = round(C1+C2+C3+C1T+C2T+C3T + D+DT +H+HT); %10 states
NonCirrhotic_age = round(I+I1+I2+I3+I1T+I2T+I3T); % %6 states
Cirrhotic_total = sum(Cirrhotic_age');
NonCirrhotic_total = sum(NonCirrhotic_age');

Phase_Cirr_Report = round([ sum(I')', sum(I1'+I1T')', sum(I2'+I2T')', sum(I3'+I3T')', sum(C1'+C1T')', sum(C2'+C2T')', sum(C3'+C3T')']);
Phase_Report = round([ sum(I')', sum(I1'+I1T')', sum(I2'+I2T')', sum(I3'+I3T')', sum(C1'+C1T')', sum(C2'+C2T')', sum(C3'+C3T')',(sum(H') + sum(HT'))',(sum(D') + sum(DT'))']);

% Mortality by age
HCC_deaths_age_all = HCC_deaths + HCC_deaths_T;
HCC_incident_deaths_age = round([zeros(1,18); diff(HCC_deaths_age_all)]);
DC_deaths_age_all = DC_deaths + DC_deaths_T;
DC_incident_deaths_age = round([zeros(1,18); diff(DC_deaths_age_all)]);

% Cumulative number of deaths
Cumulative_HCC_deaths_T = sum(HCC_deaths_T')';
Cumulative_HCC_deaths = sum(HCC_deaths')';
Cumulative_HCC_deaths_all = Cumulative_HCC_deaths_T + Cumulative_HCC_deaths;

% Cumulative number of deaths
Cumulative_DC_deaths_T = sum(DC_deaths_T')';
Cumulative_DC_deaths = sum(DC_deaths')';
Cumulative_DC_deaths_all = Cumulative_DC_deaths_T + Cumulative_DC_deaths;


%% Number of deaths per year
incident_HCC_deaths_T = [0; diff(Cumulative_HCC_deaths_T)];
incident_HCC_deaths= [0; diff(Cumulative_HCC_deaths)];

incident_DC_deaths_T = [0, sum(diff(DC_deaths_T)')];
incident_DC_deaths = [0, sum(diff(DC_deaths)')];

incident_bg_deaths =[0, sum(diff(bg_deaths)')];

%% Total mortality by year
incident_DC_deaths_all = [0; diff(Cumulative_DC_deaths_all)];
incident_HCC_deaths_all = [0; diff(Cumulative_HCC_deaths_all)];
HBV_incident_deaths_all = incident_DC_deaths_all + incident_HCC_deaths_all; 


%% Incidence of HCC and DC each year
% Incident number of HCC cases
HCC_incidence_age = [HCC_incidence_cum(1,:); diff(HCC_incidence_cum)];
HCC_incidence_onT = sum(HCC_incidence_age');
HCC_T_incidence_age = [HCC_T_incidence_cum(1,:); diff(HCC_T_incidence_cum)];
HCC_T_incidence_all = sum(HCC_T_incidence_age');

HCC_incidence_total = HCC_incidence_onT + HCC_T_incidence_all;

% Incident number of DC cases
DC_incidence_age = [DC_incidence_cum(1,:); diff(DC_incidence_cum)];
DC_incidence_age_T = [DC_T_incidence_cum(1,:); diff(DC_T_incidence_cum)];
DC_incidence_noT = sum(DC_incidence_age');
DC_incidence_T = sum(DC_incidence_age_T'); 

DC_incidence_total = DC_incidence_noT + DC_incidence_T;


% Cumulative number of acute deaths:
Cumulative_acute_deaths = sum(acute_deaths')';
% Number of acute HBV deaths per year
Acute_deaths_overtime = [0; diff(Cumulative_acute_deaths)];
% non-HBV related deaths
Cumulative_bg_deaths = sum(bg_deaths')';
bg_deaths_overtime = [0; diff(Cumulative_bg_deaths)];

Total_Pop_Age = (S+V+A+I+I1+I2+I3+C1+C2+C3+I1T+I2T+I3T+C1T+C2T+C3T+H+HT+D+DT+R);
SumTotal_Pop= sum(Total_Pop_Age,2);

% Total HBsAg prevalence (prop of pop chronicly infected with HBV)
% HBsAg = HBV surface antigen
Total_sAg_prev = SumTotal_Chronic./SumTotal_Pop;
%size(Total_sAg_prev)
sAg_prev_age= Total_Chronic./Total_Pop_Age;

%% Number of people ever chronicly infected up to each year:
%Years1(21:66) = 1971 to 2016
% Prop diag including deaths (used in 2016 & BC thesis):
Ever_Chronic = SumTotal_Chronic + sum((HCC_deaths + HCC_deaths_T)')' + sum((DC_deaths + DC_deaths_T)')' + sum(bg_deaths')';
% Same as Ever_Chronic1: Ever_Chronic3 =  SumTotal_Chronic' + Cumulative_DC_deaths_all  + Cumulative_HCC_deaths_all;

EvChronic_Est = cumsum(SumTotal_Chronic) - cumsum( sum((HCC_deaths + HCC_deaths_T)')' + sum((DC_deaths + DC_deaths_T)')' + sum(bg_deaths')');
Ever_Chronic_report = round(Ever_Chronic); 

%% Proportion Diagnosed for years 1971 to 2016
%params.notifs1971to2017; %46 by 1
prop_diag = params.cum_notifs1971to2021./Ever_Chronic(21:71); % 68% in 2016
prop_not_diag = 1 - prop_diag;

%% save vao file "National2022_futuretx_to2030_report.mat"

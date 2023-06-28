%% Determine the age- and time- specific migration numbers
%  each year entering the population for HBV2018 Model.
% THEN determine which infection state they enter in.
function [MIG_S, MIG_I,MIG_I1,MIG_I2, MIG_I3, MIG_R] = MigrationEstimates_1951to1990(params)
% Migration Age Distribution is in 10 year blocks up to 60+
% 0 to 9; 10 to 19; 20 to 29; 30 to 39; 40 to 49; 50 to 59; 60+ 
%years =1951:2015; %65 years

idx = 9; % 1 - ACT; 2- NSW; 3-NT; 4-QLD; 5-SA; 6-TAS; 7-VIC; 8-WA; 9-AUS
%% National/STT NOM by age = age distribution of ALL migrants entering population
% Total NOM into Australia from 1951 to 2050: (incl ABS Series B projections)
NOM_tot = params.NOM_tot;%dlmread('Inputs/NOM_all.csv',',',[1 idx 99 idx]);

% NOM_agedist(age, year)
% ages: 0-4;5-9;10-14;15-19;20-24;25-29;30-34;35-39;40-44;45-49;50-54;
% 55-59; 60-64; 65+
NOM_agedist = zeros(14,100); 
NOM_agedist(:,54:65) = dlmread('Inputs/NOM_AgeDist_2004to2015_clean.csv',',',[1 1 12 14])'; %years(54:65) = 2004:2015
%sum(NOM_agedist) should = 1!


% Age distribution derived from Victorian migrants 1975/76 to 2006/07 
BC_age_dist = [11, 9, 7, 7, 11, 15, 13, 9, 5, 3, 3, 2, 2, 3]'./100;


% Apply the historical Vic migrant age dist (or 2004 NOM age distribution)
% to 1951-2003:
for i=1:53
NOM_agedist(:,i) = BC_age_dist;% NOM_agedist(:,54);
end

% Take the average age distribution of 2014-2015 and apply it forward
for i=66:100
    for j=1:14
    NOM_agedist(j,i) = mean(NOM_agedist(j,64:65));
    end
end


%bar(NOM_agedist(:,54:100))
NOM_age=zeros(1,length(NOM_agedist)); 
for i=1:2:13 
NOM_age=[NOM_age; NOM_agedist(i,:)+NOM_agedist(i+1,:)];
end
NOM_age=NOM_age(2:end,:);


%--------------------------------------------------------------------------
%% 1951 to 1974
%% Federation to Centurys End Data to get Permanent Settlers by COB
% rows = countries; cols = years 1951:1974 
Arr_COB_prop = dlmread('Inputs/PermArrCOB_1951to1974.csv',',',[2 1 39 24])/100;
Arr_COB_tot = zeros(size(Arr_COB_prop));
for i = 1:24 
Arr_COB_tot(:,i) = Arr_COB_prop(:,i)*NOM_tot(i);
end

% Arr_China_1951to1974=Arr_COB_tot(1,:); %China
% Arr_Vietnam_1951to1974=Arr_COB_tot(2,:); %Vietnam
% Arr_Philippines_1951to1974=Arr_COB_tot(3,:); %Philippines

Arr_Others_1951to1974 = Arr_COB_tot(4:end,:); %Remaining Countries

% Prevalence Est (excl China, Viet & Philippines)
%prev_COB_1951to1974 = dlmread('Inputs/Prev_ArrCOB_1951to1974.csv',',',[1 1 35 1])/100;
prev_COB_1951to1974 = params.prev_COB_1951to1974; %dlmread('Inputs/Prev_ArrCOB_1951to1974_BC.csv',',',[1 1 35 1]);
region_COB_1951to1974 = dlmread('Inputs/Prev_ArrCOB_1951to1974.csv',',',[1 2 35 2]);


% Determine propoprtions resolved and susceptible (assuming no vaccination)
r = zeros(length(prev_COB_1951to1974),1);
s_1951to1974 = zeros(length(prev_COB_1951to1974),1);
for k=1:length(prev_COB_1951to1974)
%k
r(k) = 7*prev_COB_1951to1974(k); %7 was the ratio of prop resolved to prop chronic from Cowie2009 serosurveys.     

if r(k)>0.75
    r(k)=0.75;
end
%r(k)
s_1951to1974(k) = 1 - prev_COB_1951to1974(k) - r(k);
end

%countries by years
inf_other_1951to1974 = zeros(35 , length(1951:1974));
sus_other_1951to1974 = zeros(35 , length(1951:1974));
res_other_1951to1974 = zeros(35 , length(1951:1974));
for i = 1:length(1951:1974) %24
    inf_other_1951to1974(:,i) = Arr_Others_1951to1974(:,i).*prev_COB_1951to1974;
    sus_other_1951to1974(:,i) = Arr_Others_1951to1974(:,i).*(s_1951to1974);
    res_other_1951to1974(:,i) = Arr_Others_1951to1974(:,i).*(r);
end


% Derive number of chronic infections from each region per year
% to apply phase distribution to later:
s=size(inf_other_1951to1974);
I_Africa1951 = zeros(1,s(2));
I_Asia1951 = zeros(1,s(2));
I_Other1951 = zeros(1,s(2));

for g=1:s(1)
if region_COB_1951to1974(g)==1
    I_Africa1951(1,:) = I_Africa1951(1,:) + inf_other_1951to1974(g,:);
elseif region_COB_1951to1974(g)==2
    I_Asia1951(1,:) = I_Asia1951(1,:) + inf_other_1951to1974(g,:);
else
    I_Other1951(1,:) = I_Other1951(1,:) + inf_other_1951to1974(g,:);
end
end
%I_Africa+I_Asia+I_Other %should = sum(inf_other_1951to1974)


% row 1 = S, row 2 = I, row 3 = Res
total_other_mig1951to1974 = [sum(sus_other_1951to1974); sum(inf_other_1951to1974); sum(res_other_1951to1974)];

%--------------------------------------------------------------------------
%% 1975 to 1990
%% 3412.0 Migration Data from ABS by COB
% rows = countries; cols = years 1951:1974 
MigCOB_tot = dlmread('Inputs/MigCOB_1975to1990.csv',',',[4 1 72 16]);

% MigCOB_tot_China_1975to1990=MigCOB_tot(1,:); %China
% MigCOB_tot_Vietnam_1975to1990=MigCOB_tot(2,:); %Vietnam
% MigCOB_tot_Philippines_1975to1990=MigCOB_tot(3,:); %Philippines

MigCOB_tot_Others_1975to1990 = MigCOB_tot(5:end,:); %Remaining Countries

% Prevalence Est (excl China, Viet & Philippines)
%prev_COB_1975to1990 = dlmread('Inputs/Prev_MigCOB_1975to1990.csv',',',[1 1 66 1])/100;
prev_COB_1975to1990 = params.prev_COB_1975to1990; %dlmread('Inputs/Prev_MigCOB_1975to1990_BC.csv',',',[1 1 66 1]);
region_COB_1975to1990 = dlmread('Inputs/Prev_MigCOB_1975to1990.csv',',',[1 2 65 2]);

r1 = zeros(length(prev_COB_1975to1990),1);
s_1975to1990 = zeros(length(prev_COB_1975to1990),1);
for k=1:length(prev_COB_1975to1990)

r1(k) = 7*prev_COB_1975to1990(k); %7 was the ratio of prop resolved to prop chronic from Cowie2009 serosurveys.     

if r1(k)>0.75
    r1(k)=0.75;
end
%r(k)
s_1975to1990(k) = 1 - prev_COB_1975to1990(k) - r1(k);

end


% Countries by years
inf_other_1975to1990 = zeros(65 , length(1975:1990));
sus_other_1975to1990 = zeros(65, length(1975:1990));
res_other_1975to1990 = zeros(65 , length(1975:1990));
for i = 1:length(1975:1990) %16
inf_other_1975to1990(:,i) = MigCOB_tot_Others_1975to1990(:,i).*prev_COB_1975to1990;
sus_other_1975to1990(:,i) = MigCOB_tot_Others_1975to1990(:,i).*(s_1975to1990);
res_other_1975to1990(:,i) = MigCOB_tot_Others_1975to1990(:,i).*(r1);
end

% Derive number of chronic infections from each region per year
% to apply phase distribution to later:
ss=size(inf_other_1975to1990);
I_Africa1975 = zeros(1,ss(2));
I_Asia1975 = zeros(1,ss(2));
I_Other1975 = zeros(1,ss(2));

for g=1:ss(1)
if region_COB_1975to1990(g)==1
    I_Africa1975(1,:) = I_Africa1975(1,:) + inf_other_1975to1990(g,:);
elseif region_COB_1975to1990(g)==2
    I_Asia1975(1,:) = I_Asia1975(1,:) + inf_other_1975to1990(g,:);
else
    I_Other1975(1,:) = I_Other1975(1,:) + inf_other_1975to1990(g,:);
end
end
%I_Africa1975 + I_Asia1975 + I_Other1975 = sum(inf_other_1975to1990)


% row 1 = S, row 2 = I, row 3 = Res
total_other_mig1975to1990 = [sum(sus_other_1975to1990); sum(inf_other_1975to1990); sum(res_other_1975to1990)];


%% Gather totals from 1951 to 1974 
I_Africa1951to1990 = [I_Africa1951, I_Africa1975];
I_Asia1951to1990 = [I_Asia1951, I_Asia1975];
I_Other1951to1990 = [I_Other1951, I_Other1975];

%I_Africa1951to1990 + I_Asia1951to1990 + I_Other1951to1990

%% Total NOM coming in states Sus, I (chronic) and Resolved 1951 to 2016.
total_other_mig_SIR_1951to1990 = [total_other_mig1951to1974, total_other_mig1975to1990];
%Sus, Chronic & Res migrants by age:
%Age distribution: 7 by 100: NOM_age
other_mig_S = zeros(7,40);
other_mig_I = zeros(7,40);
other_mig_R = zeros(7,40);
I_Africa1951to1990_age = zeros(7,40);
I_Asia1951to1990_age = zeros(7,40);
I_Other1951to1990_age = zeros(7,40);

for i=1:40
other_mig_S(:,i)= total_other_mig_SIR_1951to1990(1,i)*NOM_age(:,i);
other_mig_I(:,i)= total_other_mig_SIR_1951to1990(2,i)*NOM_age(:,i);

I_Africa1951to1990_age(:,i) = I_Africa1951to1990(i)*NOM_age(:,i); 
I_Asia1951to1990_age(:,i) = I_Asia1951to1990(i)*NOM_age(:,i); 
I_Other1951to1990_age(:,i) = I_Other1951to1990(i)*NOM_age(:,i); 

other_mig_R(:,i)= total_other_mig_SIR_1951to1990(3,i)*NOM_age(:,i);
end

%% ASSIGN DISEASE PHASES by proportion to get

%% June 12th update uses BC_estimates2.csv (prev = BC_estimates.csv)

% Africa, Asia, US, Indigenous, Combined.
% Ages 0 to 9; 10 to 19; 20 to 29; 30 to 39; 40 to 49; 50 to 59; 60+ 
phase_dist_Africa = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[2 1 8 4]);
phase_dist_Asia = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[12 1 18 4]);
%phase_dist_US = dlmread('Inputs/Phase_a_gram_BC_estimates.csv',',',[22 1 28 4])/100;

% Age group 60+ for Indigenous not known:
%phase_dist_Indigenous = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[32 1 37 4]);
phase_dist_combined = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[45 1 51 4]);

I_1951to1990=zeros(7,length(1951:1990));
I1_1951to1990=zeros(7,length(1951:1990));
I2_1951to1990=zeros(7,length(1951:1990));
I3_1951to1990=zeros(7,length(1951:1990));

for i=1:length(1951:1990)
I_1951to1990(:,i) = (I_Africa1951to1990_age(:,i).*phase_dist_Africa(:,1) + I_Asia1951to1990_age(:,i).*phase_dist_Asia(:,1) + I_Other1951to1990_age(:,i).*phase_dist_combined(:,1));
I1_1951to1990(:,i) = (I_Africa1951to1990_age(:,i).*phase_dist_Africa(:,2) + I_Asia1951to1990_age(:,i).*phase_dist_Asia(:,2) + I_Other1951to1990_age(:,i).*phase_dist_combined(:,2));
I2_1951to1990(:,i) = (I_Africa1951to1990_age(:,i).*phase_dist_Africa(:,3) + I_Asia1951to1990_age(:,i).*phase_dist_Asia(:,3) + I_Other1951to1990_age(:,i).*phase_dist_combined(:,3));
I3_1951to1990(:,i) = (I_Africa1951to1990_age(:,i).*phase_dist_Africa(:,4) + I_Asia1951to1990_age(:,i).*phase_dist_Asia(:,4) + I_Other1951to1990_age(:,i).*phase_dist_combined(:,4));
end 


%% Need to derive migration for Top 3 for 1951 to 1990 then add to above. 

% Age Distribution of Migrants 1991 to 2016:
% 0 to 9; 10 to 19; 20 to 29; 30 to 39; 40 to 49; 50 to 59; 60+ 

%% CHINA, VIETNAM & PHILIPPINES ESTIMATES (top 3 COB with hep B) AND TAIWAN %%
% Total Numbers Migrating over time:
%1951 to 1974
Arr_China_1951to1974=Arr_COB_tot(1,:); %China
Arr_Vietnam_1951to1974=Arr_COB_tot(2,:); %Vietnam
Arr_Philippines_1951to1974=Arr_COB_tot(3,:); %Philippines
Arr_Taiwan_1951to1974= zeros(1,24); %Taiwan not listed as COB during this time for this data.

%1975 to 1990
MigCOB_tot_China_1975to1990=MigCOB_tot(1,:); %China
MigCOB_tot_Vietnam_1975to1990=MigCOB_tot(2,:); %Vietnam
MigCOB_tot_Philippines_1975to1990=MigCOB_tot(3,:); %Philippines
MigCOB_tot_Taiwan_1975to1990=MigCOB_tot(4,:); %Taiwan

Mig_tot_China_1951to1990 = [Arr_China_1951to1974, MigCOB_tot_China_1975to1990];
Mig_tot_Vietnam_1951to1990 = [Arr_Vietnam_1951to1974, MigCOB_tot_Vietnam_1975to1990];
Mig_tot_Phil_1951to1990 = [Arr_Philippines_1951to1974, MigCOB_tot_Philippines_1975to1990];
Mig_tot_Taiwan_1951to1990 = [Arr_Taiwan_1951to1974, MigCOB_tot_Taiwan_1975to1990]; %Taiwan


% Assign Age Distribution:
% Apply NOM age distribution to China, Viet & Philippines for 1951 - 1990:
Mig_age_China_1951to1990 = zeros(7 ,length(1951:1990));
Mig_age_Vietnam_1951to1990= zeros(7 ,length(1951:1990));
Mig_age_Phil_1951to1990= zeros(7 ,length(1951:1990));
Mig_age_Taiwan_1951to1990 = zeros(7 ,length(1951:1990));

%NOM_age(:,1:40) %7 by t
for i = 1:length(1951:1990)
Mig_age_China_1951to1990(:,i) = NOM_age(:,i)*Mig_tot_China_1951to1990(i);
Mig_age_Vietnam_1951to1990(:,i) = NOM_age(:,i)*Mig_tot_Vietnam_1951to1990(i);
Mig_age_Phil_1951to1990(:,i) = NOM_age(:,i)*Mig_tot_Phil_1951to1990(i);
Mig_age_Taiwan_1951to1990(:,i) = NOM_age(:,i)*Mig_tot_Taiwan_1951to1990(i);
end


% Determine numbers Chronicly infected (using prev data), Susceptible &
% Resolved:
% age_index_l = [4,14,29,39,49,59,60];
% age_index_u = [0,5,15,30,40,50,60];
age_index_l = [9:10:59, 60];
age_index_u = (0:10:60);
years = 1951:1990;

% China:
% Razavi: est 6.1% overall and 0.2% in those <5
% Jenn/ Sys Rev Estimates: 7.56%
% Serosurveys: Cui et al. 2017

% China: 
I_age_China = zeros(7, length(1951:1990));
R_age_China = zeros(7, length(1951:1990));
S_age_China = zeros(7, length(1951:1990));
prev_China=zeros(1,length(1891:1990));
prev_China(:,88:100) = [0.1017, 0.0983, 0.095, 0.094, 0.093, 0.092, 0.091, 0.09, 0.08, 0.075, 0.07, 0.06, 0.05];
for i = 1:length(1891:1977)
prev_China(i) = 10.5/100;
end

for y= 1:length(years)
        for j=1:length(age_index_l)

            a=length(1891:years(y)-age_index_l(j));
            b=length(1891:years(y)-age_index_u(j));

            %p(j) = mean( prev_top30_yob(z, a:b ) ); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            p_China(j) = (sum( prev_China(a:b) )/length(a:b)); 
     
res_china = 7* p_China;
% assuming no vaccination
sus_china = 1 - p_China - res_china;
I_age_China= Mig_age_China_1951to1990(:,y).*p_China;
R_age_China= Mig_age_China_1951to1990(:,y).*res_china;
S_age_China= Mig_age_China_1951to1990(:,y).*sus_china;
        end
end

for y= 1:length(years)
        for j=1:length(age_index_l)

            a=length(1891:years(y)-age_index_l(j));
            b=length(1891:years(y)-age_index_u(j));

            %p(j) = mean( prev_top30_yob(z, a:b ) ); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            p_China(j) = (sum( prev_China(a:b) )/length(a:b)); 
     
res_china = 7* p_China';
% assuming no vaccination
sus_china = 1 - p_China' - res_china;
I_age_China(:,y)= Mig_age_China_1951to1990(:,y).*p_China';
R_age_China(:,y)= Mig_age_China_1951to1990(:,y).*res_china;
S_age_China(:,y)= Mig_age_China_1951to1990(:,y).*sus_china;
        end
end

% Vietnam: vacc in 1997

% Note: 0 migrants recorded from Vietnam 1951 to 1974
I_age_Viet = zeros(7,length(1951:1990));
R_age_Viet = zeros(7,length(1951:1990));
S_age_Viet = zeros(7,length(1951:1990));
% Moi 330 - 384
prev_Viet=zeros(1,length(1891:1990));

for i=1:length(1891:1951)
    prev_Viet(i) = 0.12; %prevalence in birth year in VN
end

for i=length(1891:1952):length(1891:1970)
    prev_Viet(i) = prev_Viet(i-1)-((0.12-0.1)/(1970-1952+1));
end

for i=length(1891:1971):length(1891:1990)
    prev_Viet(i) = prev_Viet(i-1)-((0.1-0.085)/(1990-1971+1));
end

% p_Viet = zeros(length(age_index_l),1);

for y= 1:length(years)
        for j=1:length(age_index_l)

            a=length(1891:years(y)-age_index_l(j));
            b=length(1891:years(y)-age_index_u(j));

            %p(j) = mean( prev_top30_yob(z, a:b ) ); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            p_Viet(j) = (sum( prev_Viet(a:b) )/length(a:b)); 
     
res_viet = 7* p_Viet;
% assuming no vaccination
sus_viet = 1 - p_Viet - res_viet;
I_age_Viet= Mig_age_Vietnam_1951to1990(:,y).*p_Viet;
R_age_Viet= Mig_age_Vietnam_1951to1990(:,y).*res_viet;
S_age_Viet= Mig_age_Vietnam_1951to1990(:,y).*sus_viet;
        end
end

for y= 1:length(years)
        for j=1:length(age_index_l)

            a=length(1891:years(y)-age_index_l(j));
            b=length(1891:years(y)-age_index_u(j));

            %p(j) = mean( prev_top30_yob(z, a:b ) ); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            p_Viet(j) = (sum( prev_Viet(a:b) )/length(a:b)); 
     
res_viet = 7* p_Viet';
% assuming no vaccination
sus_viet = 1 - p_Viet' - res_viet;
I_age_Viet(:,y)= Mig_age_Vietnam_1951to1990(:,y).*p_Viet';
R_age_Viet(:,y)= Mig_age_Vietnam_1951to1990(:,y).*res_viet;
S_age_Viet(:,y)= Mig_age_Vietnam_1951to1990(:,y).*sus_viet;
        end
end

% Philippines: little data.
% Razavi: est 9.8% overall and 2.2% in those <5
% Jenn/ Sys Rev Estimates: 3.02 (highest ~ 4.5%)
% Assume prior to 2005 prev = 4.5%. From 2005 onwards = 3% with 2% <5?
prev_phil_high=9.8/100; %4.5/100;
res_phil = 7*prev_phil_high;
sus_phil = 1 - res_phil -prev_phil_high;

I_age_Phil = zeros(7, length(1951:1990));
R_age_Phil = zeros(7, length(1951:1990));
S_age_Phil = zeros(7, length(1951:1990));

for i=1:length(1951:1990)
I_age_Phil(:,i) = Mig_age_Phil_1951to1990(:,i).*prev_phil_high;
R_age_Phil(:,i) = Mig_age_Phil_1951to1990(:,i).*res_phil;
S_age_Phil(:,i) = Mig_age_Phil_1951to1990(:,i).*sus_phil;
end

% Taiwan

I_age_Taiwan = zeros(7, length(1951:1990));
R_age_Taiwan = zeros(7, length(1951:1990));
S_age_Taiwan = zeros(7, length(1951:1990));
prev_Taiwan=zeros(1,length(1891:1990));
prev_Taiwan(:,73:100) = [14.8	15	14.8	15.3	15	15.3	15	14.7	14.5	14.3	14	13.3	13.3	12.8	12.4	11.8	11.6	11.2	10.1	9.2	8.1	7	5.7	4.5	2.8	2.8	2.3	2.3	]./100;
for i = 1:length(1891:1962)
prev_Taiwan(i) = 14.8/100;
end

for y= 1:length(years)
        for j=1:length(age_index_l)

            a=length(1891:years(y)-age_index_l(j));
            b=length(1891:years(y)-age_index_u(j));

            %p(j) = mean( prev_top30_yob(z, a:b ) ); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            p_Taiwan(j) = (sum( prev_Taiwan(a:b) )/length(a:b)); 
     
res_taiwan = 7* p_Taiwan;
% assuming no vaccination
sus_taiwan = 1 - p_Taiwan - res_taiwan;
I_age_Taiwan= Mig_age_Taiwan_1951to1990(:,y).*p_Taiwan;
R_age_Taiwan= Mig_age_Taiwan_1951to1990(:,y).*res_taiwan;
S_age_Taiwan= Mig_age_Taiwan_1951to1990(:,y).*sus_taiwan;
        end
end

for y= 1:length(years)
        for j=1:length(age_index_l)

            a=length(1891:years(y)-age_index_l(j));
            b=length(1891:years(y)-age_index_u(j));

            %p(j) = mean( prev_top30_yob(z, a:b ) ); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            p_Taiwan(j) = (sum( prev_Taiwan(a:b) )/length(a:b));  
     
res_taiwan = 7* p_Taiwan';
% assuming no vaccination
sus_taiwan = 1 - p_Taiwan' - res_taiwan;
I_age_Taiwan(:,y)= Mig_age_Taiwan_1951to1990(:,y).*p_Taiwan';
R_age_Taiwan(:,y)= Mig_age_Taiwan_1951to1990(:,y).*res_taiwan;
S_age_Taiwan(:,y)= Mig_age_Taiwan_1951to1990(:,y).*sus_taiwan;
        end
end


china_migI =zeros(7,length(1951:1990));
china_migI1 =zeros(7,length(1951:1990));
china_migI2 =zeros(7,length(1951:1990));
china_migI3 =zeros(7,length(1951:1990));

phil_migI =zeros(7,length(1951:1990));
phil_migI1 =zeros(7,length(1951:1990));
phil_migI2 =zeros(7,length(1951:1990));
phil_migI3 =zeros(7,length(1951:1990));

viet_migI =zeros(7,length(1951:1990));
viet_migI1 =zeros(7,length(1951:1990));
viet_migI2 =zeros(7,length(1951:1990));
viet_migI3 =zeros(7,length(1951:1990));

Taiwan_migI =zeros(7,length(1951:1990));
Taiwan_migI1 =zeros(7,length(1951:1990));
Taiwan_migI2 =zeros(7,length(1951:1990));
Taiwan_migI3 =zeros(7,length(1951:1990));

% Age Distribution of phases 1951 to 1990:
for i=1:length(1951:1990)
%
china_migI(:,i) =I_age_China(:,i).*phase_dist_Asia(:,1);
china_migI1(:,i)=I_age_China(:,i).*phase_dist_Asia(:,2);
china_migI2(:,i)=I_age_China(:,i).*phase_dist_Asia(:,3);
china_migI3(:,i)=I_age_China(:,i).*phase_dist_Asia(:,4);
%
phil_migI(:,i) =I_age_Phil(:,i).*phase_dist_Asia(:,1);
phil_migI1(:,i)=I_age_Phil(:,i).*phase_dist_Asia(:,2);
phil_migI2(:,i)=I_age_Phil(:,i).*phase_dist_Asia(:,3);
phil_migI3(:,i)=I_age_Phil(:,i).*phase_dist_Asia(:,4);
%
viet_migI(:,i) =I_age_Viet(:,i).*phase_dist_Asia(:,1);
viet_migI1(:,i)=I_age_Viet(:,i).*phase_dist_Asia(:,2);
viet_migI2(:,i)=I_age_Viet(:,i).*phase_dist_Asia(:,3);
viet_migI3(:,i)=I_age_Viet(:,i).*phase_dist_Asia(:,4);
%
Taiwan_migI(:,i) =I_age_Taiwan(:,i).*phase_dist_Asia(:,1);
Taiwan_migI1(:,i)=I_age_Taiwan(:,i).*phase_dist_Asia(:,2);
Taiwan_migI2(:,i)=I_age_Taiwan(:,i).*phase_dist_Asia(:,3);
Taiwan_migI3(:,i)=I_age_Taiwan(:,i).*phase_dist_Asia(:,4);
end


%-------------------------------------------------------------------------------
%% Combine to get total migration input for model:
%-------------------------------------------------------------------------------
% Will need to multiple by total settlers/NOM ratio:
 MIG_S  = S_age_China  + S_age_Viet  + S_age_Phil + S_age_Taiwan  + other_mig_S; 
 MIG_I  = china_migI  + viet_migI  + phil_migI  + Taiwan_migI  + I_1951to1990; 
 MIG_I1 = china_migI1 + viet_migI1 + phil_migI1 + Taiwan_migI1   + I1_1951to1990;
 MIG_I2 = china_migI2 + viet_migI2 + phil_migI2 + Taiwan_migI2   + I2_1951to1990; 
 MIG_I3 = china_migI3 + viet_migI3 + phil_migI3 + Taiwan_migI3   + I3_1951to1990; 
 MIG_R  = R_age_China  + R_age_Viet  + R_age_Phil + R_age_Taiwan  + other_mig_R; 

 MIG_Chron = MIG_I + MIG_I1 + MIG_I2 + MIG_I3;
 

%% Average over refined age-groups for model input:
% NOM Age Groups: 0-4;5-9;10-14;15-19;20-24;25-29;30-34;35-39;40-44;45-49;
% 50-54;55-59;60-64;65-69;70-74;75-79;80-84;85+
newNOM_agedist = [NOM_agedist(:,1:length(1951:1990)); zeros(4,length(1951:1990))];
%NOM_agedist(14,1) - divide into 65-69;70-74;75-79;80-84;85+
ad = [1.5; 1; 0.4; 0.05; 0.05]./100; 

for i=1:length(1951:1990)
newNOM_agedist(14:18,i) = ad;
end

MIG_S = newNOM_agedist(:,1:40)*diag(sum(MIG_S));
MIG_I = newNOM_agedist(:,1:40)*diag(sum(MIG_I)); 
MIG_I1 = newNOM_agedist(:,1:40)*diag(sum(MIG_I1));
MIG_I2 = newNOM_agedist(:,1:40)*diag(sum(MIG_I2));
MIG_I3 = newNOM_agedist(:,1:40)*diag(sum(MIG_I3));
MIG_R = newNOM_agedist(:,1:40)*diag(sum(MIG_R));

total_mig_ratio = (sum(MIG_S + MIG_I + MIG_I1 + MIG_I2 + MIG_I3 + MIG_R)'./NOM_tot(1:40))';
% 
for i=1:18
    MIG_S(i,:) = MIG_S(i,:)./total_mig_ratio;
    MIG_I(i,:) = MIG_I(i,:)./total_mig_ratio; 
    MIG_I1(i,:)= MIG_I1(i,:)./total_mig_ratio;
    MIG_I2(i,:) = MIG_I2(i,:)./total_mig_ratio; 
    MIG_I3(i,:) = MIG_I3(i,:)./total_mig_ratio; 
    MIG_R(i,:) = MIG_R(i,:)./total_mig_ratio;
end

MIG_S = round(MIG_S);
MIG_I = round(MIG_I);
MIG_I1  = round(MIG_I1);
MIG_I2  = round(MIG_I2);
MIG_I3 = round(MIG_I3);
MIG_R = round(MIG_R);

MIG_Chron_a = MIG_I + MIG_I1 + MIG_I2 + MIG_I3;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATED APRIL 2021 - generate total number of CHB by age migrating 1991 -
% 2003 for input into HBV model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mig_S, mig_I, mig_I1, mig_I2, mig_I3, mig_R] = ABS_MigrationEstimates_1991to2003(params)

y_end = 2003; %2050;
years=1991:y_end;

% 1 - ACT; 2- NSW; 3-NT; 4-QLD; 5-SA; 6-TAS; 7-VIC; 8-WA; 9-AUS
idx=9;

% Total NOM into Australia from 1951 to 2050: (incl ABS Series B projections)
NOM_tot = params.NOM_tot; %dlmread('Inputs/NOM_STT_all.csv',',',[1 idx 100 idx]);

%---------------------------------------------------------------------------

%NOM_agedist = dlmread('Inputs/NOM_agedist_DSSMig.csv',',',[0 0 13 99]); 
NOM_age = params.NOM_age; % dlmread('Inputs/NOM_age_DSSMig.csv',',',[0 0 6 99]); 

%% Read in DSS Settlement Data Top30 COB Data:
% China = row 4; Vietname = row 30; Philippines = row 20
DSSYears = (1991:2003); %years we have age dist of migrants
%ABSYears = (2004:2019);
% Top30_2017Prev: Pt Est(x = 6); Low CI Est (x=4); Upper CI Est (x=5)
Top30_2017Prev = dlmread('Inputs/Top30_COB_Prev.csv',',',[1 params.DSS_MigPrev 30 params.DSS_MigPrev]); %current prev estimates for COBs (Jenns sys rev, 2017)
Top30_2019Prev = dlmread('Inputs/Top30_COB_Prev.csv',',',[1 params.DSS_MigPrev2019 30 params.DSS_MigPrev2019]); %current prev estimates for COBs (Jenns sys rev, 2017)
Top30_Prev_high = params.Top30_Prev_high; %dlmread('Inputs/Top30_COB_Prev.csv',',',[1 2 30 2]); %higher prev estimates for COBs (Cowie PhD)
regions_top30_DSS = params.regions_top30_DSS; %dlmread('Inputs/Top30_COB_Prev.csv',',',[1 3 30 3]); %current prev estimates for COBs (Jenns sys rev, 2017) 
%DSS_Top30Age = dlmread('Inputs/National_Top30COB_AgeDist_1991_to_2016.csv',',',[1 2 780 8]);
DSS_Top30Age = params.DSS_Top30Age; %dlmread('Inputs/National_Top30COB_AgeDist_1991_to_2016_nz_data.csv',',',[1 2 780 8]); %With NZ UPDATE
%ABS_Top30Age = params.DSS_Top30Age2004to2019; %480 x 18

% Remaining COB's of migrants & corresponding prevalence ests. 1991 to 2016
DSSbyCOB_other_1991to2019 = params.DSSbyCOB_other_1991to2003; %dlmread('Inputs/DSS_other_mig_1991to2016.csv',',',[1 1 240 26]); 
PrevDSSbyCOB_other_1991to2003 = params.PrevDSSbyCOB_other_1991to2003; %dlmread('Inputs/DSS_other_prev_Nat.csv',',',[1 1 240 1])/100; %MacLachlan Ests
PrevDSSbyCOB_other2019 = params.PrevDSSbyCOB_other2019;
regions_other_DSS = params.regions_other_DSS; %dlmread('Inputs/DSS_other_prev_Nat.csv',',',[1 2 240 2]);

%% Years of birth and indices
%yob = 1931:2016; %years of birth for migrants entering Aus between 1991 - 2016
yob1 = 1931:y_end;
n_cob = length(regions_top30_DSS); %number of countries = top 30
%n_other_cob = length(regions_other_DSS); %240
age_index_l = [9:10:59, 60];
age_index_u = (0:10:60);


%% Estimate proportion of NOM entering through 'other' DSS COB
%-------------------------------------------------------------------------
DSS_top30_totals = zeros(30,length(DSSYears));
DSSYearsLabel =  {zeros(1,length(DSSYears))};
for y= 1:length(DSSYears) %loop through years we have DSS settlement data for
  DSSYearsLabel{y}=sprintf('y%d',1990+y); %Labels to be used to store matrices
  %filename = sprintf('Inputs/DSSTop30_%d.CSV',1990+ii) % Create string for each file
  DSS_raw.(DSSYearsLabel{y}) = DSS_Top30Age( (y-1)*n_cob+1:y*n_cob , :); % Save data for each year  
  DSS_top30_totals(:,y) = sum(DSS_raw.(DSSYearsLabel{y}),2); % (i,j) = total for cob i in year j
end
DSS_ALL = sum(DSS_top30_totals) + sum(DSSbyCOB_other_1991to2019); %Total DSS Settlers in each year 1991 to 2016 incl.

DSS_COB_other_Prop = zeros(size(DSSbyCOB_other_1991to2019));
for y= 1:length(DSSYears) 
DSS_COB_other_Prop(:,y) = DSSbyCOB_other_1991to2019(:,y)./DSS_ALL(y);
end
Mean_DSS_COB_other_Prop = mean(DSS_COB_other_Prop,2); %Determine the mean proportion by COB for DSS data over 1991 to 2016

DSSbyCOB_other_1991to2003 = zeros(length(Mean_DSS_COB_other_Prop), length(1991:y_end));
DSSbyCOB_other_1991to2003(:,1:length(DSSYears)) = DSSbyCOB_other_1991to2019;

%-----------------------------------------------
%% All COB after Top 30 excluded:
%-----------------------------------------------

% Determine propoprtions resolved and susceptible (assuming no vaccination)
r = zeros(length(PrevDSSbyCOB_other_1991to2003),2);
s_DSS = zeros(length(PrevDSSbyCOB_other_1991to2003),2);
for k=1:length(PrevDSSbyCOB_other_1991to2003)
%k
r(k,1) = 7*PrevDSSbyCOB_other_1991to2003(k); %7 was the ratio of prop resolved to prop chronic from Cowie2009 serosurveys.     
r(k,2) = 7*PrevDSSbyCOB_other2019(k);

if r(k,1)>0.75
    r(k,1)=0.75;
end

if r(k,2)>0.75
    r(k,2)=0.75;
end
%r(k)
s_DSS(k,1) = 1 - PrevDSSbyCOB_other_1991to2003(k) - r(k,1);
s_DSS(k,2) = 1 - PrevDSSbyCOB_other2019(k) - r(k,2);
end

%countries by years
inf_other_DSS = zeros(length(PrevDSSbyCOB_other_1991to2003) , length(1991:y_end));
sus_other_DSS = zeros(length(PrevDSSbyCOB_other_1991to2003), length(1991:y_end));
res_other_DSS = zeros(length(PrevDSSbyCOB_other_1991to2003), length(1991:y_end));
for i = 1:length(1991:y_end)%y_end) %60
    inf_other_DSS(:,i) = DSSbyCOB_other_1991to2003(:,i).*PrevDSSbyCOB_other_1991to2003;
    sus_other_DSS(:,i) = DSSbyCOB_other_1991to2003(:,i).*(s_DSS(:,1));
    res_other_DSS(:,i) = DSSbyCOB_other_1991to2003(:,i).*(r(:,1));
end

% Total numbers in Resolved & Susceptible states immigrating
R_all_other = zeros(7, length(1991:y_end));
sumRes = sum(res_other_DSS)';  %will end up been a 7 by 26 = summing up those in R from each COB
S_all_other = zeros(7, length(1991:y_end));
sumSus = sum(sus_other_DSS)'; 


% Derive number of chronic infections from each region per year
% to apply phase distribution to later:
s=size(inf_other_DSS);
I_AfricaDSS = zeros(1,s(2));
I_AsiaDSS = zeros(1,s(2));
I_OtherDSS = zeros(1,s(2));

for g=1:s(1)
if regions_other_DSS(g)==1
    I_AfricaDSS(1,:) = I_AfricaDSS(1,:) + inf_other_DSS(g,:);
elseif regions_other_DSS(g)==2
    I_AsiaDSS(1,:) = I_AsiaDSS(1,:) + inf_other_DSS(g,:);
else
    I_OtherDSS(1,:) = I_OtherDSS(1,:) + inf_other_DSS(g,:);
end
end

%test_NOM_age=NOM_age'

I_AfricaOther = zeros(7,length(1991:y_end));
I_AsiaOther = zeros(7,length(1991:y_end));
I_OtherOther = zeros(7,length(1991:y_end));
% Apply NOM age distribution to each year to get total # in each age
% group in each region 1991 to 2050

for ll = 1:length(1991:y_end)
    R_all_other(:,ll) = sumRes(ll)*NOM_age(:,40+ll);
    S_all_other(:,ll) = sumSus(ll)*NOM_age(:,40+ll);
    
    I_AfricaOther(:,ll) = I_AfricaDSS(ll)*NOM_age(:,40+ll);
    I_AsiaOther(:,ll) = I_AsiaDSS(ll)*NOM_age(:,40+ll);
    I_OtherOther(:,ll) = I_OtherDSS(ll)*NOM_age(:,40+ll);    
end
%R_all_other_test = R_all_other(:,1:26)'



%-----------------------------------------------
%% Top 30 COB Estimates for 1991 to 2050:
%-----------------------------------------------
%  Generate prev_yob matrix
%% China = row 4;  Prevalence estimates based on birth year from Cui2017
%prev_top30_yob(a,j) = prev in country a for people born in year j.
prev_top30_yob = zeros(n_cob,length(yob1)); %estimated prevalence in each country over time
% Assign 2017 estimated prevalence to all years for each COB until 2019
% then use 2019 prev estimates.
for i=1:n_cob
    for j=1:length(1931:1990)
    prev_top30_yob(i,j) = Top30_Prev_high(i); %Top30_Prev_high(i);
    prev_top30_yob(17,j) = Top30_2017Prev(17); %NZ
    end
    for j=(length(1931:1990)+1):length(1931:2018) %set all COB prev to be same as 2017 est.
    prev_top30_yob(i,j) = Top30_2017Prev(i);
    end
    for j=(length(1931:2018)+1):length(1931:y_end) %set all COB prev to be same as 2017 est.
    prev_top30_yob(i,j) = Top30_2019Prev(i);
    end
end


%% China prevalence: those born during 1931 to 1978
for i=1:length(1931:1978)
    prev_top30_yob(4,i) = 0.105; %prevalence in birth year in China
end
% China prevalence: those born during 1979 to 1985
for i=length(1931:1978)+1:length(1931:1985)
    prev_top30_yob(4,i) = 0.095; %prevalence in birth year in China
end
% Smooth out prevalence estimates a bit from 0.105 in 1978 to 0.0950 in 1985
prev_top30_yob(4,47:50)=(0.1050:-(0.1050-0.0950)/3:0.0950);
prev_top30_yob(4,50:55) = (0.0950:-(0.0950-0.09)/5:0.09);
% China prevalence: those born during 1986 to 1990
prev_top30_yob(4,length(1931:1985)+1:length(1931:1990)) = [0.08, 0.075, 0.07, 0.06, 0.05] ; %prevalence in birth year in China
% China prevalence: those born during 1991 to 1996
prev_top30_yob(4,length(1931:1990)+1:length(1931:1996)) = (5:-(5-3)/5:3)./100 ; %prevalence in birth year in China
%1991:1996 = 5:-(5-3)/5:3
% China prevalence: those born during 1997 to 2013
prev_top30_yob(4,length(1931:1996)+1:length(1931:2013)) = (3:-(3-0.3)/16:0.3)./100 ; %prevalence in birth year in China
% China prevalence: those born during 2014 to 2016
prev_top30_yob(4,length(1931:2013)+1:length(1931:2016)) = [0.3, 0.2, 0.2]./100 ; %prevalence in birth year in China
%plot(prev_top30_yob(4,:),'ro')
%%China prevalence: those born during 2017 to 2050
prev_top30_yob(4,length(1931:2016)+1:length(1931:2050)) = 0.2/100;

%% Philippines = row 20; Little data indicates vaccine hasn't had much effect
% Taking higher est prev = 9.8% pre-2005 and 3.02% from 2005 onwards.
for i=1:length(1931:2004)
prev_top30_yob(20,i) = 0.098;
end

%% Vietnam = row 30; 
% Use age- and time- dependent prevalence for China, Vietnam & Philippines
% Vietnam: vacc in 1997

for i=1:length(1931:1951)
   prev_top30_yob(30,i) = 0.12; %prevalence in birth year in VN
end

for i=length(1931:1952):length(1931:1970)
    prev_top30_yob(30,i) = prev_top30_yob(30,i-1)-((0.12-0.1)/(1970-1952+1));
end

for i=length(1931:1971):length(1931:1990)
    prev_top30_yob(30,i) = prev_top30_yob(30,i-1)-((0.1-0.085)/(1990-1971+1));
end

prev_top30_yob(30,length(1931:1991):length(1931:y_end)) = [8.48  8.4  8.27  8.09  7.89  7.54  6.97  6.28  5.55  4.90  4.23  3.50  2.79]/100;


%% Taiwan prevalence: row 26. %New for 2019 Model/2018 report. 
% From Su et al. 2018 ('Impact of Universal Infant Hepatitis B Immunisation
% on Reducing the Hepatitis B Carrier Rate in Pregnant Women')
for i=1:length(1931:1962)
prev_top30_yob(26,i) = 0.148;
end
prev_top30_yob(26,length(1931:1963):length(1931:1993))=[14.8	15	14.8	15.3	15	15.3	15	14.7	14.5	14.3	14	13.3	13.3	12.8	12.4	11.8	11.6	11.2	10.1	9.2	8.1	7	5.7	4.5	2.8	2.8	2.3	2.3	2.3	2	1.8]./100;
for i=length(1931:1994):length(1931:1999)
prev_top30_yob(26,i) = 0.018;
end
for i=length(1931:2000):length(1931:2050)
prev_top30_yob(26,i) = 0.01; %was 1%
end

%NOM_tot(41:41+25) - sum(DSSbyCOB_other_1991to2018)
k = zeros(n_cob,7);
rr = zeros(n_cob,7);
%ss = zeros(n_cob,7);
%COB_dist_top30 = zeros(n_cob,length(DSSYears));
%COB_prop_ofNOM_top30 = zeros(n_cob,length(DSSYears));


DSS_COB_top30_Prop = zeros(size(DSS_top30_totals));
for y= 1:length(DSSYears) 
DSS_COB_top30_Prop(:,y) = DSS_top30_totals(:,y)./DSS_ALL(y);
end

age_lengths = length(age_index_l);
r1 = zeros(age_lengths,1);
p = zeros(age_lengths,1);
DSSYearsFuture = 1991:2003;

for y= 1:length(DSSYearsFuture) %loop through years we have DSS settlement data for
  DSSYearsLabel{y}=sprintf('y%d',1990+y); %Labels to be used to store matrices
  %filename = sprintf('Inputs/DSSTop30_%d.CSV',1990+ii) % Create string for each file
  %DSS_raw2.(DSSYearsLabel{y}) = DSS_Top30Age_w_proj( (y-1)*n_cob+1:y*n_cob , :); % Save data for each year   
  DSS_raw2.(DSSYearsLabel{y}) = DSS_Top30Age( (y-1)*n_cob+1:y*n_cob , :); % Save data for each year   
  
    for z=1:n_cob
        for j=1:age_lengths %length(age_index_l) %loop through each age group
           
            a=length(1931:years(y)-age_index_l(j));
            b=length(1931:years(y)-age_index_u(j));

            %p(j) = mean( prev_top30_yob(z, a:b ) ); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            p(j) = sum( prev_top30_yob(z, a:b ) )/length(a:b); 
            
            %k(z,j) = number indivs from COB z entering in age group j that
            %         have CHB
            k(z,j) = DSS_raw2.(DSSYearsLabel{y})(z,j)*p(j); %mult # migrants in age group by corresp prev when born
            
                r1(j) = 7*p(j);
                if r1(j)>0.75
                r1(j)=0.75;
                end
                
            rr(z,j) = DSS_raw2.(DSSYearsLabel{y})(z,j)*r1(j); %mult # migrants in age group by corresp est in resolved group when born
            %ss(z,j) = DSS_raw2.(DSSYearsLabel{y})(z,j)*(1-r1(j)-p(j)); %number susceptibles in each age group.
            %v(j); %mult # migrants in age group by corresp est in vacc group when born
        
        end
    end

%% Get Total numbers of immigrants by COB by Time    
%COB_dist_top30(:,y) =   sum(DSS_raw.(DSSYearsLabel{y})');

%COB_prop_ofNOM_top30(:,y) =  COB_dist_top30(:,y)./NOM_tot(40+y);

% Create structure containing filled in k(z,j) matrix for each year
% DSS_CHB.y1991 = k(z,j) for 1991 to DSS_CHB.y2016 = k(z,j) for 2016
DSS_CHB.(DSSYearsLabel{y}) = k;
DSS_Res.(DSSYearsLabel{y}) = rr; % number in resolved class %check serosurvey data from Cowie 2009
%DSS_V.(DSSYearsLabel{y}) = v; % number in immune/vacc class

ss = DSS_raw2.(DSSYearsLabel{y}) - k - rr; 
DSS_S.(DSSYearsLabel{y}) = ss; %remaining pop = Susceptible

end

% debgugging: fine up to here

% Total number chronic infections by age coming from African region:
I_AfricaDSS_age = zeros(7,length(1991:y_end)); % age by 26 matrix 
I_AsiaDSS_age = zeros(7,length(1991:y_end));  % age by 26 matrix 
I_OtherDSS_age = zeros(7,length(1991:y_end)); % age by 26 matrix 
R_all = zeros(7, length(DSSYearsFuture));
S_all = zeros(7, length(DSSYearsFuture));
%
for y= 1:length(DSSYearsFuture) %loop through years we have DSS settlement data for
    DSSYearsLabel{y}=sprintf('y%d',1990+y); %Labels to be used to store matrices
    
    CHB = DSS_CHB.(DSSYearsLabel{y}); % Get matrix of cob by age with # inf
    
    for g=1:n_cob
       
        if regions_top30_DSS(g)==1
            I_AfricaDSS_age(:,y) = I_AfricaDSS_age(:,y) + CHB(g,:)';
        elseif regions_top30_DSS(g)==2 
            I_AsiaDSS_age(:,y) = I_AsiaDSS_age(:,y) + CHB(g,:)';
        else
            I_OtherDSS_age(:,y) = I_OtherDSS_age(:,y) + CHB(g,:)';
        end
    end
    
    R_all(:,y) = sum(DSS_Res.(DSSYearsLabel{y}))'; %7 by 26
    S_all(:,y) = sum(DSS_S.(DSSYearsLabel{y}))'; %7 by 26
    
end

%R_test = [R_all(:,1:26)', R_all_other(:,1:26)'];

% Add the total in R and S from 'other COB' from DSS:
R_all = R_all + R_all_other;
S_all = S_all + S_all_other;

I_AfricaDSS_age = I_AfricaDSS_age+ I_AfricaOther;
I_AsiaDSS_age = I_AsiaDSS_age + I_AsiaOther;
I_OtherDSS_age = I_OtherDSS_age+ I_OtherOther;


%% Disease phase by regions:
% Distribute the total number with chronic infection across the 4 phases by
% age


%% June 12th update uses BC_estimates2 (prev = BC_estimates.csv)
% Ages 0 to 9; 10 to 19; 20 to 29; 30 to 39; 40 to 49; 50 to 59; 60+ 
phase_dist_Africa = params.phase_dist_Africa; % dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[2 1 8 4]);
phase_dist_Asia = params.phase_dist_Asia; %dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[12 1 18 4]);
phase_dist_combined = params.phase_dist_combined; %dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[45 1 51 4]);

% Assign each COB a region.

I_all=zeros(7,length(1991:y_end));
I1_all=zeros(7,length(1991:y_end));
I2_all=zeros(7,length(1991:y_end));
I3_all=zeros(7,length(1991:y_end));

for i=1:length(1991:y_end)
I_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,1) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,1) + I_OtherDSS_age(:,i).*phase_dist_combined(:,1));
I1_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,2) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,2) + I_OtherDSS_age(:,i).*phase_dist_combined(:,2)); 
I2_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,3) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,3) + I_OtherDSS_age(:,i).*phase_dist_combined(:,3));
I3_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,4) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,4) + I_OtherDSS_age(:,i).*phase_dist_combined(:,4));
end


%test2 = [I_all'; I1_all'; I2_all'; I3_all']

%% Convert age groups from: 0-9; 10-19; 20-29; 30-39; 40-49; 50-59; 60+
% To ages: 0-4;5-9;  10-14;15-19;  20-24;25-29;  30-34;35-39; 40-44;45-49;
% 50-54;55-59;   60-64;65-69;70-74;75-79;80-84;85+
mig_I = ([I_all(1,:)/2;I_all(1,:)/2;  I_all(2,:)/2; I_all(2,:)/2;   I_all(3,:)/2; I_all(3,:)/2;   I_all(4,:)/2; I_all(4,:)/2;   I_all(5,:)/2; I_all(5,:)/2;    I_all(6,:)/2; I_all(6,:)/2;            I_all(7,:)/6; I_all(7,:)/6;I_all(7,:)/6; I_all(7,:)/6;I_all(7,:)/6; I_all(7,:)/6;]);              
mig_I1 = ([I1_all(1,:)/2;I1_all(1,:)/2;  I1_all(2,:)/2; I1_all(2,:)/2;   I1_all(3,:)/2; I1_all(3,:)/2;   I1_all(4,:)/2; I1_all(4,:)/2;   I1_all(5,:)/2; I1_all(5,:)/2;    I1_all(6,:)/2; I1_all(6,:)/2;            I1_all(7,:)/6; I1_all(7,:)/6;I1_all(7,:)/6; I1_all(7,:)/6;I1_all(7,:)/6; I1_all(7,:)/6;]);
mig_I2 = ([I2_all(1,:)/2;I2_all(1,:)/2;  I2_all(2,:)/2; I2_all(2,:)/2;   I2_all(3,:)/2; I2_all(3,:)/2;   I2_all(4,:)/2; I2_all(4,:)/2;   I2_all(5,:)/2; I2_all(5,:)/2;    I2_all(6,:)/2; I2_all(6,:)/2;            I2_all(7,:)/6; I2_all(7,:)/6;I2_all(7,:)/6; I2_all(7,:)/6;I2_all(7,:)/6; I2_all(7,:)/6;]);
mig_I3 = ([I3_all(1,:)/2;I3_all(1,:)/2;  I3_all(2,:)/2; I3_all(2,:)/2;   I3_all(3,:)/2; I3_all(3,:)/2;   I3_all(4,:)/2; I3_all(4,:)/2;   I3_all(5,:)/2; I3_all(5,:)/2;    I3_all(6,:)/2; I3_all(6,:)/2;            I3_all(7,:)/6; I3_all(7,:)/6;I3_all(7,:)/6; I3_all(7,:)/6;I3_all(7,:)/6; I3_all(7,:)/6;]);
mig_R = ([R_all(1,:)/2;R_all(1,:)/2;  R_all(2,:)/2; R_all(2,:)/2;   R_all(3,:)/2; R_all(3,:)/2;   R_all(4,:)/2; R_all(4,:)/2;   R_all(5,:)/2; R_all(5,:)/2;    R_all(6,:)/2; R_all(6,:)/2;            R_all(7,:)/6; R_all(7,:)/6;R_all(7,:)/6; R_all(7,:)/6;R_all(7,:)/6; R_all(7,:)/6;]);
mig_S = ([S_all(1,:)/2;S_all(1,:)/2;  S_all(2,:)/2; S_all(2,:)/2;   S_all(3,:)/2; S_all(3,:)/2;   S_all(4,:)/2; S_all(4,:)/2;   S_all(5,:)/2; S_all(5,:)/2;    S_all(6,:)/2; S_all(6,:)/2;            S_all(7,:)/6; S_all(7,:)/6;S_all(7,:)/6; S_all(7,:)/6;S_all(7,:)/6; S_all(7,:)/6;]);

%mig_check = [mig_I(:,1:26);mig_I1(:,1:26);mig_I2(:,1:26);mig_I3(:,1:26);mig_R(:,1:26);mig_S(:,1:26)];

%mig_chronic = (mig_I + mig_I1 + mig_I2 + mig_I3)'; 

%mig_all_age = (mig_I + mig_I1 + mig_I2 + mig_I3+mig_R + mig_S)';
mig_all = sum(mig_I + mig_I1 + mig_I2 + mig_I3+mig_R + mig_S)';


% If using ABS NOM by COB for 2005 to 2018 - shouldn't scale further so,
% for 2005 onwards scale = 1. Scale from 1991 to 2004.
ratio1 = NOM_tot(41:53)./mig_all;
%ratio=ones(length(NOM_tot(41:100)),1);
%ratio= [ratio1(1:14);  ones(60-14,1)];


for i = 1:length(ratio1)
mig_I(:,i) = round( mig_I(:,i)*ratio1(i));
mig_I1(:,i) =round( mig_I1(:,i)*ratio1(i));
mig_I2(:,i) =round( mig_I2(:,i)*ratio1(i));
mig_I3(:,i) =round(mig_I3(:,i)*ratio1(i));
mig_R(:,i)= round(mig_R(:,i)*ratio1(i));
mig_S(:,i)= round(mig_S(:,i)*ratio1(i));
end


end



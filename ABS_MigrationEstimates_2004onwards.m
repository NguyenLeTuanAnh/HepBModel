%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATED JUNE 2022 - ABS NOM by COB and AGE for 2004 - 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mig_Sa, mig_Ia, mig_I1a, mig_I2a, mig_I3a, mig_Ra] = ABS_MigrationEstimates_2004onwards(params)

y_end = 2050;
years=2004:y_end;

%% Read in DSS Settlement Data Top30 COB Data:
% China = row 4; Vietname = row 30; Philippines = row 20
ABSYears = (2004:2021); %2004 is year we started to have age dist of migrants
Years_to_Average = 12:length(ABSYears)-2; %2015 - 2019
Years_to_Average_ages = 12:length(ABSYears)-3; %2015 - 2018 %excludes extremes of 2019 in China (-ve under 40s, very high +ves over 40s)

%ABSYears = (2004:2019);
% Top30_2017Prev: Pt Est(x = 6); Low CI Est (x=4); Upper CI Est (x=5)
Top30_2017Prev = dlmread('Inputs/Top30_COB_Prev_from2004.csv',',',[1 params.DSS_MigPrev 30 params.DSS_MigPrev]); %current prev estimates for COBs (Jenns sys rev, 2017)
Top30_2019Prev = dlmread('Inputs/Top30_COB_Prev_from2004.csv',',',[1 params.DSS_MigPrev2019 30 params.DSS_MigPrev2019]); %current prev estimates for COBs (Jenns sys rev, 2017)
Top30_Prev_high = params.Top30_Prev_high_from2004; %dlmread('Inputs/Top30_COB_Prev.csv',',',[1 2 30 2]); %higher prev estimates for COBs (Cowie PhD)
regions_top30_DSS = params.regions_top30_DSS_from2004; %dlmread('Inputs/Top30_COB_Prev.csv',',',[1 3 30 3]); %current prev estimates for COBs (Jenns sys rev, 2017) 
%DSS_Top30Age = dlmread('Inputs/National_Top30COB_AgeDist_1991_to_2016.csv',',',[1 2 780 8]);
%DSS_Top30Age = params.DSS_Top30Age(1:(length(ABSYears)*30),:); %dlmread('Inputs/National_Top30COB_AgeDist_1991_to_2016_nz_data.csv',',',[1 2 780 8]); %With NZ UPDATE
ABS_Top30Age = params.ABS_Top30Age2004to2021; %510 x 18

% Remaining COB's of migrants & corresponding prevalence ests. 2004 to 2020
ABSbyCOB_other_from2004 = params.ABSbyCOB_other_from2004;
PrevABSbyCOB_other_from2004 = params.PrevABSbyCOB_other_from2004; %dlmread('Inputs/DSS_other_prev_Nat.csv',',',[1 1 240 1])/100; %MacLachlan Ests
PrevABSbyCOB_other2019 = params.PrevABSbyCOB_other2019;
regions_other_ABS = params.regions_other_ABS; %dlmread('Inputs/DSS_other_prev_Nat.csv',',',[1 2 240 2]);

%% Years of birth and indices
%yob1 = years of birth for migrants entering Aus between 2004 -
%2020
yob1 = 1919:y_end;
n_cob = length(regions_top30_DSS); %number of countries = top 30
%n_other_cob = length(regions_other_DSS); %240
age_index_l = [4:5:85, 85]; %[9:10:59, 60];
%age_index_u = (0:5:85); %(0:10:60);


%% Estimate proportion of NOM entering through 'other' DSS COB
%-------------------------------------------------------------------------
ABS_top30_totals = zeros(30,length(ABSYears));
ABSYearsLabel =  {zeros(1,length(ABSYears))};
for y= 1:length(ABSYears) %loop through years we have DSS settlement data for
  ABSYearsLabel{y}=sprintf('y%d',2003 +y); %Labels to be used to store matrices
  %filename = sprintf('Inputs/DSSTop30_%d.CSV',1990+ii) % Create string for each file
  ABS_raw.(ABSYearsLabel{y}) = ABS_Top30Age( (y-1)*n_cob+1:y*n_cob , :); % Save data for each year  
  ABS_top30_totals(:,y) = sum(ABS_raw.(ABSYearsLabel{y}),2); % (i,j) = total for cob i in year j
end
ABS_ALL = sum(ABS_top30_totals) + sum(ABSbyCOB_other_from2004); %Total ABS NOM

ABS_COB_other_Prop = zeros(size(ABSbyCOB_other_from2004));
for y= 1:length(ABSYears) 
ABS_COB_other_Prop(:,y) = ABSbyCOB_other_from2004(:,y)./ABS_ALL(y);
end
Mean_ABS_COB_other_Prop = mean(ABS_COB_other_Prop(:,Years_to_Average),2); %Determine the mean proportion by COB for ABS data over 5 yrs to 2019

ABSbyCOB_other_2004to2050 = zeros(length(Mean_ABS_COB_other_Prop), length(2004:y_end));
ABSbyCOB_other_2004to2050(:,1:length(ABSYears)) = ABSbyCOB_other_from2004;

for yy = (length(ABSYears)+1):length(years) 
    ABSbyCOB_other_2004to2050(:,yy) = (Mean_ABS_COB_other_Prop)*params.NOM_tot(53+yy); %project numbers by 'other' cob forward to 2050
end    
%-------------------------------------------------------------------------


%-----------------------------------------------
%% All COB after Top 30 excluded:
%-----------------------------------------------
% Determine propoprtions resolved and susceptible (assuming no vaccination)
r = zeros(length(PrevABSbyCOB_other_from2004),2);
s_DSS = zeros(length(PrevABSbyCOB_other_from2004),2);
for k=1:length(PrevABSbyCOB_other_from2004)
%k
r(k,1) = 7*PrevABSbyCOB_other_from2004(k); %7 was the ratio of prop resolved to prop chronic from Cowie2009 serosurveys.     
r(k,2) = 7*PrevABSbyCOB_other2019(k);

if r(k,1)>0.75
    r(k,1)=0.75;
end

if r(k,2)>0.75
    r(k,2)=0.75;
end
%r(k)
s_DSS(k,1) = 1 - PrevABSbyCOB_other_from2004(k) - r(k,1);
s_DSS(k,2) = 1 - PrevABSbyCOB_other2019(k) - r(k,2);
end

%countries by years
inf_other_ABS = zeros(length(PrevABSbyCOB_other_from2004) , length(2004:y_end));
sus_other_ABS = zeros(length(PrevABSbyCOB_other_from2004), length(2004:y_end));
res_other_ABS = zeros(length(PrevABSbyCOB_other_from2004), length(2004:y_end));
for i = 1:length(2004:2018)%y_end) %60
    inf_other_ABS(:,i) = ABSbyCOB_other_2004to2050(:,i).*PrevABSbyCOB_other_from2004;
    sus_other_ABS(:,i) = ABSbyCOB_other_2004to2050(:,i).*(s_DSS(:,1));
    res_other_ABS(:,i) = ABSbyCOB_other_2004to2050(:,i).*(r(:,1));
end

% Include updated Prevalence estimates for 2019 onwards:
for i = length(2004:2019):length(2004:y_end) %60
    inf_other_ABS(:,i) = ABSbyCOB_other_2004to2050(:,i).*PrevABSbyCOB_other2019;
    sus_other_ABS(:,i) = ABSbyCOB_other_2004to2050(:,i).*(s_DSS(:,2));
    res_other_ABS(:,i) = ABSbyCOB_other_2004to2050(:,i).*(r(:,2));
end

%s_DSS+r+PrevDSSbyCOB_other_1991to2018 %check = 1 for all COB. 

% Total numbers in Resolved & Susceptible states immigrating
R_all_other = zeros(18, length(2004:y_end));
sumRes = sum(res_other_ABS)';  
S_all_other = zeros(18, length(2004:y_end));
sumSus = sum(sus_other_ABS)'; 


% Derive number of chronic infections from each region per year
% to apply phase distribution to later:
s=size(inf_other_ABS);
I_AfricaABS = zeros(1,s(2));
I_AsiaABS = zeros(1,s(2));
I_OtherABS = zeros(1,s(2));

for g=1:s(1)
if regions_other_ABS(g)==1
    I_AfricaABS(1,:) = I_AfricaABS(1,:) + inf_other_ABS(g,:);
elseif regions_other_ABS(g)==2
    I_AsiaABS(1,:) = I_AsiaABS(1,:) + inf_other_ABS(g,:);
else
    I_OtherABS(1,:) = I_OtherABS(1,:) + inf_other_ABS(g,:);
end
end

%test_NOM_age=NOM_age'

I_AfricaOther = zeros(18,length(2004:y_end));
I_AsiaOther = zeros(18,length(2004:y_end));
I_OtherOther = zeros(18,length(2004:y_end));

%% Apply derived average age distribution of all 'other COB' combined
%% to each other COB: (average across 2004 - 2019 is applied forward)
for ll = 1:length(2004:y_end)
    R_all_other(:,ll) = sumRes(ll)*params.ABS_Other_agedist(:,ll);% NOM_age(:,40+ll);
    S_all_other(:,ll) = sumSus(ll)*params.ABS_Other_agedist(:,ll);%NOM_age(:,40+ll);
    
    I_AfricaOther(:,ll) = I_AfricaABS(ll)*params.ABS_Other_agedist(:,ll);%NOM_age(:,40+ll);
    I_AsiaOther(:,ll) = I_AsiaABS(ll)*params.ABS_Other_agedist(:,ll);%NOM_age(:,40+ll);
    I_OtherOther(:,ll) = I_OtherABS(ll)*params.ABS_Other_agedist(:,ll);%NOM_age(:,40+ll);    
end
%R_all_other_test = R_all_other(:,1:26)'


%-----------------------------------------------
%% Top 30 COB Estimates for 2004 to 2050:
%-----------------------------------------------
%  Generate prev_yob matrix

%% China = row 4;  Prevalence estimates based on birth year from Cui2017
%prev_top30_yob(a,j) = prev in country a for people born in year j.
prev_top30_yob = zeros(n_cob,length(yob1)); %estimated prevalence in each country over time
% Assign 2017 estimated prevalence to all years for each COB until 2019
% then use 2019 prev estimates.
for i=1:n_cob
    for j=1:length(1919:1990)
    prev_top30_yob(i,j) = Top30_Prev_high(i); %Top30_Prev_high(i);
    prev_top30_yob(16,j) = Top30_2017Prev(16); %NZ
    end
    for j=(length(1919:1990)+1):length(1919:2018) %set all COB prev to be same as 2017 est.
    prev_top30_yob(i,j) = Top30_2017Prev(i);
    end
    for j=(length(1919:2018)+1):length(1919:y_end) %set all COB prev to be same as 2017 est.
    prev_top30_yob(i,j) = Top30_2019Prev(i);
    end
end


%% China prevalence: those born during 1931 to 1978
for i=1:length(1919:1978)
    prev_top30_yob(4,i) = 0.105; %prevalence in birth year in China
end
% China prevalence: those born during 1979 to 1985
for i=length(1919:1978)+1:length(1919:1985)
    prev_top30_yob(4,i) = 0.095; %prevalence in birth year in China
end
% Smooth out prevalence estimates a bit from 0.105 in 1978 to 0.0950 in 1985
prev_top30_yob(4,59:62)=(0.1050:-(0.1050-0.0950)/3:0.0950);
prev_top30_yob(4,62:67) = (0.0950:-(0.0950-0.09)/5:0.09);
% China prevalence: those born during 1986 to 1990
prev_top30_yob(4,length(1919:1985)+1:length(1919:1990)) = [0.08, 0.075, 0.07, 0.06, 0.05] ; %prevalence in birth year in China
% China prevalence: those born during 1991 to 1996
prev_top30_yob(4,length(1919:1990)+1:length(1919:1996)) = (5:-(5-3)/5:3)./100 ; %prevalence in birth year in China
%1991:1996 = 5:-(5-3)/5:3
% China prevalence: those born during 1997 to 2013
prev_top30_yob(4,length(1919:1996)+1:length(1919:2013)) = (3:-(3-0.3)/16:0.3)./100 ; %prevalence in birth year in China
% China prevalence: those born during 2014 to 2016
prev_top30_yob(4,length(1919:2013)+1:length(1919:2016)) = [0.3, 0.2, 0.2]./100 ; %prevalence in birth year in China
%plot(prev_top30_yob(4,:),'ro')
%%China prevalence: those born during 2017 to 2050
prev_top30_yob(4,length(1919:2016)+1:length(1919:2050)) = 0.2/100;

%% Philippines = row 19; Little data indicates vaccine hasn't had much effect
% Taking higher est prev = 9.8% pre-2005 and 3.02% from 2005 onwards.
for i=1:length(1919:2004)
prev_top30_yob(19,i) = 0.098;
end

%% Vietnam = row 30; 
% Use age- and time- dependent prevalence for China, Vietnam & Philippines
% Vietnam: vacc in 1997

for i=1:length(1919:1951)
   prev_top30_yob(30,i) = 0.12; %prevalence in birth year in VN
end

for i=length(1919:1952):length(1919:1970)
    prev_top30_yob(30,i) = prev_top30_yob(30,i-1)-((0.12-0.1)/(1970-1952+1));
end

for i=length(1919:1971):length(1919:1990)
    prev_top30_yob(30,i) = prev_top30_yob(30,i-1)-((0.1-0.085)/(1990-1971+1));
end

prev_top30_yob(30,length(1919:1991):length(1919:2009)) = [8.48  8.4  8.27  8.09  7.89  7.54  6.97  6.28  5.55  4.90  4.23  3.50  2.79   2.23   1.93   1.82   1.75   1.69   1.64]/100;
prev_top30_yob(30,length(1919:2010):length(1919:y_end)) = 1.57/100;

%% Taiwan prevalence: row 25. %New for 2021 Model/2020 report. 
% From Su et al. 2018 ('Impact of Universal Infant Hepatitis B Immunisation
% on Reducing the Hepatitis B Carrier Rate in Pregnant Women')
for i=1:length(1919:1962)
prev_top30_yob(25,i) = 0.148;
end
prev_top30_yob(25,length(1919:1963):length(1919:1993))=[14.8	15	14.8	15.3	15	15.3	15	14.7	14.5	14.3	14	13.3	13.3	12.8	12.4	11.8	11.6	11.2	10.1	9.2	8.1	7	5.7	4.5	2.8	2.8	2.3	2.3	2.3	2	1.8]./100;
for i=length(1919:1994):length(1919:1999)
prev_top30_yob(25,i) = 0.018;
end
for i=length(1919:2000):length(1919:2050)
prev_top30_yob(25,i) = 0.01; %was 1%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Matrix with Top30 COB by age that includes 2022 - 2050:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ABS_COB_top30_Prop = zeros(size(ABS_top30_totals));
for y= 1:length(ABSYears) 
ABS_COB_top30_Prop(:,y) = ABS_top30_totals(:,y)./ABS_ALL(y);
end

%Determine the mean proportion of ABS settlers due to the Top30 COB for the
%last 10years (looked at the mean prop over 5, 10, 15 years. 
%yrs_to_av = [8     9    10    11   12   13    14    15    16    17]; %last 10 years excl 2015
%Mean_ABS_COB_top30_Prop = mean(ABS_COB_top30_Prop(:,length(ABSYears)-4:length(ABSYears)),2); % Average over last 15 years.
%Mean_ABS_COB_top30_Prop = mean(ABS_COB_top30_Prop(:,yrs_to_av),2); % Average over last 10 years.
Mean_ABS_COB_top30_Prop = mean(ABS_COB_top30_Prop(:,Years_to_Average),2); % av over 5 years to 2019


% Estimate the total number migrating by COB for the top30 2022 to 2050
ABSbyCOB_Top30_2022to2050 = zeros(n_cob,length(years)); 
for yy = length(ABSYears)+1:length(2004:2050)
ABSbyCOB_Top30_2022to2050(:,yy) = (Mean_ABS_COB_top30_Prop')*params.NOM_tot(53+yy);
end    


% Determine the average age distribution for each COB to use for
% projections
COB_labels={'AFGHANISTAN','BANGLADESH','CAMBODIA','CHINA','ETHIOPIA','GERMANY','GREECE','HKSAR','INDIA',...
    'INDONESIA','ITALY','KENYA','KOREA_SOUTH','MALAYSIA','MYANMAR','NEW_ZEALAND','NIGERIA','PAKISTAN','PHILIPPINES','POLAND', ...
    'SAMOA','SINGAPORE','SOMALIA','SUDAN','TAIWAN','THAILAND','TONGA','TURKEY','UK','VIETNAM'};
for yy= 1:length(ABSYears) %loop through years we have DSS settlement data for
  ABSYearsLabel{yy}=sprintf('y%d',2003+yy); %Labels to be used to store matrices
  %filename = sprintf('Inputs/DSSTop30_%d.CSV',1990+ii) % Create string for each file
  %DSS_raw.(DSSYearsLabel{yy}) = DSS_Top30Age( (yy-1)*n_cob+1:yy*n_cob , :); % Save data for each year   
 
    for z=1:n_cob

    COB_agedist.(ABSYearsLabel{yy})(z,:) =  ABS_raw.(ABSYearsLabel{yy})(z,:)./ABS_top30_totals(z,yy); %cob by age
    
    if ABS_top30_totals(z,yy)<0
        COB_agedist.(ABSYearsLabel{yy})(z,:) = (COB_agedist.(ABSYearsLabel{yy})(z,:))*(-1);
    end
    
    end
    
end

for i=1:length(ABSYears)
    for n=1:n_cob
    IndivCOB_agedist.(COB_labels{n})(:,i) =  COB_agedist.(ABSYearsLabel{i})(n,:)'; %age x year
    %IndivCOB_agedist_numbers.(COB_labels{n})(:,i) =  DSS_raw.(DSSYearsLabel{i})(n,:)';
    end
     
end

%ages = [1:7];
COB_mean_age = zeros(n_cob, 18);
%for t=1:length(DSSYears)
%yrs_to_av = [7     8     9    10    11   13    14    15    16]; %length(2004:2019) - 9:length(2004:2019) excludes 2015
for n=1:n_cob
  %for a=1:7
%COB_mean_age(n,:) = mean( IndivCOB_agedist.(COB_labels{n})(:,length(ABSYears)-9:length(ABSYears)),2); %9 -> 14
COB_mean_age(n,:) = mean( IndivCOB_agedist.(COB_labels{n})(:,Years_to_Average_ages),2);  %av over 4 years to 2018

if n==16
    COB_mean_age(n,:) = mean( IndivCOB_agedist.(COB_labels{n})(:,Years_to_Average_ages+1),2);  %av over 4 years to 2019 for NZ - extreme data in 2012 excluded.
end

%end 
%end
end

% COB_mean_age = params.COB_mean_age; %Average age dist of each Top30 COB during 2014 - 2018.

ABS_Top30Age_w_proj = zeros(length(years)*30,18);
ABS_Top30Age_w_proj(1:length(ABSYears)*30,:) = ABS_Top30Age;
ABS_top30_f = zeros(n_cob ,18);
for kj=(length(ABSYears)+1):length(years)
ABSYearsLabel{kj}=sprintf('k%d',2003+kj);
    for a=1:18
    ABS_top30_f(:,a) = COB_mean_age(:,a).*ABSbyCOB_Top30_2022to2050(:,kj);
    end
ABS_Top30Age_w_proj( (kj-1)*30 +1 : kj*30,: ) = ABS_top30_f; %next
end

age_lengths = length(age_index_l);
r1 = zeros(age_lengths,1);
p = zeros(age_lengths,1);
%DSSYearsFuture = ABSYears; %2004 - 2020;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for y= 1:length(years) %loop through years we have ABS data for
  ABSYearsLabel{y}=sprintf('y%d',2003+y); %Labels to be used to store matrices
  %filename = sprintf('Inputs/DSSTop30_%d.CSV',1990+ii) % Create string for each file 
  %ABS_raw2.(ABSYearsLabel{y}) = ABS_Top30Age( (y-1)*n_cob+1:y*n_cob , :); % Save data for each year   
  ABS_raw2.(ABSYearsLabel{y}) = ABS_Top30Age_w_proj( (y-1)*n_cob+1:y*n_cob , :); % Save data for each year 
  
    for z=1:n_cob
        for j=1:age_lengths %length(age_index_l) %loop through each age group
           
%             a=length(1919:years(y)-age_index_l(j));
%             b=length(1919:years(y)-age_index_u(j));

            %p(j) = sum( prev_top30_yob(z, a:b ) )/length(a:b); %derive average prevalence in country z, in birth decade to corresp age group. p(j) = prev for age group j in decade they were born.
            %p(j) = sum( prev_top30_yob(z, params.yj_a(y,j):params.yj_b(y,j) ) )/length(params.yj_a(y,j):params.yj_b(y,j)); %to make more efficient
            
            p(z,j) = sum( prev_top30_yob(z, params.yj_a(y,j):params.yj_b(y,j) ) )/length(params.yj_a(y,j):params.yj_b(y,j)); %to make more efficient. 
            % p = params.p_matrix.(ABSYearsLabel{y});
             
            %k(z,j) = number indivs from COB z entering in age group j that
            %         have CHB
            %kk(z,j) = ABS_raw2.(ABSYearsLabel{y})(z,j)*p(j); %mult # migrants in age group by corresp prev when born
            
                r1(z,j) = 7*p(z,j);
                if r1(z,j)>0.75
                r1(z,j)=0.75;
                end
                
            %rr(z,j) = ABS_raw2.(ABSYearsLabel{y})(z,j)*r1(j); %mult # migrants in age group by corresp est in resolved group when born  
        end
        
    end
    %% Improving simulation efficiency:                    
     kk = ABS_raw2.(ABSYearsLabel{y}).*(p); %p(j); %mult # migrants in age group by corresp prev when born     
     rr = ABS_raw2.(ABSYearsLabel{y}).*(r1); %mult # migrants in age group by corresp est in resolved group when born    
     %kk = ABS_raw2.(ABSYearsLabel{y}).*(ones(30,1)*(p')); %p(j); %mult # migrants in age group by corresp prev when born     
     %rr = ABS_raw2.(ABSYearsLabel{y}).*(ones(30,1)*(r1')); %mult # migrants in age group by corresp est in resolved group when born   

     %p_matrix.(ABSYearsLabel{y}) = p;
     %save('p_mat.mat','-struct','p_matrix')
     
%% Get Total numbers of immigrants by COB by Time    
%COB_dist_top30(:,y) =   sum(DSS_raw.(DSSYearsLabel{y})');

%COB_prop_ofNOM_top30(:,y) =  COB_dist_top30(:,y)./NOM_tot(40+y);

% Create structure containing filled in k(z,j) matrix for each year
% DSS_CHB.y1991 = k(z,j) for 1991 to DSS_CHB.y2016 = k(z,j) for 2016
ABS_CHB.(ABSYearsLabel{y}) = kk;
ABS_Res.(ABSYearsLabel{y}) = rr; 
%DSS_V.(DSSYearsLabel{y}) = v; % number in immune/vacc class

ss = ABS_raw2.(ABSYearsLabel{y}) - kk - rr; 
ABS_S.(ABSYearsLabel{y}) = ss; %remaining pop = Susceptible

end


% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Total number chronic infections by age coming from African region:
I_AfricaDSS_age = zeros(18,length(2004:y_end)); % age by 26 matrix 
I_AsiaDSS_age = zeros(18,length(2004:y_end));  % age by 26 matrix 
I_OtherDSS_age = zeros(18,length(2004:y_end)); % age by 26 matrix 
R_all = zeros(18, length(2004:y_end));
S_all = zeros(18, length(2004:y_end));
%
for y= 1:length(years) %loop through years we have ABS settlement data for
    ABSYearsLabel{y}=sprintf('y%d',2003+y); %Labels to be used to store matrices
    
    CHB = ABS_CHB.(ABSYearsLabel{y}); % Get matrix of cob by age with # inf
    
    for g=1:n_cob
       
        if regions_top30_DSS(g)==1
            I_AfricaDSS_age(:,y) = I_AfricaDSS_age(:,y) + CHB(g,:)';
        elseif regions_top30_DSS(g)==2 
            I_AsiaDSS_age(:,y) = I_AsiaDSS_age(:,y) + CHB(g,:)';
        else
            I_OtherDSS_age(:,y) = I_OtherDSS_age(:,y) + CHB(g,:)';
        end
    end
    
    R_all(:,y) = sum(ABS_Res.(ABSYearsLabel{y}))'; %7 by 26
    S_all(:,y) = sum(ABS_S.(ABSYearsLabel{y}))'; %7 by 26
    
end

%R_test = [R_all(:,1:26)', R_all_other(:,1:26)'];

% Add the total in R and S from 'other COB' from DSS:
R_all = R_all + R_all_other;
S_all = S_all + S_all_other;

I_AfricaDSS_age = I_AfricaDSS_age+ I_AfricaOther;
I_AsiaDSS_age = I_AsiaDSS_age + I_AsiaOther;
I_OtherDSS_age = I_OtherDSS_age+ I_OtherOther;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Disease phase by regions:
% Distribute the total number with chronic infection across the 4 phases by
% age

% June 12th update uses BC_estimates2 (prev = BC_estimates.csv)
% Ages 0 to 9; 10 to 19; 20 to 29; 30 to 39; 40 to 49; 50 to 59; 60+ 
phase_dist_Africa = params.phase_dist_Africa_2004onwards; % dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[2 1 8 4]);
phase_dist_Asia = params.phase_dist_Asia_2004onwards; %dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[12 1 18 4]);
%phase_dist_US = dlmread('Inputs/Phase_a_gram_BC_estimates.csv',',',[22 1 28 4])/100;
% Age group 60+ for Indigenous not known:
%phase_dist_Indigenous = dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[32 1 37 4]);
phase_dist_combined = params.phase_dist_combined_2004onwards; %dlmread('Inputs/Phase_a_gram_BC_estimates2.csv',',',[45 1 51 4]);

% Assign each COB a region.

% I_Africa = age by 26yrs matrix
% I_Asia
% I_other

I_all=zeros(18,length(2004:y_end));
I1_all=zeros(18,length(2004:y_end));
I2_all=zeros(18,length(2004:y_end));
I3_all=zeros(18,length(2004:y_end));

for i=1:length(2004:y_end)
I_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,1) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,1) + I_OtherDSS_age(:,i).*phase_dist_combined(:,1));
I1_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,2) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,2) + I_OtherDSS_age(:,i).*phase_dist_combined(:,2)); % sai
I2_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,3) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,3) + I_OtherDSS_age(:,i).*phase_dist_combined(:,3));
I3_all(:,i) = (I_AfricaDSS_age(:,i).*phase_dist_Africa(:,4) + I_AsiaDSS_age(:,i).*phase_dist_Asia(:,4) + I_OtherDSS_age(:,i).*phase_dist_combined(:,4));
end


% To ages: 0-4;5-9;  10-14;15-19;  20-24;25-29;  30-34;35-39; 40-44;45-49;
% 50-54;55-59;   60-64;65-69;70-74;75-79;80-84;85+
mig_Sa = S_all; 
mig_Ia = I_all; 
mig_I1a = I1_all; 
mig_I2a = I2_all; 
mig_I3a = I3_all; 
mig_Ra = R_all; 


end



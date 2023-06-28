function state = Nat_init_state(params)
% state = Nat_init_state(params) contains initial conditions for population.

%S0=zeros(params.nA,1); %col vector prop of N0 in susceptible class in age group params.g
V0=zeros(params.nA,1);
A0=zeros(params.nA,1);
%I0=zeros(params.nA,1); %prop of N0 in Immune Tolerant to HBV due to vacc in age group params.g
%I10=zeros(params.nA,1);
%I20=zeros(params.nA,1);
%I30=zeros(params.nA,1);
C10=zeros(params.nA,1);
C20=zeros(params.nA,1);
C30=zeros(params.nA,1);
I1T0=zeros(params.nA,1);
I2T0=zeros(params.nA,1);
I3T0=zeros(params.nA,1);
C1T0=zeros(params.nA,1);
C2T0=zeros(params.nA,1);
C3T0=zeros(params.nA,1);
H0=zeros(params.nA,1);
HT0=zeros(params.nA,1);
D0=zeros(params.nA,1);
DT0=zeros(params.nA,1);
%R0=zeros(params.nA,1); %prop of N0 with cleared HBV in age group params.g

%params.population(1)*params.AgeDist1951 %gives number of people in each age group in 1951. %18x1
S0=params.population(1)*params.AgeDist(1,:).*params.init0(1,:);
I0=params.population(1)*params.AgeDist(1,:).*params.init0(2,:);
I10=params.population(1)*params.AgeDist(1,:).*params.init0(3,:);
I20=params.population(1)*params.AgeDist(1,:).*params.init0(4,:);
I30=params.population(1)*params.AgeDist(1,:).*params.init0(5,:);
R0=params.population(1)*params.AgeDist(1,:).*params.init0(6,:);

%Initial values for disease states:
state(1:params.nA)=S0; % Susceptible
state(params.nA+1:2*params.nA)=V0; % Vaccinated/Immune
state(2*params.nA+1:3*params.nA)=A0; % Acute Infection
state(3*params.nA+1:4*params.nA)=I0; % Immune Tolerant 
state(4*params.nA+1:5*params.nA)=I10; % Immune Clearance
state(5*params.nA+1:6*params.nA)=I20; % Immune Control
state(6*params.nA+1:7*params.nA)=I30; % Immune Escape
state(7*params.nA+1:8*params.nA)=C10; % Cirrhotic Immune Clearance
state(8*params.nA+1:9*params.nA)=C20; % Cirrhotic Immune Control
state(9*params.nA+1:10*params.nA)=C30; % Cirrhotic Immune Escape
state(10*params.nA+1:11*params.nA)=I1T0; % Immune Clearance on Treatment
state(11*params.nA+1:12*params.nA)=I2T0; % Immune Control on Treatment
state(12*params.nA+1:13*params.nA)=I3T0; % Immune Escape on Treatment
state(13*params.nA+1:14*params.nA)=C1T0; % Cirrhotic Immune Clearance on Treatment
state(14*params.nA+1:15*params.nA)=C2T0; % Cirrhotic Immune Control on Treatment
state(15*params.nA+1:16*params.nA)=C3T0; % Cirrhotic Immune Escape on Treatment
state(16*params.nA+1:17*params.nA)=H0;  % Hepatocellular Carcinoma
state(17*params.nA+1:18*params.nA)=HT0; % Hepatocellular Carcinoma on Treatment
state(18*params.nA+1:19*params.nA)=D0;  % Decompensated Cirrhosis
state(19*params.nA+1:20*params.nA)=DT0; % Decompensated Cirrhosis on Treatment
state(20*params.nA+1:21*params.nA)=R0;  % Recovered/Resolved HBV infection

state(21*params.nA+1:22*params.nA)=zeros(params.nA,1); %# new acute infections at time t
state(22*params.nA+1:23*params.nA)=zeros(params.nA,1); %# new Im Tol infections at time t
state(23*params.nA+1:24*params.nA)=zeros(params.nA,1); %# HCC deaths
state(24*params.nA+1:25*params.nA)=zeros(params.nA,1); %# HCC deaths w Treatment
state(25*params.nA+1:26*params.nA)=zeros(params.nA,1); %# DC deaths
state(26*params.nA+1:27*params.nA)=zeros(params.nA,1); %# DC deaths w Treatment
state(27*params.nA+1:28*params.nA)=zeros(params.nA,1); %# bg deaths
state(28*params.nA+1:29*params.nA)=zeros(params.nA,1); %# acute deaths
state(29*params.nA+1:30*params.nA)=zeros(params.nA,1); % HCC_incidence
state(30*params.nA+1:31*params.nA)=zeros(params.nA,1); % HCC_T_incidence
state(31*params.nA+1:32*params.nA)=zeros(params.nA,1); % DC_incidence
state(32*params.nA+1:33*params.nA)=zeros(params.nA,1); % DC_T_incidence
state(33*params.nA+1:34*params.nA)=zeros(params.nA,1); % treat uptake cumulative
% 

end
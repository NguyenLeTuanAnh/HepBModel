function dists = default_dists_HBV(params)
% dists = def_dists_pertussis(params)sets up the default parameter
% distributions for the HBV LHS experiment
%dists =struct();

% DISEASE PARAMETERS

% Proportion of migrants with cirrhosis
%dists.dis_phase_cirr1 = dist_beta(params.dis_phase_cirr1+0.001,0.5*(params.dis_phase_cirr1+0.001) ,0.001,1); %prop cirrhotic  <20yrs

%% Add back in:
dists.dis_phase_cirr2 = dist_beta(params.dis_phase_cirr2,0.5*params.dis_phase_cirr2 ,0.001,1); %prop cirrhotic  20-40yrs
dists.dis_phase_cirr3 = dist_beta(params.dis_phase_cirr3,0.5*params.dis_phase_cirr3 ,0.001,1); %prop cirrhotic  >40yrs

% Force of infection & acute hbv rates
dists.foi_scale = dist_flat(1,5); %dist_beta(1,0.8,0.9,100); *factor to multiply notifs by (baseline assumes 2x notifications)
dists.ac_duration = dist_flat(2,5); %duration of acute infection. 

% Risk of acute infection by age (Acute - Imm tolerant)
params.sig = num2cell(params.sig);
dists.sig{1} = dist_flat(0.2,0.95);
% for i =2:3
% dists.sig{i} = dist_flat(0.1, 0.25);
% end
% for i=4:9
% dists.sig{i} = dist_flat(0.03,0.1);
% end
% for i = 10:params.nA
% dists.sig{i} = dist_flat(0.01,0.06);
% end

% sig ranges for publication:
dists.sig{2} = dist_flat(0.1, 0.25);
dists.sig{3} = dist_flat(0.1, 0.25);
dists.sig{4} = dist_flat(0.03,0.15);
dists.sig{5} = dist_flat(0.03,0.1);
dists.sig{6} = dist_flat(0.03,0.08);
dists.sig{7} = dist_flat(0.02,0.08);
dists.sig{8} = dist_flat(0.01,0.07);
dists.sig{9} = dist_flat(0.01,0.06);
dists.sig{10} = dist_flat(0.01,0.05);
dists.sig{11} = dist_flat(0.01,0.05);
dists.sig{12} = dist_flat(0.008,0.04);
dists.sig{13} = dist_flat(0.008,0.04);
dists.sig{14} = dist_flat(0.008,0.03);
dists.sig{15} = dist_flat(0.008,0.03);
dists.sig{16} = dist_flat(0.006,0.02);
dists.sig{17} = dist_flat(0.006,0.02);
dists.sig{18} = dist_flat(0.004,0.02);


% Rate cirrhotic cases that revert back to non-cirrhotic
dists.psi_fact1 = dist_flat(0.13,0.16); 
dists.psi_fact2 = dist_flat(0.001,0.005); 

%% Estimates with capacity to breakdown by age: 
% Disease progression parameter estimates
% params.alpha2 = num2cell(params.alpha2);
% params.alpha3 = num2cell(params.alpha3);
% params.alpha4 = num2cell(params.alpha4);
% params.alpha5 = num2cell(params.alpha5);
% params.alpha6 = num2cell(params.alpha6);
% params.alpha7 = num2cell(params.alpha7);
% params.alpha9 = num2cell(params.alpha9);
% 
% params.gamma_a = num2cell(params.gamma_a);
% params.gamma_b = num2cell(params.gamma_b);
% 
% params.kappa1 = num2cell(params.kappa1);
% params.kappa2 = num2cell(params.kappa2);
% params.kappa3 = num2cell(params.kappa3);
% params.kappa4 = num2cell(params.kappa4);
% params.kappa5 = num2cell(params.kappa5);
% params.kappa6 = num2cell(params.kappa6);
% params.kappa7 = num2cell(params.kappa7);
% params.kappa8 = num2cell(params.kappa8);
% params.kappa9 = num2cell(params.kappa9);
% params.kappa10 = num2cell(params.kappa10);
% params.kappa11 = num2cell(params.kappa11);
% params.kappa12 = num2cell(params.kappa12);
% params.kappa13 = num2cell(params.kappa13);
% params.kappa14 = num2cell(params.kappa14);
% params.kappa15 = num2cell(params.kappa15);
% 
% params.omega2 = num2cell(params.omega2);
% params.omega3 = num2cell(params.omega3);
% params.omega6 = num2cell(params.omega6);
% 
% params.l1 = num2cell(params.l1);

% for i=1:params.nA
%     dists.alpha2{i} = dist_flat(0.8*params.alpha2{i},1.2*params.alpha2{i});
%     dists.alpha3{i} = dist_flat(0.8*params.alpha3{i},1.2*params.alpha3{i});
%     dists.alpha4{i} = dist_flat(0.8*params.alpha4{i},1.2*params.alpha4{i});
%     dists.alpha5{i} = dist_flat(0.8*params.alpha5{i},1.2*params.alpha5{i});
%     dists.alpha6{i} = dist_flat(0.8*params.alpha6{i},1.2*params.alpha6{i});
%     dists.alpha7{i} = dist_flat(0.8*params.alpha7{i},1.2*params.alpha7{i});
%     dists.alpha9{i} = dist_flat(0.8*params.alpha9{i},1.2*params.alpha9{i});
%     
%     dists.gamma_a{i} = dist_flat(0.8*params.gamma_a{i},1.2*params.gamma_a{i});
%     dists.gamma_b{i} = dist_flat(0.8*params.gamma_b{i},1.2*params.gamma_b{i});
%     
%     dists.kappa1{i} = dist_beta(params.kappa1{i},0.5*params.kappa1{i},0,1); 
%     dists.kappa2{i} = dist_beta(params.kappa2{i},0.5*params.kappa2{i},0,1);
%     dists.kappa3{i} = dist_beta(params.kappa3{i},0.5*params.kappa3{i},0,1);
%     dists.kappa4{i} = dist_beta(params.kappa4{i},0.5*params.kappa4{i},0,1);
%     dists.kappa5{i} = dist_beta(params.kappa5{i},0.5*params.kappa5{i},0,1);
%     dists.kappa6{i} = dist_beta(params.kappa6{i},0.5*params.kappa6{i},0,1);
%     dists.kappa7{i} = dist_beta(params.kappa7{i},0.5*params.kappa7{i},0,1);
%     dists.kappa8{i} = dist_beta(params.kappa8{i},0.5*params.kappa8{i},0,1);
%     dists.kappa9{i} = dist_beta(params.kappa9{i},0.5*params.kappa9{i},0,1);
%     dists.kappa10{i} = dist_beta(params.kappa10{i},0.5*params.kappa10{i},0,1);
%     dists.kappa11{i} = dist_beta(params.kappa11{i},0.5*params.kappa11{i},0,1);
%     dists.kappa12{i} = dist_beta(params.kappa12{i},0.5*params.kappa12{i},0,1);
%     dists.kappa13{i} = dist_beta(params.kappa13{i},0.5*params.kappa13{i},0,1);
%     dists.kappa14{i} = dist_beta(params.kappa14{i},0.5*params.kappa14{i},0,1);
%     dists.kappa15{i} = dist_beta(params.kappa15{i},0.5*params.kappa15{i},0,1);
%     
%     dists.omega2{i} = dist_flat(0.8*params.omega2{i},1.2*params.omega2{i});
%     dists.omega3{i} = dist_flat(0.8*params.omega3{i},1.2*params.omega3{i});
%     dists.omega6{i} = dist_flat(0.8*params.omega6{i},1.2*params.omega6{i});
%     
%     dists.l1{i} = dist_flat(0.8*params.l1{i}, 1.2*params.l1{i});
%     
% end


% Progression to cirrhosis:
dists.omega1 = dist_flat(0.013, 0.019);
dists.omega4 = dist_flat(0.0061, 0.012);


% Progression to Decompensated Cirrhosis (not stratefied by age)
%dist_fn = dist_beta(mean, variance, min_v, max_v);
dists.theta1 = dist_beta(params.theta1,0.5*params.theta1 ,0.2*params.theta1,0.046);
dists.theta2 = dist_beta(params.theta2,0.5*params.theta2 ,0.5*params.theta2,0.02);
dists.theta3 = dist_flat(0.00726,0.0294);
dists.theta4 = dist_flat(0.0022,0.0089);



% Rate of going off treatment:
dists.eta_scale = dist_flat(0.1,0.2);
dists.eta_scale2 = dist_flat(0.05,0.15);


% Mortality due to HBV
dists.l2 = dist_flat(0.25,0.4);
dists.l3 = dist_flat(0.18,0.3);
dists.l4 = dist_flat(0.16,0.314);

%%%%

%dist_norm(params.alpha(:,2))

% dists.FoI_mult = dist_flat(params.FoI_mult*0.8, params.FoI_mult*5);
% dists.ChronicL = dist_beta(0.35, 0.2, 0, 2); %pt_est=0.005; 
% dists.ChronicI = dist_beta(3.8, 1.1, 2, 8); %pt est = 0.05;
% dists.ChronicH = dist_normal(11,0.9); %pt_est=0.10.

%Progression rates through HBV phases
% dists.p_alpha = dist_flat(0.1*params.p_alpha,2*params.p_alpha);
 
% dists.ms_mult = dist_flat(0.9*params.ms_mult,1.2*params.ms_mult);
% dists.mi_mult = dist_flat(0.9*params.mi_mult,1.1*params.mi_mult);
% dists.mi1_mult = dist_flat(0.9*params.mi1_mult,1.1*params.mi1_mult);
% dists.mi2_mult = dist_flat(0.9*params.mi2_mult,1.1*params.mi2_mult);
% dists.mi3_mult = dist_flat(0.9*params.mi3_mult,1.1*params.mi3_mult);


%dists.kappa_scale = dist_flat(0.9,1.1);
%dists.omega_scale = dist_flat(0.9,1.2); 
%dists.lambda_scale = dist_const(0.8); 
% dists.mr_mult = dist_flat(8,12);

%% Migration Params:
dists.DSS_MigPrev = dist_flatint(4,6); %choose between Lower CI (4), Upper CI (5) and Jenns for prev ests in Top30 COB
dists.DSS_MigPrev2019 = dist_flatint(7,9);  %2019 syst review update. choose between Lower CI (8), Upper CI (9) and Jenns for prev ests in Top30 COB

dists.NOM_Series = dist_flatint(0,2);
% params.prev_COB_1951to1974 = num2cell(params.prev_COB_1951to1974);
% params.prev_COB_1975to1990 = num2cell(params.prev_COB_1975to1990);
% 
 %dist_beta(params.theta2,0.5*params.theta2 ,0.5*params.theta2,0.02);
%  for ii=1:length(params.prev_COB_1951to1974)
% %params.prev_COB_1951to1974{ii}
%  dists.prev_COB_1951to1974{ii} = dist_beta(params.prev_COB_1951to1974{ii}, params.prev_COB_1951to1974{ii}, 0.5*params.prev_COB_1951to1974{ii}, 1.5*params.prev_COB_1951to1974{ii});
%  end
%  for jj=1:length(params.prev_COB_1975to1990)
%  dists.prev_COB_1975to1990{jj} = dist_beta(params.prev_COB_1975to1990{jj}, params.prev_COB_1975to1990{jj}, 0.5*params.prev_COB_1975to1990{jj}, 1.5*params.prev_COB_1975to1990{jj});
%  end


%% Add back in:
 dists.low_prev = dist_beta(params.low_prev, 0.001, 0.001, 1); 
 dists.int_prev = dist_beta(params.int_prev, 0.5*params.int_prev, 0.001, 1); 
 dists.high_prev = dist_beta(params.high_prev, 0.5*params.high_prev, 0.001, 1); 


% Think about how we could do uncertainty on phase distribution too...



end




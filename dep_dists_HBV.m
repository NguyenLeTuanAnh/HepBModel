function dists = dep_dists_HBV(params)
% dists = dep_dists_pertussis(params) calculates the distributions of
% parameters which are dependent on other parameters

dists=struct();

%dists.MIG_S = 1 - (dists.MIG_I + dists.MIG_I1 + dists.MIG_I2 + dists.MIG_I3 + dists.MIG_R);


% Allow prog rates to vary up/down by p_alpha between +/-10% %
% dists.alpha(:,2)= dist_const([ones(2,1)*0.03; ones(4,1)*0.02; ones(2,1)*0.03; ones(10,1)*0.06]./12);  % I -> I1 **
% dists.alpha(:,2)= dist_const( params.p_alpha*([ones(2,1)*0.03; ones(4,1)*0.02; ones(2,1)*0.03; ones(2,1)*0.06; 0.07*ones(2,1); ones(6,1)*0.08]./12) );  % I -> I1 **
% dists.alpha(:,3)= dist_const( params.p_alpha*(1.1*[0.05; 0.12*ones(3,1); 0.184*ones(3,1); 0.136*ones(3,1); 0.15*ones(8,1) ]./12) ); % I1 -> I2  
% %dists.alpha(:,3)= [ones(4,1)*0.05; ones(4,1)*0.1; ones(10,1)*0.15]./12; % I1 -> I2 
% dists.alpha(:,4)= dist_const( params.p_alpha*( 0.8*[0.0087*ones(6,1); 0.0149*ones(2,1); 0.0278*ones(2,1); 0.0202*ones(8,1)]./12) ); % I2 -> I3 
% dists.alpha(:,5)= dist_const( params.p_alpha*([0.05; 0.12*ones(3,1); 0.184*ones(3,1); 0.136*ones(3,1); 0.15*ones(8,1) ]./12) ); % C1 -> C2 1.1*
% dists.alpha(:,6)= dist_const( params.p_alpha*([0.0087*ones(6,1); 0.0149*ones(2,1); 0.0278*ones(2,1); 0.0202*ones(8,1)]./12) );  % C2 -> C3 
% dists.alpha(:,7)= dist_const( params.p_alpha*(1.1*[0.11; 0.2664*ones(3,1); 0.4085*ones(3,1); 0.3019*ones(3,1); ones(8,1)*0.33]./12) ); % I1T -> I2T
% %dists.alpha(:,8)= zeros(18,1)./12; % I2T -> I3T
% dists.alpha(:,9)= dist_const( params.p_alpha*([0.11; 0.2664*ones(3,1); 0.4085*ones(3,1); 0.3019*ones(3,1); ones(8,1)*0.33]./12) );% C1T -> C2T 1.1*
% dists.alpha(:,10)= dist_const( params.p_alpha*((zeros(18,1)./12) ) ); % C2T -> C3T
% 

end
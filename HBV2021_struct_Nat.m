function model=HBV2021_struct_Nat()
% Age-structured HBV model in format suitable for Latin Hypercube
% Sampling

% Written by Trish Campbell, adapted by Karen McCulloch & Anh Nguyen for HBV 

% Collect the various parts of the model into a single structure
model = struct(...
    'params', @params, 'unsaved_params',@unsaved_params,'state', @state, 'events', @init_events,'statistics', @stats,...
    'default_dists', @def_dists, 'dependent_dists', @dep_dists,...
    'equations', @rhs...
    );
end

%% PARAMS
function params=params()

%Define a parameter file for: National, TAS, VIC, ACT, NSW, QLD, NT, WA, SA
% and National-ATSI.
params = Nat_Params(); %setup_HBV2();
%params = NSW_Params();


end

%% UNSAVED_PARAMS
function unsaved_params = unsaved_params()

unsaved_params = {'age_groups'};

end
    
%% STATE
function state = state(params)


%% Specifiy which state/territory or national model you are running:
state = Nat_init_state(params);
%state = NSW_init_state(params);


end

%% STATS
function output = stats(time,params,state,events)

%% Specifiy which state/territory or national model you are running:
if nargin == 2
    output = Nat_stats(time,params);
    %output = NSW_stats(time,params);
return;

elseif nargin ~= 4
    fprintf(1,'ERROR: invalid # parameters to statistics\n');
    return;
end

%% Specifiy which state/territory or national model you are running:
output = Nat_stats(time,params,state,events);

end

%% DEFAULT DISTRIBUTIONS
function dists = def_dists(params)

dists = default_dists_HBV(params);

end

%% DEPENDENT DISTRIBUTIONS
function dists = dep_dists(params)

dists = dep_dists_HBV(params);

end

%% EVENTS
function events = init_events(initial_state,params)

events = init_events_HBV(initial_state,params);

end

%% EQUATIONS
function rates = rhs(time,state,params)

rates = HBV2021_Eqs_report(time,state,params); %stable tx 2021/2022 then inc based on 2016 - 2019 average increase
% rates = HBV2021_Eqs_report_opt_tx(time,state,params); %stable tx 2021/2022 then inc to reach GHSS 2030 target.

%rates = HBV2021_Eqs_stablefuturetx(time,state,params); %future tx numbers remain stable


end

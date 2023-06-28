%SOLVE_DISCRETELY   Solve a stochastic model using Euler's method.
%   Given a set of model equation, an initial model state and a set of model
%   parameters, SOLVE_DISCRETELY uses fixed-step solvers (eg, ODE1) to solve
%   the model equations. The fixed-step solvers are generally much slower than
%   the built-in MATLAB ODE solvers used by SOLVE_CONTINUOUSLY, which should
%   always be used for deterministic models.
%
%   [time state] = SOLVE_DISCRETELY(rhs, initial_state, params, solver, change_fn);
%
%   RHS is a function that calculates the rate of change for each state
%   variable in the model, given the time, the current state of the model and
%   the values of the model parameters.
%
%   INITIAL_STATE is a vector that specifies the initial value of each state
%   variable in the model.
%
%   PARAMS is a record (structure array) that defines the value of each
%   independent and dependent  model parameter.
%
%   SOLVER is an optional function handle that specifies which fixed-step
%   solver to use. The default value is @ode1_nneg; valid options are
%   @ode1_nneg, @ode1, @ode2, @ode3, @ode4 and @ode5.
%
%   CHANGE_FN is an optional function handle that is called to make
%   arbitrary changes to the state and parameters at the times specified in
%   params.time_splits. It must be defined if params.time_splits is
%   defined. Furthermore, if you want to specify CHANGE_FN, you must also
%   specficy SOLVER must be defined for technical reasons.
%
%   TIME is the vector of time steps taken by the solver.
%
%   STATE is a matrix that contains the evolution of each state variable over
%   the time steps specified by TIME.
%
%   Example:
%      % Use the ode2 solver instead of the (default) ode1_nneg solver
%      >> [time state] = solve_discretely(rhs, init_state, params, @ode2);
%
%See also solve_continuously, solve, ode1
function [steps, final_state] = solve_discretely(rhs, initial_state, params, ...
    solver, change_fn)
    global events;
    
    if nargin < 4
        solver = @ode1_nneg; % _nneg options: ode1_nneg, ode1, ode2, ode3, ode4, ode5
    end
    
    % The random stream used for the stochastic model
    %events.stochStream = RandStream('mt19937ar', 'Seed', 2010);
    events.stochStream = RandStream.getGlobalStream;
    oldStream = RandStream.setGlobalStream(events.stochStream);
    
    % Size of time steps
    events.timescale = ...
        params.time_entries / (params.time_end-params.time_start);
    
    steps = 0:1:params.time_entries;
    
    if ~isfield(params,'time_splits')
        
        segment_steps = cell(1,1);
        % Simulation time steps
        segment_steps{1} = 0:1:params.time_entries;
        
    else
        
        % First integration period
        % We need to evenly split time_entries across the
        key_times = [params.time_start params.time_splits' params.time_end];
        
        segment_steps = cell(1,length(key_times)-1);
        
        lengths_of_segments = key_times(2:end)-key_times(1:end-1);
        steps_per_day = events.timescale;
        
        segment_steps{1} = 0:1:lengths_of_segments(1)*steps_per_day;
        for i = 2:(length(key_times)-1)
            % The 2nd (and further) segment steps begin at the same time as
            % the preceeding segment finished so that the initial state
            % is the last recorded state in the previous segment.
            % We write states 2:end into the final_state matrix.
            segment_steps{i} = segment_steps{i-1}(end):1:segment_steps{i-1}(end)+lengths_of_segments(i)*steps_per_day;
        end
        
    end
    
    % Solve the equations and return the state evolution
    do_convert = isstruct(initial_state);
    if do_convert
        % Convert between vector and record layouts
        layout = layout_of(initial_state);
        initial_state = layout_vec(layout, initial_state);
    end
    
    final_state = zeros(length(steps),length(initial_state));
    final_state(1,:) = initial_state; % Write the initial state into the
                                      % first row
    
    % Loop over each sub-part of the time integration
    for i = 1:length(segment_steps)
    
        if do_convert
            % Convert between vector and record layouts
            rhs_fn = @(time, values) ...
                layout_vec(layout, rhs(time, layout_rec(layout, values), params));
        else
            rhs_fn = @(time, values) rhs(time, values, params);
        end
        
        returned_states = solver(rhs_fn, segment_steps{i}, initial_state);
        final_state(segment_steps{i}(2:end)+1,:) = returned_states(2:end,:);
        
        initial_state = returned_states(end,:)';
        
        % Call the "change_state" function here to do funky things to the
        % initial state and the parameters for the next time segment.
        if isfield(params,'time_splits') && (i < length(segment_steps))
            [initial_state params] = change_fn(segment_steps{i}(end), initial_state, params);
        end
    end
    
    if do_convert
        % Convert the final state from vector layout to record layout
        final_state = layout_time_rec(layout, final_state);
    end
    
    % Restore the previous random stream
    RandStream.setGlobalStream(oldStream);
end

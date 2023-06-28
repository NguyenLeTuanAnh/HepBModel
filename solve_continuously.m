%SOLVE_CONTINUOUSLY   Solve a deterministic model using a continuous solver.
%   Given a set of model equation, an initial model state and a set of model
%   parameters, SOLVE_CONTINUOUSLY uses built-in MATLAB ODE solvers (eg, ODE23)
%   to solve the model equations. The MATLAB ODE solvers are generally much
%   faster than the fixed-step (Euler method) solvers used by SOLVE_DISCRETELY,
%   which should only be used for stochastic models.
%
%   [time state] = SOLVE_CONTINUOUSLY(rhs, initial_state, params, solver);
%
%   RHS is a function that calculates the rate of change for each state
%   variable in the model, given the time, the current state of the model and
%   the values of the model parameters.
%
%   INITIAL_STATE is a vector that specifies the initial value of each state
%   variable in the model.
%
%   PARAMS is a record (structure array) that defines the value of each
%   independent and dependent model parameter.
%
%   SOLVER is an optional function handle that specifies which ODE solver to
%   use; the default value is @ode23. See the MATLAB documentation for a list
%   of ODE solvers.
%
%   TIME is the vector of times at which the model state was recorded.
%
%   STATE is a matrix that contains the evolution of each state variable over
%   the time period specified by TIME.
%
%   Example:
%      % Use the ode45 solver instead of the (default) ode23 solver
%      >> [time state] = solve_continuously(rhs, init_state, params, @ode45);
%
%See also solve_discretely, solve, ode23
function [time, final_state] = solve_continuously(rhs, initial_state, params, solver)
    if nargin < 4
        solver = @ode45;%@ode45;
    end

    % Set any specified options for the ODE solver
    %options = odeset(varargin{:});
%     if ~ any( strcmp({'NonNegative'}, varargin(1:2:end-1)) )
%        % Ensure that all state variables are kept non-negative,
%        % unless the 'NonNegative' option has already been set
%        state_vec = 1:length(initial_state);
%        options = odeset(options, 'NonNegative', state_vec);
%     end

    % Ensure that all state variables are kept non-negative
    state_vec = 1:length(initial_state);
    options = odeset('NonNegative', state_vec, 'RelTol',1.0000e-08);

    % Simulation time steps
    steps = linspace(params.time_start, params.time_end, params.time_entries);

    % Solve the equations and return the state evolution
    do_convert = isstruct(initial_state);
    if do_convert
        % Convert between vector and record layouts
        layout = layout_of(initial_state);
        rhs_fn = @(time, values) ...
            layout_vec(layout, rhs(time, layout_rec(layout, values), params));
        initial_state = layout_vec(layout, initial_state);
    else
        rhs_fn = @(time, values) rhs(time, values, params);
    end
    [time, final_state] = solver(rhs_fn, steps, initial_state, options);
    if do_convert
        % Convert the final state from vector layout to record layout
        final_state = layout_time_rec(layout, final_state);
    end
end

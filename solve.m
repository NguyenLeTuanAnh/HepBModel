%SOLVE   Perform a single simulation of some model.
%   Given a defined model and values for the independent and dependent model
%   parameters, SOLVE uses SOLVE_CONTINUOUSLY or SOLVE_DISCRETELY to solve the
%   model equations, depending on whether the parameter 'solve_discretely' is
%   set to 0 or 1, respectively. If this paramter is not specified, the default
%   solver is SOLVE_CONTINUOUSLY.
%
%   [time state events] = SOLVE(model, params);
%
%   MODEL is a record (structure array) that contains the model details (see
%   EXAMPLE_SIR_MODEL for an example of such a model).
%
%   PARAMS is a record (structure array) that defines the value of each
%   independent and dependent model parameter.
%
%   TIME is the vector of times at which the model state was recorded.
%
%   STATE is a matrix that contains the evolution of each state variable over
%   the time period specified by TIME.
%
%   EVENTS is a private copy of the global variable events, which records any
%   events of interest that occurred during the simulation.
%
%See also solve_continuously, solve_discretely, single_simulation
function [time, state, my_events] = solve(model, params)
    global events;

    % initialise the model state
    initial_state = model.state(params);
    % initialise the global events record, to track simulation events
    if isfield(model, 'events')
        events = model.events(initial_state, params);
    else
        events = struct();
    end
    
    % generate the once-off unstored additional parameters
    if isfield(model, 'gen_params')
        additional_params = model.gen_params(params);
        
        names = fieldnames(additional_params);
        for i = 1:length(names)
            if isfield(params, names(i))
                error('solve:invalid', ...
                    '"%s" already defined in params.', names{i})
            else
                params.(names{i}) = additional_params.(names{i});
            end
        end
    end

    rhs = model.equations;
    if isfield(params, 'solve_discretely') && params.solve_discretely
        % Fixed step solver
        if isfield(model,'change_s_and_p')
            change_fn = model.change_s_and_p;
            [time, state] = solve_discretely(rhs, initial_state, params, ...
                @ode1_nneg, change_fn);
        else
            [time, state] = solve_discretely(rhs, initial_state, params);
        end
        time = time+params.time_start;
    else
        % ODE solver
        if isfield(params,'time_splits') || isfield(model,'change_s_and_p')
            fprintf(1,'\n\nWARNING: The model has defined split times and/or the change_s_and_p function, but is using the continuous solver. The splits will be ignored.\n\n')
        end
        
        params.solve_discretely = 0;
        [time, state] = solve_continuously(rhs, initial_state, params, @ode45); %ode45 orig @ode15s
    end

    % return a copy of the global events record
    my_events = struct(events);
end

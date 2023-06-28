%SINGLE_SIMULATION   Perform a single simulation of some model.
%   Given a defined model and values for the independent model parameters,
%   SINGLE_SIMULATION selects values for the dependent parameters and passes
%   the complete set of parameters to SOLVE, in order to solve the model
%   equations.
%
%   [time params state events] = SINGLE_SIMULATION(model, parameters, rs);
%
%   MODEL is a record (structure array) that contains the model details (see
%   EXAMPLE_SIR_MODEL for an example of such a model).
%
%   PARAMETERS is a record (structure array) that defines the values of each
%   independent model parameter.
%
%   RS is an (optional) RandStream object that is used to randomly select the
%   value of each dependent parameter; the default value is to use the default
%   stream (RandStream.getDefaultStream).
%
%   TIME is the vector of times at which the model state was recorded.
%
%   PARAMS is a record (structure array) of the model parameter values used for
%   this simulation, including the randomly-chosen dependent parameters.
%
%   STATE is a matrix that contains the evolution of each state variable over
%   the time period specified by TIME.
%
%   EVENTS is a private copy of the global variable events, which records any
%   events of interest that occurred during the simulation.
%
%See also example_SIR_model, solve, solve_continuously, solve_discretely
function [time, params, state, events] = single_simulation(model, parameters, rs)
    % use the provided RandStream, if provided, otherwise use the default stream
    if nargin < 3
        rs = RandStream('mt19937ar', 'Seed', 3001);
    end
    % set the value of all dependent parameters
    if isfield(model, 'dependent_dists')
        dep_dists = model.dependent_dists(parameters);
        dep_names = fieldnames(dep_dists);
        for param = 1:length(dep_names)
            pname = dep_names{param};
            pdist = dep_dists.(pname);
            % choose each value at random from the given distribution          
            if isa(pdist,'function_handle')
                parameters.(pname) = pdist(rs.rand());
            else
                % Handle vector and matrix parameter distributions
                for j = 1:size(pdist,2)
                    for k = 1:size(pdist,1)
                        parameters.(pname)(k,j) = pdist{k,j}(rs.rand());
                    end
                end
            end   
        end
    end
    % perform the simulation and return the results
    [time, state, events] = solve(model, parameters);
    params = struct(parameters);
end

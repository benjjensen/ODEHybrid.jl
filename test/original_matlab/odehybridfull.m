function varargout = odehybridfull(solver, ode, de, dt, ...
                                   ts, yc0, yd0, varargin)

    % Let the user pass in cells or anything else, but always make sure we
    % work with cells. This simplifies life when we have to dynamically
    % pass these into functions.
    if ~iscell(yc0), yc0 = {yc0}; end;
    if ~iscell(yd0), yd0 = {yd0}; end;
    if ~iscell(de),  de  = {de};  end;
    
    % Create the initial continuous state vector.
    yc0v = state_to_vector(yc0);

    % Create a function to expand the vector into the state, pass the state
    % to the ODE, and turn the result back into a vector.
    ode_v = @(t,yc,yd,varargin) run_ode(ode, t, yc, yd, yc0, varargin{:});
    for k = 1:length(dt)
        de_v{k} = @(t, yc, yd, varargin) run_de(de{k}, t, yc, yd, yc0, ...
                                                varargin{:});
    end

    % Determine what outputs we need.
    if nargout == 1                                     % Structure
        n_outputs = 8;
    elseif nargout <= 1 + numel(yc0)                    % Continuous
        n_outputs = 2;
    elseif nargout <=   1 + numel(yc0) ...              % Continuous
                      + 1 + numel(yd0)                  % Discrete
        n_outputs = 4;
    elseif nargout <=   1 + numel(yc0) ...              % Continuous
                      + 1 + numel(yd0) ...              % Discrete
                      + 1 + numel(yc0) + numel(yd0) + 1 % Events
        n_outputs = 8;
    else
        error('Too many outputs are requested from odoehybrid.');
    end
    
    % Propagate the vector versions of the ODE and DE.
    outputs = cell(1, n_outputs);
    [outputs{:}] = odehybridcore(solver, ...
                                 ode_v, de_v, dt, ts, yc0v, yd0, ...
                                 varargin{:});

	% For states at events, convert to the appropriate type from the big
	% cell array of stored values.
    if n_outputs >= 8
        
        % For the continuous part...
        continuous_states = cell(1, numel(yc0));
        [continuous_states{:}] = vectors_to_states(outputs{6}, yc0{:});

        % Make the output.
        outputs = [outputs(1:5), ...
                   continuous_states, ...
                   output_discrete_states(outputs{7}, yd0), ...
                   outputs(8:end)];
        
    end

	% For discrete states, convert to the appropriate type from the big
	% cell array of stored values.
    if n_outputs >= 4
        outputs = [outputs(1:3), ...
                   output_discrete_states(outputs{4}, yd0), ...
                   outputs(5:end)];
    end
    
    % For continuous state, convert back to original types from the state
    % vectors.
    if n_outputs >= 2
        continuous_states = cell(1, numel(yc0));
        [continuous_states{:}] = vectors_to_states(outputs{2}, yc0{:});
        outputs = [outputs(1), ...
                   continuous_states, ...
                   outputs(3:end)];
    end
    
	% Return separately the states as cell arrays.
    varargout = outputs;

    % If there's only one output, use the structure format.
    if nargout == 1
        
        sol.t   = outputs{1};
        sol.yc  = outputs(1+(1:numel(yc0)));
        sol.td  = outputs{1+numel(yc0)+1};
        sol.yd  = outputs(1+numel(yc0)+1+(1:numel(yd0)));
        sol.te  = outputs{1+numel(yc0)+1+numel(yd0)+1};
        sol.yce = outputs(1+numel(yc0)+1+numel(yd0)+1+(1:numel(yc0)));
        sol.yde = outputs(  1+numel(yc0)+1+numel(yd0)+1+numel(yc0) ...
                          + (1:numel(yd0)));
        sol.ie  = outputs{end};
        varargout{1} = sol;
    end

end
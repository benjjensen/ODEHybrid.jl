function dyvdt = run_ode(ode, t, ycv, yd, yc0, varargin)

    % Pull out the continuous state first (this is how we store it).
    yc = vector_to_state(ycv, yc0);
    
    % Get the state differences.
    [state_difference{1:numel(yc0)}] = ode(t, yc{:}, yd{:}, varargin{:});
    
    % Convert the derivative to a vector.
    dyvdt = state_to_vector(state_difference);
    
end
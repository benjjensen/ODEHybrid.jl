function [ycv, yd] = run_de(de, t, ycv, yd, yc0, varargin)
%    Formats discrete equation to be used in ODE Hybrid Core 
%       Called in ODE Hybrid Full when inputs are not pure numerics
%    
%    Inputs:
%       
%       de:    Function of form @(t, x, u)  (...)                       |  Function handle
%       t:      Time                                                                  |  scalar
%       ycv:  Values of continuous-time ode at given time     |  column vector
%       yd:    Values of discrete-time system at given time    | Cell of Array {[...]}
%       yc0:  Initial values/format for continuous-time eqn    |  Cell 
%    
   
    % Pull out the continuous state first (this is how we store it).
    yc = vector_to_state(ycv, yc0);

    % Run the update.
    [yc{:}, yd{:}] = de(t, yc{:}, yd{:}, varargin{:});
    
    % Convert to vectors.
    ycv = state_to_vector(yc);
    
end
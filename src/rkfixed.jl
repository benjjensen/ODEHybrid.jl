# [src/rkfixed.jl]

"""    
Runge-Kutta fixed-step integration using the specified weights, nodes,
and Runge-Kutta matrix (or the Runge-Kutta 4th order "3/8" method by
default).

Implements numerical propagation of an ordinary differential equation
from some initial value over the desired range. This function is similar
to MATLAB's variable-step ODE propagators (e.g., ode45), but uses a
fixed step method. This is useful either when one knows an appropriate 
step size or when a process is interrupted frequently (ode45 and the
similar functions in MATLAB will always make at least a certain number of
steps between ts(1) and ts(2), which may be very many more than are
necessary).

This function is generic for all fixed-step Runge-Kutta methods. That is,
any fixed-step Runge-Kutta propagator can be created by passing the
weightes, nodes, and Runge-Kutta matrix (together, the Butcher tableau) 
into this function. See the example below.

    t, x = rkfixed(ode, ts, x0, dt);
    t, x = rkfixed(ode, ts, x0, dt, a, b, c);
    t, x = rkfixed(ode, ts, x0, options);
    t, x = rkfixed(ode, ts, x0, options, a, b, c);

Inputs:

    ode     Ordinary differential equation function
    ts      Time span, [t_start, t_end]
    x0      Initial state (column vector)
    options One can specify an options structure instead of dt      
                so that this function is compatible with ode45 and its ilk. The
                only valid fields are MaxStep (the time step) and OutputFcn
    a       Runge-Kutta matrix
    b       Weights
    c       Nodes

Outputs:

    t       Time history
    x       State history, with each row containing the state corresponding to
                the time in the same row of t.

Example:

    # Simulate an undamped harmonic oscillator for 10s with a 0.1s time
    #    step, starting from an initial state of [1; 0] using RK 4th order
    #    integration (via the Butcher tableau specified by a, b, and c). This
    #    is exactly the same as the rk4 function

    a = [ 0    0    0  0; 
          0.5  0    0  0; 
          0    0.5  0  0; 
          0    0    1  0];

    b = [ 1  2  2  1]/6;
    c = [ 0  0.5  0.5  1];

    t, x = rkfixed( (t,x) -> [-x[2]; x[1]], [0 10], [1; 0], 0.1, a, b, c);
    plot(t, x);

See "Runge-Kutta methods" on Wikipedia for discussion of the Butcher
tableau (a, b, and c).

See also: odehybrid, ode45, odeset, rk4.

Online doc: http://www.anuncommonlab.com/doc/odehybrid/rkfixed.html

Copyright 2014 An Uncommon Lab
"""
function rkfixed(ode, ts, x0, options::ODE_SET, a, b, c)

    dt = options.MaxStep
    if isempty(dt)
        @error "RKFixed: Specify the time step with the 'MaxStep' option."
    end
    
    # Time history
    t = ts[1]:dt:ts[end]

    # Go up to the end
    if abs(t[end] - ts[end]) > abs(dt * 0.0000001) # != ts[end]  <-  '≂' led to some numerical issues
        t = [t; ts[end]]
    end
    
    ns = maximum(size(t));            # Number of samples
    nx = length(x0);                  # Number of states
    x  = [x0[:]'; zeros(ns-1, nx)];   # State history
    xk = x0[:];                       # Current state
    
    s = maximum(size(b));             # Length of weights
    d = zeros(nx, s);                 # Matrix of derivatives
        
    # # If the user provided an OutputFcn, use it.
    if (options.OutputFcn != [])
        options.OutputFcn(ts, x0, "init")
    end
        
    # Propagate.
    for k = 1:ns-1
        
        # The last sample may be cut short.
        if k == ns-1
            dt = t[k+1] - t[k];
        end
        
        # Current time
        tk = t[k];
            
        # Calculate derivatives.
        d[:, 1] = dt * ode(tk, xk);
        for z = 2:s
            dxk = d[:, 1:z-1] * a[z, 1:z-1]
            d[:, z] = dt * ode(tk + c[z] * dt, xk + dxk);
        end
        
        # Update the state.
        for z = 1:s
            xk = xk + b[z] * d[:, z];
        end
            
        # Store.
        x[k+1, :] = xk[:]'
        
        # # If the user provided an OutputFcn, use it.
        if (options.OutputFcn != [])

            if (options.OutputFcn(t[k+1], xk, ""))
                t = t[1:k+1]
                x = x[1:k+1, :]
                break
            end
        end
        
    end
    
    # If the user provided an OutputFcn, use it.
    if (options.OutputFcn != [])
        options.OutputFcn([], [], "done")
    end

    return t, x    
end



function rkfixed(ode, ts, x0, options::ODE_SET)
    """
        Wrapper function for when an RK tableau is not provided. Defaults to RK 3/8 
    """
    # Default to Runge and Kutta's 3/8 formulation (4th order).
    a = [ 0    0 0 0; 
          1/3  0 0 0; 
         -1/3  1 0 0; 
          1   -1 1 0];

    b = [1 3 3 1]/8;

    c = [0 1/3 2/3 1];

    return rkfixed(ode, ts, x0, options, a, b, c)
end

function rkfixed(ode, ts, x0, dt, a, b, c)
    """
        Wrapper function for when a time step (dt) is provided rather than an ODE SET. 
            Generates the default ODE SET and updates the MaxStep to dt
    """

    options = ODE_SET()
    options.MaxStep = dt 
    
    return rkfixed(ode, ts, x0, options, a, b, c)
end
 
function rkfixed(ode, ts, x0, dt)
    """
        Wrapper function for when no RK tableau is provided and a time step (dt) is provided rather than an ODE SET.
            Defaults to RK 3/8 and the default ODE SET with an updated MaxStep
    """

    # Default to Runge and Kutta's 3/8 formulation (4th order).
    a = [ 0    0 0 0; 
          1/3  0 0 0; 
         -1/3  1 0 0; 
          1   -1 1 0];

    b = [1 3 3 1]/8;

    c = [0 1/3 2/3 1];

    return rkfixed(ode, ts, x0, dt, a, b, c)
end


function rk4(ode, ts, x0, options::ODE_SET)
    """ rk4
     
        Runge-Kutta 4th order integration method.
    
        Implements numerical propagation of an ordinary differential equation
        from some initial value over the desired range. This function is similar
        to MATLAB's variable-step ODE propagators (e.g., ode45), but uses a
        fixed step method. This is useful either when one knows an appropriate 
        step size or when a process is interrupted frequently (ode45 and the
        similar functions in MATLAB will always make at least a certain number of
        steps between ts(1) and ts(2), which may be very many more than are
        necessary).

        t, x = rk4(ode, ts, x0, options);

        Inputs:
    
            ode     Ordinary differential equation function
            ts      Time span, [t_start, t_end]
            x0      Initial state (column vector)
            options One can specify an options structure instead ohttp://www.anuncommonlab.com/doc/odehybrid/rk4.htmlf dt
                    so that this function is compatible with ode45 and its ilk. The
                    only valid fields are MaxStep (the time step) and OutputFcn
                    (if dt is passed in rather than options, an empty ODE_SET is initialized
                    with the provided time step dt)

        Outputs:
        
            t       Time history
            x       State history, with each row containing the state corresponding to
                        the time in the same row of t.

        Example:
        
            Simulate an undamped harmonic oscillator for 10s with a 0.1s time
                step, starting from an initial state of [1; 0]:

            ode = (t, x) -> [-x[2]; x[1]]
            t, x = rk4(ode, [0 10], [1; 0], 0.1);
            plot(t, x);

        See also: odehybrid, ode45, odeset, rkfixed, rkadapt.

        Online doc: http://www.anuncommonlab.com/doc/odehybrid/rk4.html

        Copyright 2014 An Uncommon Lab
    """

    dt = options.MaxStep 
    if isempty(dt)
        @error "Rk4: Specify the time step with the MaxStep option."
    end
            
    # Time history, making sure the end time is explicity represented.
    t  = (ts[1]:dt:ts[end]);
    if t[end] != ts[end]
        t[end+1, 1] = ts[end];
    end

    ns = maximum(size(t));            # Number of samples
    nx = length(x0);                  # Number of states
    x  = [x0[:]'; zeros(ns-1, nx)];   # State history
    xk = x0[:];                       # Current state
    hdt = 0.5 * dt;                   # Half of the time step
    
    # If the user provided an OutputFcn, call it.
    if (options.OutputFcn != [])
        options.OutputFcn(ts, x0, "init");
    end
        
    # Propagate from each k to k+1.
    for k = 1:ns-1
        
        # The last sample may be cut short to have the correct end time.
        if k == ns-1
            dt = t[k+1] - t[k];
            hdt = 0.5 * dt;
        end
            
        # From k to k+1
        tk = t[k];                                        # Current time
        k1 = ode(tk,       xk);                           # RK derivatives
        k2 = ode(tk + hdt, xk + hdt * k1);                # ...
        k3 = ode(tk + hdt, xk + hdt * k2);
        k4 = ode(tk +  dt, xk +  dt * k3);
        xk = xk + (k1 + 2 * k2 + 2 * k3 + k4) * (dt / 6); # Updated state
        x[k+1, :] = xk[:];                                # Store.
            
        # If the user provided an OutputFcn, call it.
        if (options.OutputFcn != [])
            if options.OutputFcn(t[k+1], xk, "");
                t = t[1:k+1];
                x = x[1:k+1, :];
                break;
            end
        end
    
    end
    
        # If the user provided an OutputFcn, call it.
    if (options.OutputFcn != [])
        options.OutputFcn([], [], "done");
    end

    return t, x
end

function rk4(ode, ts, x0, dt)
    """
        Wrapper function for RK4 when only dt is specified, rather than a full ODE_SET
            t, x = rk4(ode, ts, x0, dt);
    """
    options = ODE_SET()
    options.MaxStep = dt

    return rk4(ode, ts, x0, options)
end
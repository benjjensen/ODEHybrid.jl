module ODEHybrid

using Interpolations

export ODE_SET
"""
Copy of MATLAB's ODE SET struct
"""
Base.@kwdef mutable struct ODE_SET
    AbsTol = 1e-6
    BDF = []
    Events = []
    InitialStep = []
    Jacobian = []
    JConstant = []
    JPattern = []
    Mass = []
    MassSingular = []
    MaxOrder = []
    MaxStep = []
    NonNegative = []
    NormControl = []
    OutputFcn = []
    OutputSel = []
    Refine = []
    RelTol = 1e-3 
    Stats = []
    Vectorized = []  
    MStateDependence = []
    MvPattern = []
    InitialSlope = []
end

include("interpd.jl");
include("rk4.jl");
include("rkadapt.jl");
include("rkfixed.jl");
include("TimeSeriesLogger.jl")
include("state_to_vector.jl")
include("vector_to_state.jl");
include("vectors_to_states.jl")


export odehybrid
""" 
Hybrid continuous and discrete propagation.

This function propagates an ordinary differential equation along with a
discrete update equation in a manner similar to MATLAB's ode45 function
(which only propagates an ordinary differential equation). This is useful
for implementing discrete-time controllers or simulating processes that
are updated at discrete intervals.

A large number of examples can be found in examples_odehybrid or by
entering the following at the command line:

# home = fileparts(which('examples_odehybrid'));
# web(fullfile(home, 'html', 'examples_odehybrid.html'));

Interfaces:

    t, yc, td, yd = odehybrid(solver, ode, de, 
                                dt, ts, yc0, yd0);

    t, yc1m, td, yd1n = odehybrid(solver, ode, de,
                                dt, ts, [yc1, yc2m] [yd1, yd2n]);

    t, ..., td, ..., te, yc1m, yd1n, ie = odehybrid(solver, ode, de, 
                                                      dt, ts, yc0, yd0);

    sol = odehybrid(solver, ode, [de1, de2, ...], [dt1, dt2, ...], 
                        ts, yc0, yd0);

    sol = odehybrid(..., [options], [log]);
    sol = odehybrid(...);

Inputs:

    solver  Continuous-time propagator to use, e.g. ODEHybrid.rkadapt
    ode     Ordinary differential equation to use with solver. The interface
              should be fun(t, xc1, xc2, ..., xcm, xd1, xd2, ..., xdn) where
              xc1, xc2, ..., xcm are the continuous states and xd1, xd2, ..., 
              xdn are the discrete states. It should return the derivatives of
              the continuous states (m outputs).
    de      Discrete update equation(s) (either a function or cell
              array of functions) with the same inputs as ode but
              outputing the updated continuous and discrete states (n+m
              outputs).
    dt      Time step(s) of discrete update equation(s). If de is an
              array of multiple functions, this should be an array of the same
              size.
    ts      Time span, [t_start, t_end]
    yc0     Array of initial continuous states
    yd0     Array of initial discrete states
    options (Optional) options structure from odeset
    log     (Optional) TimeSeriesLogger for logging in the ode and de. If a
              log is passed in, both ode and de *must* be able to accomodate
              having a log input or not as the final argument. E.g., |ode|
              will be called as: ode(t, xc1, ..., xd1, ..., xdn) and
              ode(t, xc1, ..., xd1, ..., xdn, log).

Outputs:

    t       Times corresponding to continuous states (nc-by-1)
    yc1m    Continuous state outputs (m outputs, each nc-by-1)
    td      Times corresponding to discrete state updates (nd-by-1)
    yd1n    Discrete state outputs (n outputs, each nd-by-1)
    te      Times corresponding to events (ne-by-1)
    yce1m   Continuous states at events (m outputs, each ne-by-1)
    yde1n   Discrete states at events (n outputs, each ne-by-1)
    ie      Index of triggering event. See documentation in odeset for more 
              on events.
    sol     If a single output is requested, it will be a structure with
              fields for each of the individual outputs. The various
              continuous states will be grouped in 'yc', the discrete into
              'yd', the continuous states at events into 'yce', and the
              discrete states at events into 'yde'.

    (NOTE that events are not currently supported!)

Example:

    # This is a quick example of simulating an unstable continuous system with
    #   a stabilizing discrete-time controller.

    ode = (t, x, u) -> [0 1; 2 0] * x + [0; 1] * u;  # Differential equation
    de  = (t, x, u) -> (x, -[8 4] * x);              # Discrete update equation
    dt  = 0.1;                                       # Discrete eq. time step (float or float array)
    ts  = [0 5];                                     # From 0 to 5s
    x0  = [1.0; 0.0];                                # Initial continuous state (floats, not Ints)
    u0  = 0;                                         # Initial discrete state
    t, x, tu, u = odehybrid(ODEHybrid.rkadapt, ode, de, dt, ts, x0, u0);    # Simulate!
    
    plot(t, x, xlabel = "Time", title = "Example", label = ["x₁" "x₂"])
    scatter!(tu, u, label = "u")

See also: examples_odehybrid, ode45, rkadapt, rkfixed.

Online doc: http://www.anuncommonlab.com/doc/odehybrid/odehybrid.html

Copyright 2014 An Uncommon Lab
"""
function odehybrid(solver, ode, de, dt, ts, yc0, yd0, options = ODE_SET(), log = [])

    # Check inputs
    if maximum(size(ts)) != 2
        @error "ODEHybrid: Invalid Arguments! Propagation window ts should specify BOTH start and stop times and ONLY those"
    end

    if (options == [])     # Meaning a log was passed in but no options specified
        options = ODE_SET()
    end
   
    # See if we're passing in separated states (the "full" version).
    if (!isa(yc0, Array{<:Number})) || (!isa(yd0, Array{<:Number}))
        # Get the full inputs
        output = odehybridfull(solver, ode, de, dt, ts, yc0, yd0, options, log) 

    else 
        # Otherwise, just pass everything on directly to the odehybridcore
        output = odehybridcore(solver, ode, de, dt, ts, yc0, yd0, options, log) 
    end

    return output
end

export odeplot   
""" 
Attempted port of MATLAB's odeplot that *appears* to match outputs.
Saves data and generates a plot.
"""  
# TODO currently uses globals to store, which is not ideal...
function odeplot(t, y, flag)

    error = false

    try
        if flag == "init"
            global odeplot_x_vals = [t[1]]
            global odeplot_y_vals = y

        elseif flag == "done"
            a = scatter(odeplot_x_vals, odeplot_y_vals', title = "ODE Plot")
            a = plot!(odeplot_x_vals, odeplot_y_vals', label = false, title = false)
            display(a)

        else
            global odeplot_x_vals = [odeplot_x_vals[:]; t]
            global odeplot_y_vals = [odeplot_y_vals y]
        end
    catch 
        error = false 
    end

    return error
end

export vec_to_mat  
"""
Helper function that converts a N-element vector with each element of length M into an N x M matrix 
"""
function vec_to_mat(vec)

    N = length(vec)
    M = maximum(size(vec[1]))
    mat = zeros(N, M)
    for i = 1:N
        mat[i, :] = vec[i]
    end

    return mat 
end;

# Run the continuous-discrete-input version of odehybrid.
function odehybridfull(solver, ode, de, dt, ts, yc0, yd0, options = ODE_SET(), log = [])
    """
        Formats arguments as needed and then passes them into ODE Hybrid Core. 
        Also unconverts outputs before returning.
    """

    # Let the user pass in arrays or anything else, but always make sure we
    # work with arrays. This simplifies life when we have to dynamically
    # pass these into functions.

    # UPDATED: Only make non-structs into an array
    if !isa(yc0, Array) && (isa(yc0, Number) || isa(yc0, Char))
        yc0 = [yc0]
    end
    if !isa(yd0, Array) && (isa(yc0, Number) || isa(yc0, Char))
        yd0 = [yd0]
    end
    if !isa(de, Array) 
        de = [de]
    end
    
    # Create the initial continuous state vector.
    yc0v = state_to_vector(yc0);

    # Create a function to expand the vector into the state, pass the state
    #   to the ODE, and turn the result back into a vector.
    if (log == [])
        ode_v = (_t, _yc, _yd, _varargin = [])  -> run_ode(ode, _t, _yc, _yd, yc0, _varargin...)
    else 
        ode_v = (_t, _yc, _yd, _varargin = log) -> run_ode(ode, _t, _yc, _yd, yc0, _varargin)
    end


    if isa(dt, Array)
        de_v = [];
        for k = 1:maximum(size(dt))
            if (log == [])
                temp = (_t, _yc, _yd, _varargin = []) -> run_de(de[k], _t, _yc, _yd, yc0, _varargin...) 
                push!(de_v, temp)
            else
                temp = (_t, _yc, _yd, _varargin = log) -> run_de(de[k], _t, _yc, _yd, yc0, _varargin) 
                push!(de_v, temp)
            end
        end
    else
        if (log == [])
            de_v = (_t, _yc, _yd, _varargin = []) -> run_de(de[1], _t, _yc, _yd, yc0, _varargin...) 
        else 
            de_v = (_t, _yc, _yd, _varargin = log) -> run_de(de[1], _t, _yc, _yd, yc0, _varargin)
        end 
    end

    # Propagate the vector versions of the ODE and DE.
    outputs = odehybridcore(solver, ode_v, de_v, dt, ts, yc0v, yd0, options, log);


    ## IN MATLAB, DIFFERENT RETURN OPTIONS ARE AVAILABLE BUT WE RETURN THEM ALL HERE
    return outputs 
end

function odehybridcore(solver, ode, de, dt, ts, yc0, yd0, options = ODE_SET(), log = []) 
    """ 
        Inputs: 

            solver  Continuous-time propagator to use, e.g. @ode45
            ode     Ordinary differential equation to use with solver. The interface
                      should be fun(t, xc1, xc2, ..., xcm, xd1, xd2, ..., xdn) where
                      xc1, xc2, ..., xcm are the continuous states and xd1, xd2, ..., 
                      xdn are the discrete states. It should return the derivatives of
                      the continuous states (m outputs).
            de      Discrete update equation(s) (either a function handle or cell
                      array of function handles) with the same inputs as ode but
                      outputing the updated continuous and discrete states (n+m
                      outputs).
            dt      Time step(s) of discrete update equation(s). If de is a cell
                      array of multiple functions, this should be an array of the same
                      size.
            ts      Time span, [t_start, t_end]
            yc0     Array of initial continuous states
            yd0     Array of initial discrete states
            options (Optional) options structure from odeset
            log     (Optional) TimeSeriesLogger for logging in the ode and de. If a
                      log is passed in, both ode and de *must* be able to accomodate
                      having a log input or not as the final argument. E.g., |ode|
                      will be called as: ode(t, xc1, ..., xd1, ..., xdn) and
                      ode(t, xc1, ..., xd1, ..., xdn, log).

        Outputs:

            t       Array of time values for continuous function 
            yc      Array of continuous states corresponding to times in t
            td      Array of time values for discrete function
            yd      Array of discrete states corresponding to times in td
            te      Times for each event
            yce     Array of continuous states corresponding to events in te
            yde     Array of discrete states corresponding to events in te
            ie 

            (NOTE that events are not currently supported!)

    """

    # For multiple sample rates, we'll expect the discrete updates to be
    # stored in a cell array. If the user is just passing in a single
    # discrete update, it might not be in an array. Put it into one.
    if !isa(de, Array)
        de = [de]
    end

    ts = float.(ts) # Make sure that the times are floats for the "eps" function later
    
    # # Make sure dt is a row vector.
    # if size(dt, 1) > 1
    #     dt = dt';
    # end

    # Check for the things we don't support.
    if options.Vectorized != [] 
        @error "ODEHybrid: 'Options' doesn't support vectorized functions!"
    end
    if options.OutputSel != []
        @error "ODEHybrid: 'Options' doesn't support the OutputSel option because states need not be vectors!"
    end

    # Figure out the time steps. First, calculate the resolution with which 
    # we will be able to differentiate steps. Then, make arrays for each 
    # separate time step.
    epsilon = 2 * maximum((diff(ts, dims = 2) ./ dt) .* eps.(dt) .+ eps(ts[2]));
    
    if isa(dt, Array)
        tds = Array{Any}(undef, maximum(size(dt)))
        for k = 1:maximum(size(dt))
            tds[k] = Array(ts[1]:dt[k]:ts[end])
        end
    else
        tds = Array(ts[1]:dt:ts[end])
    end

    # Sort the necessary discrete times into a single list and remove
    # doubled steps. We now have a list of all times at which to break for
    # a discrete step of one or more discrete update functions.

    # TODO not very elegant 
    times = []
    for el in tds
        for i in el 
            push!(times, i)
        end
    end

    td = sort(times); 
    doubled_step_indices = ([false; diff(td) .< epsilon] )
    deleteat!(td, doubled_step_indices) # filter!
    # td[doubled_step_indices] .= []
    td = td[:]
    
    # We're going to overwrite the output_fcn, so store the original.
    orig_outputfcn = options.OutputFcn 

    # We're going to overwrite the event_fcn, so store the original.
    orig_eventfcn = options.Events 
        
    # Set the maximum time step to be the smaller of what's already specified or dt.
    if isempty(options.MaxStep)
        options.MaxStep = minimum(dt)
    else
        options.MaxStep = min(options.MaxStep, minimum(dt))
    end
    
    # Set the initial step (Check if it is uninitialized)
    if isempty(options.InitialStep)
        options.InitialStep =  min(0.1 * minimum(dt), 0.01 * diff(ts, dims = 2)[1])
    end

    # Initialize outputs.
    t   = ts[1]; 
    yc  = yc0';
    te  = [];
    yce = []; 
    ie  = []; # Events

    if isa(yd0, Array)
        yd = zeros(maximum(size(td)), length(yd0))
        yd[1, :] = yd0[:]
        yde = []
    elseif isa(yd0, Tuple)
        yd = zeros(maximum(size(td)), length(yd0))
        yd[1, :] .= yd0
        yde = []
    else # isa Struct
        yd = Array{Any}(undef, maximum(size(td))) 
        yd[1] = yd0
        yde = [];
    end
    
    # Call the output function if there is one.
    if (orig_outputfcn != [])
        orig_outputfcn(ts, yc0, yd0, "init")
    end
    
    # Loop through all k, upda)ting from k to k+1.
    if isa(dt, Array)
        counts = zeros(size(dt))
    else
        counts = [0.0]
    end
    
    for k = 1:maximum(size(td))
        
        ####################
        # Discrete Updates # 
        ####################
        
        if k <= maximum(size(td))

            # Find out what should discrete functions should be called at this time.
            to_tick = findall(abs.(counts .* dt .+ ts[1] .- td[k]) .≤ epsilon)
            
            # Call the relevant discrete updates.
            for z = to_tick
                yc0, yd0 = de[z](td[k], yc0[:], yd0);
            end

            # We can now store the updated state for k. 
            if isa(yd0, Array)
                yd[k, :] = yd0[:];
            elseif isa(yd0, Tuple)
                yd[k, :] .= yd0;
            else # isa Struct
                yd[k] = yd0;
            end

            # Tick them.
            counts[to_tick] = counts[to_tick] .+ 1
    
            # Call the output function with the udpated discrete states.
            if (orig_outputfcn != []) && orig_outputfcn(td[k], yc0, yd0, "");
                td = td[1:k];
                if isa(yd, Array)
                    yd = yd[1:k];
                else
                    yd = yd[1:k, :];
                end
                break;
            end
            
        end
        
        ######################
        # Continuous Updates #
        ######################
        
        # If we're at the end of the discrete steps, see if there's still a
        # little more left to simulation in continuous time. If so, create
        # the correct time span. If there's no continuous time left, just
        # break (we're done). And if we're not at the end of the discrete
        # steps, just create the right time window.
        if k == maximum(size(td))
            if ts[2] - td[k] > epsilon
                tsk = [td[k], ts[2]];
            else
                break;
            end
        else
            tsk = [td[k], td[k+1]];
        end
        
        # Make the output function with the new discrete states.
        if (log != [])  || (options.OutputFcn != [])
            options.OutputFcn = (t, y, flag) -> log_frame(t, y, flag, ode, yd0, log, orig_outputfcn)
        end

        if (options.Events != [])
            options.Events = (t, y) -> orig_eventfcn(t, y, yd0);
        end
        
        # If there are events, output one way. If not, output the basic way.
        if (options.Events != [])
            @warn "WARNING! Events are not currently allowed using the ODE Hybrid RK solvers!"
            tk, yck, tek, yek, iek = solver( (_t, _y) -> ode(_t, _y, yd0), tsk, yc0, options)
            te  = [te; tek]
            yce = [yce; yek]
            ie  = [ie; iek]
            for z = 1:maximum(size(tek))
                if isa(yd0, Array)
                    yde[end+1, 1:length(yd0)] = yd0[:]
                else
                    yde = [yde; yd0[:]]
                end
            end

        else 
            if isa(yd0, Array) && length(yd0) == 1
                tk, yck = solver( (t, y) -> ode(t, y, yd0[1]), tsk, yc0, options)
            else
                tk, yck = solver( (t, y) -> ode(t, y, yd0), tsk, yc0, options)
            end
        end 
        
        # Add the outputs to the list.
        t  = [t;  tk];
        yc = [yc; yck];
        
        # See if the function terminated early.
        if abs(tk[end] - tsk[end]) > 2*eps(tsk[end])
            td = td[1:k];
            if isa(yd, Array)
                yd = yd[1:k];
            else
                yd = yd[1:k, :];
            end
            break;
        end
        
        # Set the initial step for next time to be the largest step we took
        #    during this iteration.
        options.InitialStep = maximum(diff(tk));
        
        # Get ready for the next step.
        yc0 = yc[end, :];        
    end
    
    # Call the output function if there is one.
    if (orig_outputfcn != [])
        orig_outputfcn(ts, yc0, yd0, "done");
    end


    # Reset because these are modified
    options.OutputFcn = orig_outputfcn 
    options.Events  = orig_eventfcn

    return [t, yc, td, yd, te, yce, yde, ie]
end

function output_fcn(f, t, ycv, yd, flag, yc0)
    """
        Custome output function that ignores "init" and "done" flags.

        Inputs:

            t:    Time                                             |  Scalar
            ycv:  Values of continuous-time ode at given time      |  column vector
            yd:   Values of discrete-time system at given time     |  Cell of Array {[...]}
            flag: String dictating which task to run               |  String
            yc0:  Initial values/format for continuous-time eqn    |  Cell 


        Outputs:

            status:                                                |  Scalar
   
    """

    # We're going to rely on arrays to pass these to the original
    #    output function, so make sure they're arrays.
    if !isa(yc0, Array)
        yc0 = [yc0]
    end
    if !isa(yd, Array)
        yd = [yd]
    end
        
    # For 'done', all of the states should be empty.
    if flag == "done"

        # Pull out the continuous state first (this is how we store it).
        yc = vector_to_state(ycv, yc0);

        for k = 1:maximum(size(yc))
            yc[k] = [];
        end
        for k = 1:maximum(size(yd))
            yd[k] = [];
        end
        status = f([], yc[:], yd[:], flag);

    # For init, we pass the time span along with the initial state.
    elseif flag == "init"

        # Pull out the continuous state first.
        yc = vector_to_state(ycv, yc0);
        status = f(t, yc[:], yd[:], flag);

    # Otherwise, pass along the states and the flag.
    else

        for k = 1:maximum(size(t))
            # Pull out the continuous state first.
            yc = vector_to_state(ycv[:, k], yc0);
            status = f(t[k], yc[:], yd[:], flag);
        end
    end
        
    return status
end

function event_fcn(f, t, ycv, yd, yc0)
    """
        Custome output function that ignores "init" and "done" flags.
        NOTE that events are not currently supported

        Inputs:

            f:    Function to be called at each event              |  Function
            t:    Time                                             |  Scalar
            ycv:  Values of continuous-time ode at given time      |  column vector
            yd:   Values of discrete-time system at given time     |  Cell of Array {[...]}
            yc0:  Initial values/format for continuous-time eqn    |  Cell 


        Outputs:

            status: Triggers a stop when function fails            |  Boolean

    """

    # We're going to rely on cell arrays to pass these to the original
    # output function, so make sure they're cell arrays.
    if !isa(yc0, Array)
        yc0 = [yc0]
    end
    if !isa(yd, Array)
        yd = [yd]
    end
    
    # Pull out the continuous state first.
    yc = vector_to_state(ycv, yc0);
    return f(t, yc[:], yd[:]);
end

function log_frame(t, y, flag, ode, yd, log = [], f = [])
    """
        Inputs:

            t:    Time                                                  |  Scalar
            y:    Continuous Value                                      |  
            flag: String dictating what to run                          |  String ∈ {"init", "done", ""}
            ode:  Ordinary differential equation                        |  Function 
            yd:   Discrete Value at time t                              |  
            log:  (Opt) Time Series Logger used to track and plot data  |  TimeSeriesLogger
            f:    (Opt) User's original OutputFcn                       |  Function 

        Outputs:

            status: Triggers a stop when function fails            |  Boolean
    """


    # We never tell integration to stop (but the user's function can override this).
    status = false;

    # Only call the ODE with the logger on "normal" steps (steps with 
    # neither 'init' nor 'done' flags).
    if isempty(flag)

        # If there's a TimeSeriesLogger, use it.
        if log != []  #!isempty(log)
            N = isa(t, Array) ? maximum(size(t)) : 1;
            for k = 1:N
                ode(t[k], y[:, k], yd, log);
            end
        end


        # Call the user's OutputFcn if provided.
        if (f != [])
            status = Bool(f(t, y, yd, flag));
        end
    end

    return status
end

function run_ode(ode, t, ycv, yd, yc0, varargin = [])
    """
        Formats an ordinary differential equation to be used in ODE Hybrid Core.
            Converts the provided state vector into the appropriate format and 
            then runs the ODE provided by the user, returning a vector output of the ODE.
            Called in ODE Hybrid Full when inputs are not pure numerics


        Inputs: 

            ode:  User-provided ODE function of form (t, x, u) -> ()            |  function
            t:    Current time value                                            |  Scalar
            ycv:  State of continuous system at provided time                   |
            yd:   State of discrete system at provided time                     |
            yc0:  Initial state of continuous system (used to match format)     |
            varargin: (Optional) Extra arguments that may be passed in          |  Array

        Outputs: 
            dydvt: Output of the provided ODE

    """

    # Pull out the continuous state first (this is how we store it).
    yc, _ = vector_to_state(ycv, yc0);


    # Get the state differences.

    # TODO (Clumsy) Method of allowing for tuples, arrays, and structs to all be used 
    yc = isa(yc, Tuple) ? yc : (yc,)
    yd = isa(yd, Tuple) ? yd : (yd, )

    varargin = ((varargin != []) && !isa(varargin, Tuple)) ? (varargin,) : varargin

    state_difference = ode(t, yc..., yd..., varargin...); 
    
    # Convert the derivative to a vector.
    dyvdt = state_to_vector(state_difference);

    return dyvdt
end

function run_de(de, t, ycv, yd, yc0, varargin = [])

    """
        Formats a discrete equation to be used in ODE Hybrid Core.
            Converts the provided state vector into the appropriate format and 
            then runs the discrete equation provided by the user, returning a vector output.
            Called in ODE Hybrid Full when inputs are not pure numerics

        Inputs: 

            de:   User-provided discrete function of form (t, x, u) -> ()        |  function
            t:    Current time value                                             |  Scalar
            ycv:  State of continuous system at provided time                    |
            yd:   State of discrete system at provided time                      |
            yc0:  Initial state of continuous system (used to match format)      |
            varargin: (Optional) Extra arguments that may be passed in           |  Array

        Outputs: 
            ycv:  Value of continuous function that discrete function outputs    |
            yd:   Value of discrete function at output                           |

    """

    # Pull out the continuous state first (this is how we store it).
    yc, _ = vector_to_state(ycv, yc0);

    # TODO (Clumsy) Method of allowing for tuples, arrays, and structs to all be used 
    yc = isa(yc, Tuple) ? yc : (yc,)
    yd = isa(yd, Tuple) ? yd : (yd, )

    # If varargin isn't empty AND whatever its value is isnt a tuple, make it a tuple so we can splat 
    varargin = ((varargin != []) && !isa(varargin, Tuple)) ? (varargin,) : varargin

    # Adding in varargin splat... (hopefully log when it isn't [])
    yc, yd = de(t, yc..., yd..., varargin...);
    
    # Convert to vectors.
    ycv = state_to_vector(yc);
    
    return ycv, yd
end

function output_discrete_states(states, yd0)

    """
        Turns an [n × m] array of discrete states into an [m × 1] array of
            discrete state lists each with n entries.
    """
    
    discrete_states = []
    for k = 1:length(yd0)
        push!(discrete_states, states[:, k])
    end
    return discrete_states
end

end # module

""" TODO:
        - Figure out "unique figure" 
        - Figure out how to safely overload "plot" for a TSL
"""

export TimeSeriesLogger
Base.@kwdef mutable struct TimeSeriesLogger

    """ TimeSeriesLogger struct
        A struct for keeping track of numerous pieces of data over time throughout
        a process. It's designed to be easy to use and relatively quickly. It
        logs numeric types of any size (scalars, vectors, or matrices), as long
        as the size is consistent from sample to sample.
         
        Example:
         
            log = TimeSeriesLogger();    % Make a new logger.
            x   = [0; 0];
            for t = 1:100                % Make a random walk.
                x = x + randn(2, 1);     % Take a single step.
                add!(log, 'walk', t, x); % Log a single sample, x, at time t.
            end

            tsl_plot(log)                % Plot everything.
         

        We can also access specific logs by their names. In the above, we only
            have one log ('walk'). Let's get that log from the logger.
            
            t_log, x_log = get_log(log, "walk");    % Get a specific log.
            plot(x_log[:, 1], x_log[:, 2]);         % Do something with it.
        
        We can make sure a log exists before trying to do anything with it:
         
            if contains(log,"walk")
                x_log = get_log(log, "walk");
                plot(x_log[:, 1], x_log[:, 2]);
            end
        
        If we want to log something but don't want it plotted when we call
        |plot|, then we can pass in false to the logger when we add the signal.
         
            add!(log, "var1", t, x, false);
         
        To group items into different figures when they're plotted, give them
            "groups":
         
            add!(log, "var1", t, x, true, "group1");
            add!(log, "var2", t, y, true, "group2");
            add!(log, "foo",  t, bar, true, "group1");
         
        All of the items with common group names will be added to the same
            figures as subplots.
        
        Finally, we can clear out the logger with |initialize|. This deletes all
            data and returns the logger to its initial state.
        
            initialize(log);
        
        For more information on the various methods, see the individual methods
        with |doc TimeSeriesLogger|:
         
        
        Methods
        -------
        
            initialize! Clear out all data.
            add!        Add a single data point to a time series.
            contains    See if a time series exists (by name).
            get_log     Return a log by name.
            tsl_plot    Plot the logged series.
        
        Notes
        -----
        
            TimeSeriesLogger was created as the logging mechanism for odehybrid.
            However, it's not dependent on that project and so can be used anywhere.
            
            TimeSeriesLogger is not related to the 'timeseries' class in MATLAB.
            
            Since this class never knows how much more data is going to be logged, it
            can't preallocate the appropriate amount of space. However, increasing
            the size of its data stores by one row every time a new sample is added
            is very slow. To combat this, the data store starts off small and doubles
            whenever the store is at capacity. Say we're logging a 4x1 signal. When
            the first sample is added (say it's log k), data{k}{2} will be 1x4. When
            the second signal is added, it becomes 2x4. For the third, it's 4x4, then
            8x4, 16x4, etc. A separate counter stores how much of this allocated 
            space is currently used. This reduces the number of allocations from n to
            log2(n). Practically, it saves a little time during logging without too
            much complexity.
        
        Online doc: http://www.anuncommonlab.com/doc/odehybrid/TimeSeriesLogger.html
        
        Copyright 2014 An Uncommon Lab
    """

    # Properties
    names = []   # Unique identifier (within the group)
    times = []   # Time stamp corresponding to each data point
    data  = []   # Data, stored as {t, x}
    show  = []   # True iff the data should be plotted
    sizes = []   # Amount of allocated space
    count = []   # Amount of allocated space currently used
    group = []   # Plot group identifier

end;

export add!
function add!(tsl::TimeSeriesLogger, name, t, x, show_it = true, group = "")
    """
        Add a new log or append to an existing log by name.
     
            add!(log, "var1", t, x);
            add!(log, "var1", t, x, true);  # Show log in plots (default)
            add!(log, "var1", t, x, false); # Don't show log in plots
        
        Signals can be grouped together into figures by given them
            a common group argument. Here, both var1 and var2 logs will
            be plotted together.

            add!(log, "var1", t,  x, true, "group1");
            add!(log, "var2", t2, y, true, "group1");
        
        Returns true iff a new log was created.

        NOTE that unlike in MATLAB, the signals are stored separately in times (t)
            and data (x), (e.g., the kth pair would be log.times[k] and log.data[k])
    """

    # See if we've already added this variable 
    index = isempty(tsl.names) ? [] : (tsl.names .== name)
    log_is_empty = (index == []) || !maximum(index)


    # NEW
    if !isa(x, Array)
        x = [x]
    end

    if log_is_empty 
        # Initialize everything first time through

        new_times = Array{Float64, 1}()
        push!(new_times, t)
        push!(tsl.times, new_times)

        new_data = Array{Any, 1}()
        push!(new_data, x[:])
        push!(tsl.data, new_data)
        push!(tsl.names, name)

        push!(tsl.show, show_it)

        push!(tsl.sizes, 1)
        push!(tsl.count, 1)
        push!(tsl.group, group)
    
    else # otherwise, append data to existing set 

        # If we have used up all the reserved space, double it 
        c = tsl.count[index][1] + 1

        if c > tsl.sizes[index][1]
            push!(tsl.times[index][1], t)
            for i = 1:(tsl.sizes[index][1] - 1)
                push!(tsl.times[index][1], 0.0)
            end

            push!(tsl.data[index][1], x[:])
            for i = 1:(tsl.sizes[index][1] - 1)
                push!(tsl.data[index][1], zeros(maximum(size(x))))
            end

            tsl.sizes[index]   .= 2 * tsl.sizes[index][1]
        else # If space remains, add to it 
            tsl.times[index][1][c]   = t 
            tsl.data[index][1][c, :][1] .= x[:]
        end 

        tsl.count[index] .= c

        # We should show this data if the user has ever asked to plot it 
        tsl.show[index]  .= show_it || any(tsl.show[index])
        tsl.group[index] .= group

    end   

    return log_is_empty
end;
     
export contains
function contains(tsl::TimeSeriesLogger, name) 
    # Return 'true' iff the TimeSeriesLogger contains this name 

    indices = findall(tsl.names .== name)
    return !isempty(indices)
end;

export get_log
function get_log(tsl::TimeSeriesLogger, name)
    """ 
        Return a specific log by name.
        
        If one output is request, it returns the data. If two are 
        requested, it returns the time and the data. Returns empty if
        the logs don't exist (use the 'contains' function to test for
        this).

        Example:
        
            x = get_log(log, "signal1");
            t, x = get_log(log, "signal1");
        
        NOTE that |t| will be ns-by-1 and x will be ns-by-nx, where
        nx is the number of elements in a single sample.
    """

    index = (tsl.names .== name)

    if isempty(index)
        return ([], [])
    else 
        t = tsl.times[index][1][1:tsl.count[index][1]]
        x = vec_to_mat(tsl.data[index][1][1:tsl.count[index][1]])
    end

    return (t, x)
end;

export tsl_plot
function tsl_plot(tsl::TimeSeriesLogger, x_label) 
    # Plot all of the signals, grouped appropriately into figures.
    #   A custom x label can be added to figures as well

    groups = unique(tsl.group)

    # For each group 
    for g = 1:maximum(size(groups))

        # See what is to be shown and count them 
        to_show = (tsl.group .== groups[g]) .& tsl.show[:]
        n = sum(to_show) #

        # Bail if there is no need to show anything 
        if n == 0;
            continue;
        end

        plts = Array{Plots.Plot{Plots.GRBackend}, 1}(undef, n)

        ### TODO    ###############
        # Set up the figure  
        # unique_figure("LoggerDefaultPlot$g", [groups[g] " Log"])

        # % Set up the figure.
        # unique_figure(sprintf('LoggerDefaultPlot%d', g), ...
        #               'Name', [groups{g} ' Log']);

        # clf();  
        ###########################

        to_show = Int.(to_show)
        plt_idx = 1
        for k = 1:maximum(size(to_show))
            if (to_show[k] == 0)
                continue
            end
            
            N = tsl.count[k] #to_show[k]]
            plot_times = tsl.times[k][1:N]
            plot_data  = vec_to_mat(tsl.data[k][1:N, :])
            plot(plot_times, plot_data)

            ylabel!(tsl.names[k]) #to_show[k]])
            plts[plt_idx] = xlabel!(x_label)
            plt_idx += 1
        end
        display(plot(plts..., layout = (n, 1)))
    end
end;

function tsl_plot(tsl::TimeSeriesLogger) 
    # Wrapper function that plots a TimeSeriesLogger without a label 
    return tsl_plot(tsl, " ")
end;

export initialize!
function initialize!(tsl::TimeSeriesLogger)  
    # Clear everything out and start fresh.

    tsl.names = [];
    tsl.times = []; 
    tsl.data  = [];
    tsl.show  = [];
    tsl.sizes = [];
    tsl.count = [];
    tsl.group = [];

    return Nothing
end;

# NOT CURRENTLY FINISHED
function unique_figure(id, name = " ", args = [])
    """
        h = unique_figure(id, varargin)
        
        Creates a figure with the given ID (tag) or selects it if it already
        exists, allowing one to easily reuse the same figure window identified
        with text instead of handles. This is useful when, e.g., running a script
        many times after clearing between runs without having to hard-code figure
        numbers (which can become hard to keep track of).
        
        Any additional arguments are passed along to the figure's |set| method.
        
        Example:
        
        h = unique_figure('trajectory');
        
        Example with extra arguments:
        
        h = unique_figure('trajectory', 'Name', 'Trajectory'); 
    """
    # See if there is a figure with this ID already 
    # h = findobj('Type', 'figure', 'Tag', id);
    # if isempty(h)
    h = plot("Tag", id, args);

    # else
            
    #     % If we found it, select it.
    #     figure(h);
        
    #     % If there were any additional arguments, pass them along. Note
    #     % that if there weren't and we called set(h), then it would print h
    #     % to the command window -- probably not what we want, hence the
    #     % condition.
    #     if nargin > 1
    #         set(h, varargin{:});
    #     end
    # end

    return h
end;


function [t, yc, td, yd, te, yce, yde, ie] = odehybridcore(solver, ...
                                                           ode, ...
                                                           de, dt, ...
                                                           ts, yc0, yd0,...
                                                           varargin)
                                                       
    % Ignore growth in loops. We're just going to have this problem.
    %#ok<*AGROW>
    
    % For multiple sample rates, we'll expect the discrete updates to be
    % stored in a cell array. If the user is just passing in a single
    % discrete update, it might not be in a cell array. Put it into one.
    if ~iscell(de)
        de = {de};
    end
    
    % Make sure dt is a row vector.
    if size(dt, 1) > 1
        dt = dt.';
    end
    
    % Set the solver's options if the user provided them.
    if nargin >= 8 && ~isempty(varargin{1})
        options = varargin{1};
    else
        options = odeset();
    end

    % Check for the things we don't support.
    if ~isempty(options.Vectorized)
        error('odehybrid:options', ...
              'odehybrid doesn''t support vectorized functions.');
    elseif ~isempty(options.OutputSel)
        error('odehybrid:options', ...
              ['odehybrid doesn''t support the OutputSel option ' ...
               'because states need not be vectors.']);
    end
    
    % Figure out the time steps. First, calculate the resolution with which 
    % we will be able to differentiate steps. Then, make arrays for each 
    % separate time step.
    epsilon = 2 * max((diff(ts) ./ dt) .* eps(dt) + eps(ts(2)));
    tds = cell(1, length(dt));
    for k = 1:length(dt)
        tds{k} = ts(1):dt(k):ts(end);
    end
    
    % Sort the necessary discrete times into a single list and remove
    % doubled steps. We now have a list of all times at which to break for
    % a discrete step of one or more discrete update functions.
    td = sort([tds{:}]);
    td([false diff(td) < epsilon]) = [];
    td = td(:);
    
    % We're going to overwrite the output_fcn, so store the original.
    orig_outputfcn = options.OutputFcn;

    % We're going to overwrite the output_fcn, so store the original.
    orig_eventfcn = options.Events;
    
    % If the user passed in a log, use it.
    add_logging = nargin >= 9 && ~isempty(varargin{2});
    if add_logging
        
        % Add logging to the output (preserving the user's outputfcn if
        % it's provided).
        log = varargin{2};
        for k = 1:length(dt)
            de{k} = @(t, yc, yd) de{k}(t, yc, yd, log);
        end
        
    else
        log = [];
    end
    
    % Set the maximum time step to be the smaller of what's already
    % specified or dt.
    if isempty(options.MaxStep)
        options.MaxStep = min(dt);
    else
        options.MaxStep = min(options.MaxStep, min(dt));
    end
    
    % Set the initial step.
    if isempty(options.InitialStep)
        options.InitialStep = min(0.1 * min(dt), 0.01 * diff(ts));
    end
    
    % Initialize outputs.
    t  = ts(1); yc = yc0.';
    te = []; yce = []; ie = []; % Events
    if iscell(yd0)
        yd = cell(length(td), numel(yd0));
        [yd{1, :}] = yd0{:};
        yde = {};
    else
        yd = [yd0.'; zeros(length(td)-1, length(yd0))];
        yde = [];
    end
    
    % Call the output function if there is one.
    if ~isempty(orig_outputfcn)
        orig_outputfcn(ts, yc0, yd0, 'init');
    end
    
    % Loop through all k, updating from k to k+1.
    counts = zeros(size(dt));
    for k = 1:length(td)
        
        %%%%%%%%%%%%%%%%%%%%
        % Discrete Updates %
        %%%%%%%%%%%%%%%%%%%%
        
        if k <= length(td)

            % Find out what should discrete functions should be called at 
            % this time.
            to_tick = find(abs(counts .* dt + ts(1) - td(k)) <= epsilon);

            % Call the relevant discrete updates.
            for z = to_tick
                [yc0, yd0] = de{z}(td(k), yc0(:), yd0);
            end

            % We can now store the updated state for k. 
            if iscell(yd0)
                [yd{k, :}] = yd0{:};
            else
                yd(k, :) = yd0.';
            end

            % Tick them.
            counts(to_tick) = counts(to_tick) + 1;
    
            % Call the output function with the udpated discrete states.
            if    ~isempty(orig_outputfcn) ...
               && orig_outputfcn(td(k), yc0, yd0, '');
                td = td(1:k);
                if iscell(yd)
                    yd = yd(1:k);
                else
                    yd = yd(1:k, :);
                end
                break;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Continuous Updates %
        %%%%%%%%%%%%%%%%%%%%%%
        
        % If we're at the end of the discrete steps, see if there's still a
        % little more left to simulation in continuous time. If so, create
        % the correct time span. If there's no continuous time left, just
        % break (we're done). And if we're not at the end of the discrete
        % steps, just create the right time window.
        if k == length(td)
            if ts(2) - td(k) > epsilon
                tsk = [td(k), ts(2)];
            else
                break;
            end
        else
            tsk = [td(k), td(k+1)];
        end
        
        % Make the output function with the new discrete states.
        if add_logging || ~isempty(options.OutputFcn)
            options.OutputFcn = @(t, y, flag) log_frame(t, y, flag, ...
                                                        ode, yd0, log, ...
                                                        orig_outputfcn);
        end
        if ~isempty(options.Events)
            options.Events = @(t, y) orig_eventfcn(t, y, yd0);
        end
        
        % If there are events, output one way. If not, output the basic
        % way.
        if isempty(options.Events)
            [tk, yck] = solver(@(t, y) ode(t, y, yd0), ...
                               tsk, ...
                               yc0, ...
                               options);
        else
            [tk, yck, tek, yek, iek] = solver(@(t, y) ode(t, y, yd0), ...
                                              tsk, ...
                                              yc0, ...
                                              options);
            if nargout == 1 || nargout >= 3, te  = [te; tek];  end;
            if nargout == 1 || nargout >= 4, yce = [yce; yek]; end;
            if nargout == 1 || nargout >= 6, ie  = [ie; iek];  end;
            if nargout == 1 || nargout >= 5
                for z = 1:length(tek)
                    if iscell(yd0)
                        [yde{end+1, 1:numel(yd0)}] = yd0{:};
                    else
                        yde = [yde; yd0(:).'];
                    end
                end
            end;
        end
        
        % Add the outputs to the list.
        t  = [t;  tk];
        yc = [yc; yck];
        
        % See if the function terminated early.
        if abs(tk(end) - tsk(end)) > 2*eps(tsk(end))
            td = td(1:k);
            if iscell(yd)
                yd = yd(1:k);
            else
                yd = yd(1:k, :);
            end
            break;
        end
        
        % Set the initial step for next time to be the largest step we took
        % during this iteration.
        options.InitialStep = max(diff(tk));
        
        % Get ready for the next step.
        yc0 = yc(end, :)';
        
    end
    
    % Call the output function if there is one.
    if ~isempty(orig_outputfcn)
        orig_outputfcn(ts, yc0, yd0, 'done');
    end

    % If there's only one output, use the structure format.
    if nargout == 1
        sol = struct('t',   t, ...
                     'yc',  yc, ...
                     'td',  td, ...
                     'yd',  yd, ...
                     'te',  te, ...
                     'yce', yce, ...
                     'yde', yde, ...
                     'ie',  ie);
        t = sol;
    end
    
end

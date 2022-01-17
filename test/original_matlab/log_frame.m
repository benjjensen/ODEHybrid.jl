function status = log_frame(t, y, flag, ode, yd, log, f)

    % We never tell integration to stop (but the user's function can
    % override this).
    status = 0;

    % Only call the ODE with the logger on "normal" steps (steps with 
    % neither 'init' nor 'done' flags).
    if isempty(flag)

        % If there's a TimeSeriesLogger, use it.
        if ~isempty(log)
            for k = 1:length(t)
                ode(t(k), y(:, k), yd, log);    
            end
        end

        % Call the user's OutputFcn if provided.
        if ~isempty(f)
            status = f(t, y, yd, flag);
        end
        
    end

end

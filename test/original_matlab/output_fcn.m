function status = output_fcn(f, t, ycv, yd, flag, yc0)

    % We're going to rely on cell arrays to pass these to the original
    % output function, so make sure they're cell arrays.
    if ~iscell(yc0), yc0 = {yc0}; end;
    if ~iscell(yd),  yd  = {yd};  end;
    
    switch flag
        
        % For 'done', all of the states should be empty.
        case 'done'

            % Pull out the continuous state first (this is how we store it).
            yc = vector_to_state(ycv, yc0);

            for k = 1:length(yc)
                yc{k} = [];
            end
            for k = 1:length(yd)
                yd{k} = [];
            end
            status = f([], yc{:}, yd{:}, flag);

        % For init, we pass the time span along with the initial state.
        case 'init'

            % Pull out the continuous state first.
            yc = vector_to_state(ycv, yc0);
            status = f(t, yc{:}, yd{:}, flag);
    
        % Otherwise, pass along the states and the flag.
        otherwise

            for k = 1:length(t)

                % Pull out the continuous state first.
                yc = vector_to_state(ycv(:, k), yc0);
                status = f(t(k), yc{:}, yd{:}, flag);

            end
            
    end
        
end
function varargout = event_fcn(f, t, ycv, yd, yc0)

    % We're going to rely on cell arrays to pass these to the original
    % output function, so make sure they're cell arrays.
    if ~iscell(yc0), yc0 = {yc0}; end;
    if ~iscell(yd),  yd  = {yd};  end;
    
    % Pull out the continuous state first.
    yc = vector_to_state(ycv, yc0);
    [varargout{1:nargout}] = f(t, yc{:}, yd{:});

end
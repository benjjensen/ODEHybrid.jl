""" TODO 
        - Add in additional interpolation modes equivalent to 'nearest', etc...
"""

function interpd(t, x, ti, mode)
    """         
        An interpolation function for time series which may have multiple values
        at a single timestep (e.g., discrete updates). It works like MATLAB's 
        built-in interp1.
         
        xi = interpd(t, x, ti);
        xi = interpd(t, x, ti, mode);
        
        Inputs:
         
            t     Times  (n-by-1)
            x     States (n-by-m)
            ti    Output times (p-by-1)
            mode  '-' for left-most values, '+' for right-most values (default)
            
            (NOTE that unlike the MATLAB version, no additional arguments can be passed in yet (e.g., 'nearest'))
        

        Outputs:
        
            xi  States corresponding to output times (p-by-m)
         
        Example:
         
            Create some times and states. Note that there are two states at t=3.

                t = [2.76, 2.91, 3,   3,   3.12];
                x = [0.2,  0.3,  0.4, 1.1, 1.2];
         
            Create the desired output times.

                ti = (2.8:0.1:3.1);
        
            Interpolate (t, x) linearly at ti, keeping the right-most values.
        
                xi = interpd(t, x, ti)
        
            Interpolate (t, x) linearly at ti, keeping the left-most values.

                xi = interpd(t, x, ti, '-')
         
        See also: examples_odehybrid.
        
        Online doc: http://www.anuncommonlab.com/doc/odehybrid/interpd.html
        
        Copyright 2014 An Uncommon Lab
    """
    
    # Check that t and ti are columns.
    if (length(size(t)) > 1) && (~any(size(t) == 1))
        @error "INTERPD: Invalid inputs - input t should be a column!"
    else
        t = t[:];
    end

    if (length(size(ti)) > 1) && (~any(size(ti) == 1))
        @error "INTERPD: Invalid inputs - input ti should be a column!"
    else
        ti = ti[:];
    end
    
    # Check that the rows of x match the rows of t.
    if size(x, 1) != length(t) && size(x, 2) == length(t)
        x = x';
    end

    # Set the default mode.
    if isempty(mode)
        mode = '+';
    end

    # Create the output.
    xi = zeros(maximum(size(ti)), size(x, 2));
    
    # Find where there are doubled steps (t(doubled(k)) and t(doubled(k)+1) are the same.
    doubled = findall(t[1:end-1] .== t[2:end])

    # If nothing is doubled, just pass along to interp1.
    if isempty(doubled)
        xi = interp1((t,), x, ti);
        return xi
    end
    
    # We will break the interpolation into separate pieces -- those between the doubled 
    #   (or tripled or otherwise duplicated) time steps. Then we'll interpolate that span. 
    #   We'll keep the edges according to whether we're using - or + mode.

    # Make sure we get the first and last spans.
    if doubled[end] != maximum(size(t))
        doubled = [doubled; maximum(size(t))]
    end
    
    # Start with the initial point.
    if ((mode == "-") || (mode == '-')) && doubled[1] != 1
        xi[1, :] .= interp1((t[1:2],), x[1:2, :], ti[1])
    else 
        # We may overwrite this for +
        xi[1, :] .= x[1, :]
    end
    
    # For each doubled index...
    k_last = 1;
    for k = doubled'
        
        # Get the span from the last index to this one.
        tk  = t[k_last:k];
        xk  = x[k_last:k, :];
        
        # Get the range for the requested outputs.
        out_indices = (ti .≥ t[k_last]) .& (ti .≤ t[k]);
        
        # If we're using - mode, drop the first index; we don't want to
        # overwrite the value from the end of the last span with the value
        # from the beginning of this one.

        if (mode == "-" || mode == '-')
            first_one = findall(out_indices .== 1)
            if length(first_one) > 0
                first_one = first_one[1]
                out_indices[first_one] = false
            end
        end

        
        # If there are any outputs in this span...
        if any(out_indices)
            
            # Select the subset of output times.
            tik = ti[out_indices];

            # If there's only one point, don't interp. Otherwise, use interp1.
            if k - k_last == 0
                xi[out_indices, :] .= xk;
            else
                xi[out_indices, :] .= interp1((tk,), xk, tik);
            end
            
        end
        
        # Get ready to move to the next span.
        k_last = k + 1;
        
    end

    if size(xi, 2) == 1
        xi = xi[:, 1]
    end

    return xi
end;

function interp1(x, v, xq, options = Gridded(Linear()))
    """ Simplified Interp1 (Ported from MATLAB)

        x:  Vector containing sample points (in a tuple, pass in (x,) if needs be)
        v:  Vector containing values corresponding to x 
        xq: Vector containing coordinates of the query points
    """

    squeeze = false

    if length(size(v)) == 2
        N, M = size(v)
    elseif length(size(v)) == 1
        N = size(v, 1)
        M = 1
        squeeze = true
    elseif length(size(v)) == 0
        N, M = 1, 1
    else
        @error "INTERP1: Invalid size for dependent variable v!"
    end

    if !isa(xq, Array)
        xq = [xq]
    end
    P = maximum(size(xq))
    interp_vals = zeros(P, M)
    
    # Allow v to be an [n × m] matrix
    if length(size(v)) > 1    
        for i = 1:M
            interpolation = interpolate(x, v[:, i], options)
            interp_vals[:, i] .= interpolation[xq]
        end
    elseif size(v, 2) == 1
        interpolation = interpolate(x, v[:, 1], options)
        interp_vals .= interpolation[xq]
    else
        interpolation = interpolate(x, v, options)
        interp_vals .= interpolation[xq]
    end

    # Return correct array size (or scalar) 
    if squeeze
        if size(interp_vals, 2) == 1
            interp_vals = interp_vals[:, 1]
        end
        if (size(interp_vals, 1) == 1)
            interp_vals = interp_vals[1, :]
        end
    end

    if maximum(size(interp_vals)) == 1
        interp_vals = interp_vals[1]
    end


    return interp_vals
end;

function interpd(t, x, ti)
    """
        Wrapper function for interpd that defaults "mode" to right-most ("+")
    """
    return interpd(t, x, ti, "+")
end;


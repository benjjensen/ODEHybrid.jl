# [scr/vector_to_state.jl] 

""" 
Build a state from a state vector by recursively moving through the 
elements of a more complex example state (matrix, cell array, or struct),
placing the relevant numerical values from the vector into the 
appropriate places. This is used in odehybrid to convert a state vector
into the original complex setup of states. The reverse function is
state_to_vector.

state = vector_to_state(vector, state);

Inputs:

    vector  A vector containing the numerical values from the object
    state   A numeric, cell array, or struct array type representing the
            sizes and types for the values in vector

Outputs:

    state   A numeric, cell array, or struct array type
    count   (Internal use only)

Example:

    mutable struct TestStruc
        a
        bcd
    end
    Base.:(==)(a::TestStruc, b::TestStruc) = (Base.:(==)(a.a, b.a) && Base.:(==)(a.bcd, b.bcd)) # Custom

    x = [ [1 3; 4 2], TestStruc([5; 6], [7:9; 10]),  3.14 ]
    v = state_to_vector(x)
    x2 = vector_to_state(v, x)
    x .== x2[1]

See also: state_to_vector.

Online doc: http://www.anuncommonlab.com/doc/odehybrid/vector_to_state.html

Copyright 2014 An Uncommon Lab
"""
function vector_to_state(vector, state)

    count = 0

    return vector_to_state(vector, state, count)
end


function vector_to_state(vector, state, count)
    """
        Internal function that converts vector to a struct state
    """

    state_copy = deepcopy(state) # So we can setproperty! without modifying original

	# Should be a struct - convert and store in each field.       
    fields = fieldnames(typeof(state_copy)) 
    for k = 1:length(fields)
        temp, count = vector_to_state(vector, getproperty(state_copy, fields[k]), count)
        setproperty!(state_copy, fields[k], temp)
    end

    return state_copy, count
end

function vector_to_state(vector, state::Tuple, count)
    """
        Internal function that converts vector to a tuple state
    """

    result = zeros(length(state))
    # First convert to an array, and then to a tuple
    for k = 1:length(vector)
        result[k], temp = vector_to_state(vector[k], state[k])
    end

    count += length(state)
    return tuple(result...), count  
end

function vector_to_state(vector, state::Char, count)
    """
        Internal function that converts vector to a Char state
    """

    state = Char(vector[count .+ 1]) 
    count = count + 1

    return state, count
end

function vector_to_state(vector, state::Array{Char}, count)
    """
        Internal function that converts vector to a Char array state
    """

    state[:] = vector[count .+ range(1, length(state), step = 1)]  
    count = count + length(state);

    # state = string(state...) # Return State as a string rather than char array
    return state, count
end

function vector_to_state(vector, state::Number, count)
    """
        Internal function that converts vector to a Number state
    """

    state = float(vector[count + 1])  
    count = count + 1

    return state, count
end

function vector_to_state(vector, state::Array{<:Number}, count)
    """
        Internal function that converts vector to a Number array state
    """
    
    state_copy = deepcopy(state)
    state_copy[:] = float.(vector[count .+ range(1, length(state_copy), step = 1)])
    count = count + length(state_copy);

    return state_copy, count
end

function vector_to_state(vector, state::Array, count)
    """
        Internal function that converts vector to a vector of structs state
    """

    state_copy = deepcopy(state) # So we can setproperty! without modifying original

	# If it's a struct, convert and store in each field.       
    for n = 1:length(state_copy)
        temp, count = vector_to_state(vector, state_copy[n], count)
        state_copy[n] = temp
    end

    return state_copy, count
end

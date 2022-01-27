# [src/state_to_vector.jl]

""" 
Build a state vector by recursively moving through the elements of a more
complex state (matrix, cell array, or struct), saving the numerical 
values along the way in a vector. This is used in odehybrid to convert a
complex series of states into a vector to be used for continuous state
updates. The reverse function is vector_to_state.

vector = state_to_vector(state);

Inputs:

    state:   A numeric, array, or struct array type, the contents of
                which consist of numeric, array, or struct array types, etc.
                (Structs MUST be mutable)

Outputs:

    vector:  A vector containing the numerical values from the state (doubles)

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

See also: vector_to_state.

Online doc: http://www.anuncommonlab.com/doc/odehybrid/state_to_vector.html

Copyright 2014 An Uncommon Lab
"""
function state_to_vector(state)

    # If it shows up here, it should be a struct
    vector = [];
    fields = fieldnames(typeof(state))
    for n = 1:1 # 1:length(state)
        for k = 1:length(fields)
            vector = [vector; state_to_vector(getproperty(state, fields[k]))]
        end
    end

    return vector
end


function state_to_vector(state::Number)
    """
        Converts a number state to a float vector
    """
    # Converts a number to a float and returns it 
    vector = float(state)

    return vector
end 

function state_to_vector(state::Array{<:Number})
    """
        Converts a number array state to a float array vector
    """
    return float.(state[:])
end

function state_to_vector(state::Char)
    """
        Converts a Char state to a float vector
    """
    vector = float(state)

    return vector
end

function state_to_vector(state::Array{Char})
    """
        Converts a Char arary state to a float array vector
    """
    return float.(state[:])
end

function state_to_vector(state::String)
    """
        Converts a string state to a vector
    """
    # Converts a string to a vector
    vector = []
    for i = 1:length(state)
        vector = [vector; state_to_vector(state[i])]
    end

    return vector
end

function state_to_vector(state::Array) #{Any, 1})
    """
        Converts each element of an Array to the appropriate value and returns in a vector
    """
    vector = []
    for i = 1:length(state)
        vector = [vector; state_to_vector(state[i])]
    end

    return vector
end

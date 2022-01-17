# Set of functions to convert states of various types into vectors
function state_to_vector(state)
    """ state_to_vector

        Build a state vector by recursively moving through the elements of a more
        complex state (matrix, cell array, or struct), saving the numerical 
        values along the way in a vector. This is used in odehybrid to convert a
        complex series of states into a vector to be used for continuous state
        updates. The reverse function is vector_to_state.

        vector = state_to_vector(state);
        
        Inputs:
        
            state:   A numeric, cell array, or struct array type, the contents of
                        which consist of numeric, cell array, or struct array types, etc.
        
        Outputs:

            vector:  A vector containing the numerical values from the state (doubles)

        Example:
        
            
            x = [ [1 3; 4 2], TestStruct(a = [5; 6], bcd = [7:9; 10],  π ]
            v = state_to_vector(x)
            x2 = vector_to_state(v, x)
            x == x2

        See also: vector_to_state.

        Online doc: http://www.anuncommonlab.com/doc/odehybrid/state_to_vector.html

        Copyright 2014 An Uncommon Lab
    """

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
    # Converts a number to a float and returns it 
    vector = float(state)

    return vector
end 

function state_to_vector(state::Array{<:Number})
    # Converts an Array of numbers to floats
    return float.(state[:])
end

function state_to_vector(state::Char)
    # Converts a Char to float and returns
    vector = float(state)

    return vector
end

function state_to_vector(state::Array{Char})
    # Converts an array of Chars to floats
    return float.(state[:])
end

function state_to_vector(state::String)
    # Converts a string to a vector
    vector = []
    for i = 1:length(state)
        vector = [vector; state_to_vector(state[i])]
    end

    return vector
end

function state_to_vector(state::Array) #{Any, 1})
    # Converts each element of an Array to the appropriate value and returns in a vector
    vector = []
    for i = 1:length(state)
        vector = [vector; state_to_vector(state[i])]
    end

    return vector
end

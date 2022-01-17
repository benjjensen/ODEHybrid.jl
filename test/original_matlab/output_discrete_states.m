function discrete_states = output_discrete_states(states, yd0)

    % Make space for the individual outputs.
    discrete_states = cell(1, numel(yd0));

    % For each output...
    for k = 1:numel(yd0)

        % Convert to a matrix with rows containing vectors.
        if isnumeric(yd0{k}) || ischar(yd0{k})

            dim = find([size(yd0{k}), 1] == 1, 1, 'first');
            if dim == 2
                discrete_states{k} = cat(dim, states{:, k}).';
            else
                discrete_states{k} = cat(dim, states{:, k});
            end

        % Convert to an array of structs.
        elseif isstruct(yd0{k}) && length(yd0{k}) == 1

            discrete_states{k} = cellfun(@(v) v, states(:, k));

        % Give up and just use the cell array of states.
        else

            discrete_states{k} = states(:, k);

        end

    end
        
end

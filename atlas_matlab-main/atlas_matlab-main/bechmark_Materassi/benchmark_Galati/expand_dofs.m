function [out] = expand_dofs(v)

    v_reshaped = reshape(v, 3, []).';

    % Initialize expanded result
    out = zeros(size(v_reshaped,1)*12, 1);
    
    % Fill output
    for i = 1:size(v_reshaped,1)
        base = v_reshaped(i,1) - 3;
        out((i-1)*12 + (1:12)) = base : base + 11;
    end

end
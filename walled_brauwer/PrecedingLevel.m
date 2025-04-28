function [p_pos, q_pos] = PrecedingLevel(p, q)
    % (p, q-1) <- (p, q)
    % (p-1, 0) <- (p, 0)
    if q > 0
        p_pos=p;
        q_pos=q-1;
    else
        p_pos=p-1;
        q_pos=q;
    end
end
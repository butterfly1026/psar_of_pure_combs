function [M, p_pos, q_pos, d] = TensorIdentityInGelfandTsetlinOperatorBasis(p, q, d)
    % (p, q-1) -> (p, q)
    % (p-1, 0) -> (p, 0)
    
    m2=standard_hook_content_formula(p,q,d);
    if q > 0
        p_pos=p;
        q_pos=q-1;
    else
        p_pos=p-1;
        q_pos=q;
    end
    m1 = standard_hook_content_formula(p_pos, q_pos, d);
    
    if     p == 0 && q == 0
        M=[
            1;
        ];
    elseif p == 1 && q == 0
        M=[
            1;
        ];
    elseif p == 2 && q == 0
        M=[
            1;
            1;
        ];
        if d < p+q
            M(2, :) = [];
        end
    elseif p == 1 && q == 1
        M=[
            1;
            1;
        ];
        if d < p+q
            M(2, :) = [];
        end
    elseif p == 3 && q == 0
        M=[
            1, 0;
            1, 0;
            0, 0;
            0, 0;
            0, 1;
            0, 1;
        ];
        if d < p+q
            M(6, :) = [];
        end
    elseif p == 2 && q == 1
        M=[
            1, 0;
            1, 0;
            0, 0;
            0, 0;
            0, 1;
            0, 1;
        ];
        if d < p+q
            M(6, :) = [];
        end
    elseif p == 4 && q == 0
        M=[
            1, 0, 0, 0, 0, 0;
            1, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 1, 0, 0;
            0, 0, 0, 0, 1, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            0, 0, 0, 1, 0, 0;
            0, 0, 0, 0, 1, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 1, 0, 0;
            0, 0, 0, 0, 1, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 1;
            0, 0, 0, 0, 0, 1;
        ];
        if d == 3
            M(24, :) = [];
        end
        if d == 2
            M(:, 6) = [];
            M(15:end, :) = [];
        end
    elseif p == 3 && q == 1
        M=[
            1, 0, 0, 0, 0, 0;
            1, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 1, 0, 0;
            0, 0, 0, 0, 1, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            0, 0, 0, 1, 0, 0;
            0, 0, 0, 0, 1, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 1, 0, 0;
            0, 0, 0, 0, 1, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 1;
            0, 0, 0, 0, 0, 1;
        ];
        if d <= 3
            M(24, :) = [];
        end
        if d == 2
            M(:, 6) = [];
            M([11:14,17,17+3:17+3*2], :) = [];
        end
    else
        error('指定されたp, q, dの組み合わせに対応する行列は定義されていません。');
    end
end
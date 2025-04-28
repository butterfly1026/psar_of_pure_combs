function m = standard_hook_content_formula(p, q, d)
    if     p == 0 && q == 0
        m=[1];
    elseif p == 1 && q == 0
        m=[d];
    elseif p == 0 && q == 1
        m=[d];
    elseif p == 2 && q == 0
        m=[(d^2+d)/2, (d^2-d)/2];
    elseif p == 1 && q == 1
        m=[d^2-1, 1];
    elseif p == 3 && q == 0
        m=[(d^3+3*d^2+2*d)/6, (d^3-d)/3, (d^3-3*d^2+2*d)/6];
    elseif p == 2 && q == 1
        m=[(d^3+d^2-2*d)/2, d, (d^3-d^2-2*d)/2];
    elseif p == 4 && q == 0
        m=[(d^4+6*d^3+11*d^2+6*d)/24, (d^4+2*d^3-d^2-2*d)/8, (d^4-d^2)/12, (d^4-2*d^3-d^2+2*d)/8, (d^4-6*d^3+11*d^2-6*d)/24];
    elseif p == 3 && q == 1
        m=[(1/6)*d^4+(1/2)*d^3+(-1/6)*d^2+(-1/2)*d, (1/2)*d^2+(1/2)*d, (1/3)*d^4+(-4/3)*d^2, (1/2)*d^2+(-1/2)*d, (1/6)*d^4+(-1/2)*d^3+(-1/6)*d^2+(1/2)*d];
    else
        error('指定されたp, q, dの組み合わせに対応する standard_hook_content_formula は定義されていません。');
    end
end
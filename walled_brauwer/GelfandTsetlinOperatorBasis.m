function E = GelfandTsetlinOperatorBasis(p, q, d)
    E={};
    
    if     p == 0 && q == 0
        E{1}=[1];
    elseif p == 1 && q == 0
        E{1}=speye(d);
    elseif p == 0 && q == 1
        E{1}=speye(d);
    elseif p == 2 && q == 0 && d == 2
        initialize_tsu_basis_p2_q0_d2
        E=p2_q0_d2_E;
    elseif p == 1 && q == 1 && d == 2
        initialize_tsu_basis_p1_q1_d2
        E=p1_q1_d2_E;
    elseif p == 3 && q == 0 && d == 2
        initialize_tsu_basis_p3_q0_d2
        E=p3_q0_d2_E;
    elseif p == 2 && q == 1 && d == 2
        initialize_tsu_basis_p2_q1_d2
        E=p2_q1_d2_E;
    elseif p == 4 && q == 0 && d == 2
        initialize_tsu_basis_p4_q0_d2
        E=p4_q0_d2_E;
    elseif p == 3 && q == 1 && d == 2
        initialize_tsu_basis_p3_q1_d2
        E=p3_q1_d2_E;
    elseif p == 3 && q == 3 && d == 2
        initialize_tsu_basis_p3_q3_d2
        E=p3_q3_d2_E;
    elseif p == 2 && q == 0 && d == 3
        initialize_tsu_basis_p2_q0_d3
        E=p2_q0_d3_E;
    elseif p == 3 && q == 0 && d == 3
        initialize_tsu_basis_p3_q0_d3
        E=p3_q0_d3_E;
    elseif p == 2 && q == 1 && d == 3
        initialize_tsu_basis_p2_q1_d3
        E=p2_q1_d3_E;
    elseif p == 4 && q == 0 && d == 3
        initialize_tsu_basis_p4_q0_d3
        E=p4_q0_d3_E;
    elseif p == 3 && q == 1 && d == 3
        initialize_tsu_basis_p3_q1_d3
        E=p3_q1_d3_E;
    elseif p == 2 && q == 0 && d == 4
        initialize_tsu_basis_p2_q0_d4
        E=p2_q0_d4_E;
    elseif p == 1 && q == 1 && d == 4
        initialize_tsu_basis_p1_q1_d4
        E=p1_q1_d4_E;
    elseif p == 3 && q == 0 && d == 4
        initialize_tsu_basis_p3_q0_d4
        E=p3_q0_d4_E;
    elseif p == 2 && q == 1 && d == 4
        initialize_tsu_basis_p2_q1_d4
        E=p2_q1_d4_E;
    elseif p == 4 && q == 0 && d == 4
        initialize_tsu_basis_p4_q0_d4
        E=p4_q0_d4_E;
    elseif p == 3 && q == 1 && d == 4
        initialize_tsu_basis_p3_q1_d4
        E=p3_q1_d4_E;
    elseif p == 2 && q == 0 && d == 6
        initialize_tsu_basis_p2_q0_d6
        E=p2_q0_d6_E;
    elseif p == 2 && q == 1 && d == 6
        initialize_tsu_basis_p2_q1_d6
        E=p2_q1_d6_E;
    else
        p
        q
        d
        error('指定されたp, q, dの組み合わせに対応する行列は定義されていません。');
    end
end
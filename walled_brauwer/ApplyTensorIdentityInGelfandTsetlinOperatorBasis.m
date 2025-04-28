function X_ti = ApplyTensorIdentityInGelfandTsetlinOperatorBasis(X, SYS, PQD)
    % (p, q-1) -> (p, q)
    % (p-1, 0) -> (p, 0)
    tmp_pqds = PQD;
    for cur_sys = SYS
        [p_tmp, q_tmp, d_tmp]=tmp_pqds{cur_sys,:};
        [p_pre, q_pre]=PrecedingLevel(p_tmp, q_tmp);
        tmp_pqds(cur_sys, :) = {p_pre, q_pre, d_tmp};
    end
    total_dim = prod(cellfun(@(pqd) DimesionOfGelfandTsetlinOperatorBasis(pqd{:}), num2cell(tmp_pqds, 2)));
    assert(isequal(size(X), [total_dim, 1]), 'X must have size [%d, 1]', total_dim)

    % Initialize
    tmp_matrix = 1;
    tmp_pqds = PQD;

    for cur_sys = SYS
        dims_cur = cellfun(@(pqd) DimesionOfGelfandTsetlinOperatorBasis(pqd{:}), num2cell(tmp_pqds, 2));
        dim_before = prod(dims_cur(1:cur_sys-1));
        dim_after = prod(dims_cur(cur_sys+1:end));

        [tensor_identity_matrix, p_pre, q_pre, d_pre]=TensorIdentityInGelfandTsetlinOperatorBasis(tmp_pqds{cur_sys,:});
        
        % Update
        tmp_matrix = tmp_matrix * Tensor(eye(dim_before), tensor_identity_matrix, eye(dim_after));
        tmp_pqds(cur_sys, :) = {p_pre, q_pre, d_pre};
    end

    X_ti=tmp_matrix * X;
end
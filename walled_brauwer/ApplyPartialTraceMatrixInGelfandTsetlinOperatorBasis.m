function [X_pt, PQD_pt] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(X, SYS, PQD)
    % Calculate the total dimension based on the input PQD
    total_dim = prod(cellfun(@(pqd) DimesionOfGelfandTsetlinOperatorBasis(pqd{:}), num2cell(PQD, 2)));
    assert(isequal(size(X), [total_dim, 1]), 'X must have size [%d, 1]', total_dim)

    % Initialize
    tmp_matrix = 1;
    tmp_pqds = PQD;

    for cur_sys = SYS
        dims_cur = cellfun(@(pqd) DimesionOfGelfandTsetlinOperatorBasis(pqd{:}), num2cell(tmp_pqds, 2));
        dim_before = prod(dims_cur(1:cur_sys-1));
        dim_after = prod(dims_cur(cur_sys+1:end));

        [partial_trace_matrix, p_pos, q_pos, d_pos]=PartialTraceMatrixInGelfandTsetlinOperatorBasis(tmp_pqds{cur_sys,:});
        
        % Update
        tmp_matrix = Tensor(eye(dim_before), partial_trace_matrix, eye(dim_after)) * tmp_matrix;
        tmp_pqds(cur_sys, :) = {p_pos, q_pos, d_pos};
    end

    X_pt=tmp_matrix * X;
    PQD_pt = tmp_pqds;
end
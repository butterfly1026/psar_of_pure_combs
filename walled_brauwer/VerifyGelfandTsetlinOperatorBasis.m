function VerifyGelfandTsetlinOperatorBasis(p,q,d)
    eps=10^-10;

    if q > 0
        p_pos=p;
        q_pos=q-1;
    else
        p_pos=p-1;
        q_pos=q;
    end
    
    E2=GelfandTsetlinOperatorBasis(p,q,d);
    E1=GelfandTsetlinOperatorBasis(p_pos,q_pos,d);
    
    PT=PartialTraceMatrixInGelfandTsetlinOperatorBasis(p,q,d);
    for j=1:length(E2)
        result = 0;
        for i = 1:length(E1)
            result = result + PT(i,j) * E1{i};
        end
        norm_PT_j=norm(PartialTrace(E2{j},(p+q):(p+q), repmat(d, 1, p+q)) - result);
        if norm_PT_j > eps
            error('Mismatch detected: Partial trace calculation and matrix result do not match at index j=%d. Difference norm_PT_j=%.2e exceeds tolerance of %.2e.', j, norm_PT_j, eps);
        end
    end
    
    TI=TensorIdentityInGelfandTsetlinOperatorBasis(p,q,d);
    for j=1:length(E1)
        result = 0;
        for i = 1:length(E2)
            result = result + TI(i,j) * E2{i};
        end
        norm_TI_j=norm(Tensor(E1{j}, speye(d)) - result);
        if norm_TI_j > eps
            error('Mismatch detected: Tensor identity calculation and matrix result do not match at index j=%d. Difference norm_TI_j=%.2e exceeds tolerance of %.2e.', j, norm_TI_j, eps);
        end
    end

    disp('No mismatches detected.')
end
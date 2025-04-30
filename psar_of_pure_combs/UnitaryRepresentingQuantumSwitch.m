function U_QS = UnitaryRepresentingQuantumSwitch()
    % P_c P_t A_O B_O -> A_I B_I F_c F_t
    d=2;

    id=speye(2);
    ket_0=id(:,0+1);
    ket_1=id(:,1+1);

    P_AB=PermutationOperator(d, [3 1 2], 0, 1);
    P_BA=PermutationOperator(d, [2 3 1], 0, 1);

    U_QS_PFAB=Tensor(ket_0 * ket_0', P_AB) + Tensor(ket_1 * ket_1', P_BA); % P_c P_t A_O B_O -> F_c F_t A_I B_I
    
    dim=[2*d,d^2];
    U_QS=Swap(U_QS_PFAB,[1,2],dim,1); % P_c P_t A_O B_O -> A_I B_I F_c F_t
end

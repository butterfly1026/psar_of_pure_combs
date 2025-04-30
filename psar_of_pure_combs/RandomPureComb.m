function U_comb = RandomPureComb(dims, is_ends_identity)
    nports = length(dims);
    assert(mod(nports,2)==0);
    nU=nports/2;

    dims_I = dims(1:2:end);
    dims_O = dims(2:2:end);

    D_I=prod(dims_I);
    D_O=prod(dims_O);

    assert(D_I == D_O)
    d=D_I;

    U_comb = eye(d);
    d_up = 1;
    d_down = d;
    d_a = 1;
    for i_U=1:nU
        d_U = d_a*dims_I(i_U);
        assert(mod(d_down, dims_I(i_U)) == 0);
        d_down = d_down / dims_I(i_U);

        if is_ends_identity && ((i_U == 1) || (i_U == nU))
            U = speye(d_U);
        else
            U = RandomUnitary(d_U);
        end

        assert(d_up*d_U*d_down == d)
        U_embedded = EmbedOperatorInLargerSystem(U, 2, [d_up, d_U, d_down]);
        U_comb = U_embedded * U_comb;

        d_up = d_up * dims_O(i_U);
        assert(mod(d_U, dims_O(i_U)) == 0);
        d_a = d_U / dims_O(i_U);
    end
end

function embedded_X = EmbedOperatorInLargerSystem(X, SYS, DIM)
    % Validate inputs
    assert(size(X,1) == size(X,2), 'Matrix X must be square.')
    assert(size(X,1) == prod(DIM(SYS)), 'The size of X must match the product of dimensions specified in PERM.')

    embedded_X=Tensor(speye(prod(DIM(1:SYS-1))), X, speye(prod(DIM(SYS+1:end))));
end

function d_gto = DimesionOfGelfandTsetlinOperatorBasis(p, q, d)
    if p == 0 && q == 0
        d_gto = 1;
    elseif p == 3 && q == 3
        d_gto = 132;
    else
        d_gto = size(TensorIdentityInGelfandTsetlinOperatorBasis(p, q, d), 1);
    end
end
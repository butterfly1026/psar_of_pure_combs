function result = make_output_base_matrices(dims)
    nports = length(dims);

    if nports == 2
        dims_decomposed=[dims(1)];
        perm_decomposed=[1 2];
    elseif nports == 4
        if dims(1) == dims(2)
            dims_decomposed=[dims(1) dims(3)];
            perm_decomposed=[1 2 3 4 5 6];
        else
            dims_decomposed=[dims(2) dims(1)/dims(2) dims(3)];
            perm_decomposed=[1 3 2 5 4 6];
        end
    elseif nports == 6
        % Assume comb is non-trivial
        dims_decomposed=[dims(2) dims(1)/dims(2) dims(3) dims(5)];
        perm_decomposed=[1 3 2 5 4 7 6 8];
    else
        error('Invalid dims')
    end

    nports_decomposed = length(dims_decomposed);
    base_iso_or_other=cell(nports_decomposed, 2);
    for i_dims_decomposed=1:length(dims_decomposed)
        dim_decomposed=dims_decomposed(i_dims_decomposed);
        base_iso_or_other{i_dims_decomposed, 1}=IsotropicState(dim_decomposed,1);
        base_iso_or_other{i_dims_decomposed, 2}=speye(dim_decomposed^2) - base_iso_or_other{i_dims_decomposed, 1};
    end

    num_dims = size(base_iso_or_other, 1);
    indices = repmat({1:2}, 1, num_dims);
    
    [grids{1:num_dims}] = ndgrid(indices{:});
    
    combinations = cell2mat(cellfun(@(x) x(:), grids, 'UniformOutput', false));
    
    num_combinations = size(combinations, 1);
    
    output_base_matrices=cell(num_combinations);
    for i = 1:num_combinations
        idx = num2cell(combinations(i, :));
        tensor_elements = arrayfun(@(n) base_iso_or_other{n, idx{n}}, 1:num_dims, 'UniformOutput', false);
        output_base_matrices{i} = Tensor(tensor_elements{:});
    end

    output_base_matrices_permuted = arrayfun(@(i) PermuteSystems(output_base_matrices{i}, perm_decomposed, repelem(dims_decomposed, 2)), 1:num_combinations, 'UniformOutput', false);
    result = output_base_matrices_permuted;
end


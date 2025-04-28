function S = SparseArray(varargin)
    % This function creates a sparse matrix from either multiple numeric vectors or combines multiple sparse matrices.

    if nargin == 0
        error('No input provided');
    end

    if all(cellfun(@(x) isnumeric(x) && ~issparse(x), varargin))
        % Case 1: All inputs are numeric vectors
        vector = [varargin{:}];  % Combine inputs into a single vector
        n = length(vector);  % Number of columns is the length of the vector
        nonZeroIndices = find(vector);  % Find indices of non-zero elements
        nonZeroValues = vector(nonZeroIndices);  % Get non-zero values
        S = sparse(1, nonZeroIndices, nonZeroValues, 1, n);  % Create a 1-row sparse matrix

    elseif all(cellfun(@issparse, varargin))
        % Case 2: All inputs are sparse matrices
        tmpmat = [];  % Initialize empty matrix to concatenate rows

        for i = 1:nargin
            currentMatrix = varargin{i};
            tmpmat = [tmpmat; currentMatrix];  % Append each sparse matrix as a new row
        end

        S = tmpmat;  % Resulting sparse matrix after combining all inputs
    else
        error('All inputs must be either all numeric vectors or all sparse matrices.');
    end
end

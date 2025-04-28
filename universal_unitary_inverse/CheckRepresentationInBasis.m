function [solution, difference_norm] = CheckRepresentationInBasis(BaseMatrices, TargetMatrix)
    assert(~isempty(BaseMatrices))
    
    % Initialize variables
    flattenedBasis = [];
    
    % Flatten the base matrices (assuming BaseMatrices is a cell array)
    for i = 1:length(BaseMatrices)
        currentMatrix = BaseMatrices{i};
        flattenedBasis = [flattenedBasis; currentMatrix(:)'];
    end
    
    % Transpose the flattened basis matrices to create coeffMatrix
    coeffMatrix = flattenedBasis';

    % Flatten the target matrix
    flattenedTarget = TargetMatrix(:);
    
    % Solve for the coefficients that represent TargetMatrix in terms of the base matrices
    solution = coeffMatrix \ flattenedTarget;
    
    % Reconstruct the TargetMatrix using the calculated solution
    reconstructedMatrix = zeros(size(TargetMatrix));
    for i = 1:length(BaseMatrices)
        reconstructedMatrix = reconstructedMatrix + solution(i) * BaseMatrices{i};
    end
    
    % Calculate the difference between the original and the reconstructed TargetMatrix
    differenceMatrix = TargetMatrix - reconstructedMatrix;
    difference_norm=norm(differenceMatrix);
end

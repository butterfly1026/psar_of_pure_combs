function vec = BinaryToVector(binary)
    num_qubits = length(binary);
    idx = sum(binary .* (2.^(num_qubits-1:-1:0)));
    vec = sparse(2^num_qubits, 1);
    vec(idx + 1) = 1;
end

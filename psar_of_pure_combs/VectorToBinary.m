function binary = VectorToBinary(vec)
    idx = find(vec, 1) - 1;
    num_qubits = log2(length(vec));
    binary = bitget(idx, num_qubits:-1:1);
end

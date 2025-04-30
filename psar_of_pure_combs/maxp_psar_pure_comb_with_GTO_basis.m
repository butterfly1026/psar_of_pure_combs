% Evaluates the maximum success probability of a PSAR protocol
% for pure K‐comb transformations, using a Gelfand–Tsetlin operator basis.
% 
% Inputs: (wip)
%   U_in            – Input unitary operators representing target pure combs
%   Jout            – Choi operators of target pure combs |U_in>><<U_in|
%   N               – Number of calls of the input pure combs
%   dims            – Port dimensions of target pure combs
%   PORTS_BASE_COMB – Port numbers of the base comb circuit performing PSAR
%   is_Jin_identity – Flag: if true, assume U_in is composed of only identity operator
%
%   Set the following arguments as false unless you are concern about indefinite causal order
%   is_retrival_ico – Flag: whether retrieval circuit uses an "indefinite causal order" structure
%   is_storage_2_slot_no_signal – Flag: whether storage circuit insert non-signalling channel into input pure combs
%   is_switch       – Flag: set true if input is quantum switch-like
%
% Outputs:
%   maxp – maximum success probability achieved by the SDP
%   x, y – optimal SDP variables parameterizing the base comb circuit

function [maxp, x, y]=maxp_psar_pure_comb_with_GTO_basis(U_in,Jout,N,dims,PORTS_BASE_COMB,is_Jin_identity,is_retrival_ico,is_storage_2_slot_no_signal,is_switch)
    k_in=size(U_in,2); % Count the number of Choi operators
    
    assert(isequal(size(size(dims)), [1 2]))
    assert(size(dims, 1) == 1)

    dims_I=dims(1:2:end);
    dims_O=dims(2:2:end);

    D  =prod(dims); % Whole dimension of pure combs to be stored.
    D_I=prod(dims_I);
    D_O=prod(dims_O);
    assert(D_I == D_O)

    % Storage
    dims_S=repmat(dims, 1, N);
    dims_S_port=dims .^ N;
    D_S =D^N;
    D_SI=D_I^N;
    D_SO=D_O^N;

    % Retrieval
    dims_R=dims;
    dims_R_port=dims;
    D_R =D;
    D_RI=D_I;
    D_RO=D_O;

    % Whole circuit
    dims_L=[dims_S, dims_R];
    dims_L_port=dims_S_port .* dims_R_port;
    D_L =D_S *D_R;
    D_LI=D_SI*D_RI;
    D_LO=D_SO*D_RO;

    nports = length(dims);
    input_ports=1:2:nports;
    output_ports=2:2:nports;
    inserting_swap_perm = [1:2:nports, 2:2:nports];
    ntooth=size(PORTS_BASE_COMB, 2) / 2;

    p=N;
    q=1;

    PQD=cell(nports,3); % (p,q,d) for each port number 
    E=cell(nports, 1);  % The basis matrices corresponding to Gelfand–Tsetlin basis for base comb circuit
    num_E=zeros(1, nports);
    for port=1:nports
        [PQD{port,:}]=deal(p, q, dims(port));
        E{port}=GelfandTsetlinOperatorBasis(p, q, dims(port));
        num_E(port)=size(E{port}, 2);
    end
    num_Es_cell=arrayfun(@(x) 1:x, num_E, 'UniformOutput', false);
    E_index_combinations=combinations(num_Es_cell{:});
    num_E_index_combinations=size(E_index_combinations, 1);

    eps_E=1e-13;

    % Determing filename of parameter file
    dims_str = sprintf('%d,', dims(1:end-1));
    dims_str = [dims_str, sprintf('%d', dims(end))];
    if is_Jin_identity
        if is_switch
            terms_file_path = ['./terms/' 'terms' '_dims=' dims_str '_N=' num2str(N) '_switch' '.mat'];
        else
            terms_file_path = ['./terms/' 'terms' '_dims=' dims_str '_N=' num2str(N) '.mat'];
        end
    else
        terms_file_path = ['./terms/' 'terms' '_dims=' dims_str '_N=' num2str(N) '_U_in' '.mat'];
        terms_file_path_identity = ['./terms/' 'terms' '_dims=' dims_str '_N=' num2str(N) '.mat'];
    end

    % Compute comb_ports
    comb_ports_cell=arrayfun(@(port) port:nports:nports*(N+1), 1:nports, 'UniformOutput', false);
    comb_ports=[comb_ports_cell{:}];

    % Compute output_base_matrices
    output_base_matrices = {}; % base matrices for output; This is valid only if input is identity operator where K=2 or reduced to I \otimes U \otimes I where K=3
    output_base_matrices = make_output_base_matrices(dims);
    num_combinations_base_matrices=length(output_base_matrices);

    % Load parameters if parameter file exists
    if exist(terms_file_path, 'file')
        loaded_data = load(terms_file_path);
        if isfield(loaded_data, 'k_in_terms')
            k_in_terms_before = loaded_data.k_in_terms;
        else
            k_in_terms_before = 1;
        end

        if isfield(loaded_data, 'E_index_combination_index_completed')
            E_index_combination_index_completed = loaded_data.E_index_combination_index_completed;
        else
            E_index_combination_index_completed = num_E_index_combinations;
        end

        if isfield(loaded_data, 'U_in_terms')
            U_in_terms = loaded_data.U_in_terms;
        else
            U_in_terms = U_in;
        end

        if isfield(loaded_data, 'J_out_terms')
            J_out_terms = loaded_data.J_out_terms;
        else
            J_out_terms = Jout;
        end

        terms = loaded_data.terms;

        if is_Jin_identity
            if isfield(loaded_data, 'terms_on_basis')
                terms_on_basis = loaded_data.terms_on_basis;
            else
                terms_on_basis = cell(1, num_E_index_combinations);
                for E_index_combination_index=1:size(E_index_combinations, 1)
                    temp_term=terms{1, E_index_combination_index};
                    terms_on_basis{E_index_combination_index} = CheckRepresentationInBasis(output_base_matrices, temp_term);
                end
            end
        else
            loaded_data_identity = load(terms_file_path_identity);
            terms_on_basis=loaded_data_identity.terms_on_basis;
        end

        if ~is_Jin_identity && isequal(dims, [4 2 2 2 2 4]) && ~is_switch
            dims_decomposed=[2 2 2 2];
            perm_decomposed=[1 3 2 5 4 7 6 8];

            nports_decomposed = 2;
            base_iso_or_other=cell(nports_decomposed, 2);
            for i_dims_decomposed=1:2
                dim_decomposed=2;
                base_iso_or_other{i_dims_decomposed, 1}=IsotropicState(dim_decomposed,1);
                base_iso_or_other{i_dims_decomposed, 2}=speye(dim_decomposed^2) - base_iso_or_other{i_dims_decomposed, 1};
            end

            indices = repmat({1:2}, 1, nports_decomposed);
            grids={};
            [grids{1:nports_decomposed}] = ndgrid(indices{:});
            combinations_422224 = cell2mat(cellfun(@(x) x(:), grids, 'UniformOutput', false));
            num_combinations_for_Uin = size(combinations_422224, 1);
            
            output_projectors=cell(1, num_combinations_for_Uin);
            for i_output_projectors = 1:num_combinations_for_Uin
                idx = num2cell(combinations_422224(i_output_projectors, :));
                output_projector_abcd = Tensor(base_iso_or_other{1, idx{1}}, speye((2*2)^2), base_iso_or_other{2, idx{2}});
                output_projectors{i_output_projectors} = PermuteSystems(output_projector_abcd, perm_decomposed, repelem(dims_decomposed, 2));
            end
        end

        if isfield(loaded_data, 'terms_on_basis_for_Uin')
            terms_on_basis_for_Uin = loaded_data.terms_on_basis_for_Uin;
        else
            terms_on_basis_for_Uin=cell(k_in, num_E_index_combinations);

            if ~is_Jin_identity && isequal(dims, [4 2 2 2 2 4]) && ~is_switch
                for i=1:k_in_terms_before
                    for E_index_combination_index=1:size(E_index_combinations, 1)
                        E_index_combination_index

                        temp_term=terms{i, E_index_combination_index};
                        temp_term_on_basis_for_Uin=cell(1,2^2);
                        for i_output_projectors = 1:num_combinations_for_Uin
                            global_temp_term_on_basis_for_Uin_i_output_projectors = output_projectors{i_output_projectors} * temp_term * output_projectors{i_output_projectors};
                            global_temp_term_on_basis_for_Uin_i_output_projectors_permuted=PermuteSystems(global_temp_term_on_basis_for_Uin_i_output_projectors, perm_decomposed, repelem(dims_decomposed, 2), 0, 1) / (trace(output_projectors{i_output_projectors}) / (2*2)^2);
                            temp_term_on_basis_for_Uin{i_output_projectors}=PartialTrace(global_temp_term_on_basis_for_Uin_i_output_projectors_permuted, [1 3], [2^2 4^2 2^2]);
                        end
                        terms_on_basis_for_Uin{i, E_index_combination_index} = temp_term_on_basis_for_Uin;
                    end
                end
            end
        end
    else
        k_in_terms_before = 0;
        E_index_combination_index_completed=0;
        terms=cell(k_in, num_E_index_combinations);
        U_in_terms={};
        J_out_terms={};
        terms_on_basis=cell(1, num_E_index_combinations);
        terms_on_basis_for_Uin=cell(k_in, num_E_index_combinations);
    end

    % Preprocessing in the case where input is switch-like
    if isequal(dims, [4 2 2 2 2 4]) && is_switch
        % Pc Pt Ai Ao Bi Bo Fc Ft
        %  1  2  3  4  5  6  7  8

        % Symmetry: p=3 q=3 d=2
        %           Pt Ao Bo | Ai Bi Ft
        % QSWITCH = |0><0|_PcFc ¥tensor III_PtAoBoAiBiFt + |1><1|_PcFc ¥tensor Perm312_PtAoBoAiBiFt
        
        dim_switch_control=2;
        dim_switch_target=2;
        ntooth_switch=3;
        dims_switch_decomposed=[dim_switch_control repelem(dim_switch_target, ntooth_switch)];
        perm_switch_decomposed=[1 7 2 4 6 3 5 8];
        %                       PcFcPtAoBoAiBiFt

        base_operators_switch=GelfandTsetlinOperatorBasis(ntooth_switch,ntooth_switch,dim_switch_target);
        num_base_operators_switch=size(base_operators_switch,2);

        base_operators_switch_global_permuted_adjoint=cell(1, num_base_operators_switch);
        for i_base_operators_switch = 1:num_base_operators_switch
            temp_base_operator_switch_global_PcFcPtAoBoAiBiFt_adjoint = Tensor(speye(dim_switch_control^2), base_operators_switch{i_base_operators_switch}');
            base_operators_switch_global_permuted_adjoint{i_base_operators_switch} = PermuteSystems(temp_base_operator_switch_global_PcFcPtAoBoAiBiFt_adjoint, perm_switch_decomposed, repelem(dims_switch_decomposed, 2));
        end

        dims_switch_decomposed_global=repelem(dims_switch_decomposed, 2);
        ports_switch_with_symmetry=[2 3 4 5 6 8]; % Pt Ai Ao Bi Bo Ft
    end

    % Compute (\bigotimes_l E^{\lambda_l}_{S_l T_l}) * J_out
    % (J_out is reduced into identity operator where N=2 and I \otimes U \otimes I where N=3)
    % These results are saved into 'terms'
    k_in_terms=k_in_terms_before;
    for i=k_in_terms_before+1:k_in
        if i > length(U_in_terms)
            U_in_terms{i}=sparse(U_in{i});
            J_out_terms{i}=Jout{i};
            save(terms_file_path, 'terms', 'E_index_combination_index_completed', 'k_in_terms', 'U_in_terms', "J_out_terms", "terms_on_basis", "terms_on_basis_for_Uin", '-v7.3')
        end

        identity_choi_ket_vector=MaxEntangled(D_I,1,0);
        U_in_choi_ket_vector=kron(speye(D_I),conj(U_in_terms{i})) * identity_choi_ket_vector;
        U_in_choi_ket_vector_1234=PermuteSystems(U_in_choi_ket_vector, inserting_swap_perm, dims(inserting_swap_perm), 0, 1);
        U_in_choi_ket_vector_1234_tensored=Tensor(U_in_choi_ket_vector_1234, N);
        U_in_choi_ket_vector_1234_tensored_global=Tensor(U_in_choi_ket_vector_1234_tensored, speye(D_R));
        dims_mn_12345678 = [repmat(dims            , 1, N), dims;
                            repmat(ones(size(dims)), 1, N), dims];
        U_in_choi_ket_vector_1234_tensored_global_permuted=PermuteSystems(U_in_choi_ket_vector_1234_tensored_global, comb_ports, dims_mn_12345678);
        U_in_choi_bra_vector_1234_tensored_global_permuted=U_in_choi_ket_vector_1234_tensored_global_permuted';

        batch_size=100;
        for E_index_combination_index_start=E_index_combination_index_completed+1:batch_size:num_E_index_combinations
            E_index_combination_index_stop=min(E_index_combination_index_start+batch_size-1, num_E_index_combinations);

            parfor E_index_combination_index=E_index_combination_index_start:E_index_combination_index_stop
                E_index_combination_index

                temp_term=U_in_choi_ket_vector_1234_tensored_global_permuted;
                for port=1:nports
                    local_E=EmbedOperatorInLargerSystem(E{port}{E_index_combinations{E_index_combination_index, port}}, port, dims_L_port)
                    temp_term=local_E * temp_term;
                end
                temp_term=U_in_choi_bra_vector_1234_tensored_global_permuted * temp_term;
                
                terms{i, E_index_combination_index}=temp_term;
            end

            if is_Jin_identity
                if is_switch
                    parfor E_index_combination_index=E_index_combination_index_start:E_index_combination_index_stop
                        E_index_combination_index
                        temp_term=terms{i, E_index_combination_index};
                        temp_term_on_basis_switch=cell(1,num_base_operators_switch);
                        for i_base_operators_switch = 1:num_base_operators_switch
                            temp_term_on_basis_switch{i_base_operators_switch}=PartialTrace(temp_term * base_operators_switch_global_permuted_adjoint{i_base_operators_switch}, ports_switch_with_symmetry, dims_switch_decomposed_global(perm_switch_decomposed));
                        end
                        terms_on_basis{i, E_index_combination_index} = temp_term_on_basis_switch;
                    end
                else
                    parfor E_index_combination_index=E_index_combination_index_start:E_index_combination_index_stop
                        E_index_combination_index
                        temp_term=terms{i, E_index_combination_index};
                        terms_on_basis{E_index_combination_index} = CheckRepresentationInBasis(output_base_matrices, temp_term);
                    end
                end
            else
                if isequal(dims, [4 2 2 2 2 4]) && ~is_switch
                    % input:
                    % 2 ---- 2
                    % 2 -UU- 2
                    % 2 -UU- 2
                    % 2 ---- 2
                    % の \tensor N
                    %
                    % output:
                    % 2 ---- 2
                    % 2 -AA- 2
                    % 2 -AA- 2
                    % 2 ---- 2
                    %
                    % outputについて、
                    % 一番上と一番下のワイヤーでは対称性を作れるため、base_opに分解できる（symとasym）-> base_iso_or_other
                    % 真ん中の[2 2]は対称性を作れないため、
                    % outputのChoiがbase_op \tensor A \tensor base_opの線型結合になる
                    % terms_on_basis_for_Uin は Aの部分を保存する

                    dims_decomposed=[2 2 2 2];
                    perm_decomposed=[1 3 2 5 4 7 6 8];

                    nports_decomposed = 2;
                    base_iso_or_other=cell(nports_decomposed, 2);
                    for i_dims_decomposed=1:2
                        dim_decomposed=2;
                        base_iso_or_other{i_dims_decomposed, 1}=IsotropicState(dim_decomposed,1);
                        base_iso_or_other{i_dims_decomposed, 2}=speye(dim_decomposed^2) - base_iso_or_other{i_dims_decomposed, 1};
                    end

                    indices = repmat({1:2}, 1, nports_decomposed);
                    grids={};
                    [grids{1:nports_decomposed}] = ndgrid(indices{:});
                    combinations_422224 = cell2mat(cellfun(@(x) x(:), grids, 'UniformOutput', false));
                    num_combinations_for_Uin = size(combinations_422224, 1);
                    
                    output_projectors=cell(1, num_combinations_for_Uin);
                    for i_output_projectors = 1:num_combinations_for_Uin
                        idx = num2cell(combinations_422224(i_output_projectors, :));
                        output_projector_abcd = Tensor(base_iso_or_other{1, idx{1}}, speye((2*2)^2), base_iso_or_other{2, idx{2}});
                        output_projectors{i_output_projectors} = PermuteSystems(output_projector_abcd, perm_decomposed, repelem(dims_decomposed, 2));
                    end
                    
                    parfor E_index_combination_index=E_index_combination_index_start:E_index_combination_index_stop
                        E_index_combination_index
                        temp_term=terms{i, E_index_combination_index};
                        temp_term_on_basis_for_Uin=cell(1,2^2);
                        for i_output_projectors = 1:num_combinations_for_Uin
                            global_temp_term_on_basis_for_Uin_i_output_projectors = output_projectors{i_output_projectors} * temp_term * output_projectors{i_output_projectors};
                            global_temp_term_on_basis_for_Uin_i_output_projectors_permuted=PermuteSystems(global_temp_term_on_basis_for_Uin_i_output_projectors, perm_decomposed, repelem(dims_decomposed, 2), 0, 1) / (trace(output_projectors{i_output_projectors}) / (2*2)^2);
                            temp_term_on_basis_for_Uin{i_output_projectors}=PartialTrace(global_temp_term_on_basis_for_Uin_i_output_projectors_permuted, [1 3], [2^2 4^2 2^2]);
                        end
                        terms_on_basis_for_Uin{i, E_index_combination_index} = temp_term_on_basis_for_Uin;
                    end
                end
            end

            E_index_combination_index_completed=E_index_combination_index_stop;
            save(terms_file_path, 'terms', 'E_index_combination_index_completed', 'k_in_terms', 'U_in_terms', "J_out_terms", "terms_on_basis", "terms_on_basis_for_Uin", '-v7.3')
        end

        k_in_terms=i;
        E_index_combination_index_completed=0;
        save(terms_file_path, 'terms', 'E_index_combination_index_completed', 'k_in_terms', 'U_in_terms', "J_out_terms", "terms_on_basis", "terms_on_basis_for_Uin", '-v7.3')
    end

    save(terms_file_path, 'terms', 'E_index_combination_index_completed', 'k_in_terms', 'U_in_terms', "J_out_terms", "terms_on_basis", "terms_on_basis_for_Uin", '-v7.3')


    cvx_begin SDP
        cvx_solver mosek

        variable x(num_E)
        variable y(num_E)
        
        % Comb condition
        expression z
        z=reshape(permute(x+y, nports:-1:1), [], 1);

        z_pt2=z;
        pqd_2=PQD;
        if is_retrival_ico
            % 1 2  3  4  5  6
            % P Ai Ao Bi Bo F
            port_P=1;
            port_Ai=2;
            port_Ao=3;
            port_Bi=4;
            port_Bo=5;
            port_F=6;

            [z_F, pqd_F] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z, port_F, PQD);

            [z_AoF, pqd_AoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_F, port_Ao, pqd_F);
            z_AoF = z_AoF / dims(port_Ao);
            [z_AiAoF, pqd_AiAoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_AoF, port_Ai, pqd_AoF);
            [z_AiAoBoF, pqd_AiAoBoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_AiAoF, port_Bo, pqd_AiAoF);
            z_AiAoBoF = z_AiAoBoF / dims(port_Bo);

            [z_BoF, pqd_BoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_F, port_Bo, pqd_F);
            z_BoF = z_BoF / dims(port_Bo);
            [z_BiBoF, pqd_BiBoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_BoF, port_Bi, pqd_BoF);
            [z_AoBiBoF, pqd_AoBiBoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_BiBoF, port_Ao, pqd_BiBoF);
            z_AoBiBoF = z_AoBiBoF / dims(port_Ao);

            [z_AoBoF, pqd_AoBoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_AoF, port_Bo, pqd_AoF);
            z_AoBoF = z_AoBoF / dims(port_Bo);

            [z_AiAoBiBoF, pqd_AiAoBiBoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_AiAoBoF, port_Bi, pqd_AiAoBoF);
            [z_PAiAoBiBoF, pqd_PAiAoBiBoF] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_AiAoBiBoF, port_P, pqd_AiAoBiBoF);
            z_PAiAoBiBoF = z_PAiAoBiBoF / dims(port_P);

            subject to
                z_BiBoF == ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_AoBiBoF, port_Ao, pqd_BiBoF)
                z_AiAoF == ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_AiAoBoF, port_Bo, pqd_AiAoF)
                z_AiAoBiBoF == ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_PAiAoBiBoF, port_P, pqd_AiAoBiBoF)
                z_F  == ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_BoF, port_Bo, pqd_F) + ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_AoF, port_Ao, pqd_F) - ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_AoBoF, [port_Ao port_Bo], pqd_F)

            z_pt2=z_PAiAoBiBoF;
            pqd_2=pqd_PAiAoBiBoF;
        end
        
        for tooth=ntooth:-1:1
            if is_storage_2_slot_no_signal
                if (mod(tooth, 3) == 0) && floor(tooth / 3) <= N
                    tooth_ns=tooth-1;
                    [z_pt1_ns, pqd_1_ns] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_pt2, PORTS_BASE_COMB{2*tooth_ns}, pqd_2);
                    [z_pt2_ns, pqd_2_ns] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_pt1_ns, PORTS_BASE_COMB{2*tooth_ns-1}, pqd_1_ns);
                    z_pt2_ns = z_pt2_ns / prod(dims(PORTS_BASE_COMB{2*tooth-1}));
                    subject to
                        z_pt1_ns == ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_pt2_ns, PORTS_BASE_COMB{2*tooth_ns-1}, pqd_1_ns)
                end
            end
            [z_pt1, pqd_1] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_pt2, PORTS_BASE_COMB{2*tooth}, pqd_2);
            [z_pt2, pqd_2] = ApplyPartialTraceMatrixInGelfandTsetlinOperatorBasis(z_pt1, PORTS_BASE_COMB{2*tooth-1}, pqd_1);
            z_pt2 = z_pt2 / prod(dims(PORTS_BASE_COMB{2*tooth-1}));
            subject to
                z_pt1 == ApplyTensorIdentityInGelfandTsetlinOperatorBasis(z_pt2, PORTS_BASE_COMB{2*tooth-1}, pqd_1)
        end

        subject to
            z_pt2 == 1
        
        % Positivity
        leaf_to_paths = arrayfun(@(dim) LeafToPaths(p, q, dim), dims, 'UniformOutput', false);
        leaf_indices_combination = cellfun(@(leaf_to_path) 1:size(leaf_to_path, 2), leaf_to_paths, 'UniformOutput', false);
        leaf_combinations=combinations(leaf_indices_combination{:});
        for leaf_combination_index=1:size(leaf_combinations, 1)
            paths = arrayfun(@(port) leaf_to_paths{port}{leaf_combinations{leaf_combination_index,port}}, 1:nports, 'UniformOutput', false);
            path_index_combination = cellfun(@(path) 1:size(path, 1), paths, 'UniformOutput',false);
            path_matrix_dimension_combinations=combinations(path_index_combination{:});
            num_path_matrix_dimension_combinations=size(path_matrix_dimension_combinations, 1);

            expression A(num_path_matrix_dimension_combinations, num_path_matrix_dimension_combinations)
            expression B(num_path_matrix_dimension_combinations, num_path_matrix_dimension_combinations)
            for i=1:num_path_matrix_dimension_combinations
                for j=1:num_path_matrix_dimension_combinations
                    path_combination=arrayfun(@(port) paths{port}(path_matrix_dimension_combinations{i, port},path_matrix_dimension_combinations{j, port}), 1:nports, 'UniformOutput', false);
                    A(i,j)=x(path_combination{:});
                    B(i,j)=y(path_combination{:});
                end
            end

            subject to
                A == semidefinite(num_path_matrix_dimension_combinations);
                B == semidefinite(num_path_matrix_dimension_combinations);
        end

        % PSAR condition
        variable maxp
        maximise maxp

        if isequal(dims, [4 2 2 2 2 4]) && is_switch
            i=1;
            Jout_1234=PermuteSystems(J_out_terms{i}, inserting_swap_perm, dims(inserting_swap_perm), 0, 1);

            expression output_sum_on_basis_switch(dim_switch_control^2, dim_switch_control^2, num_base_operators_switch)

            for E_index_combination_index=1:size(E_index_combinations, 1)
                indices_cell=num2cell(E_index_combinations{E_index_combination_index,:});

                term_on_basis = terms_on_basis{E_index_combination_index};
                for i_base_operators_switch = 1:num_base_operators_switch
                    term_on_basis{i_base_operators_switch}(abs(term_on_basis{i_base_operators_switch}) < eps_E) = 0; % Replace "junk" with 0
                    output_sum_on_basis_switch(:,:,i_base_operators_switch) = output_sum_on_basis_switch(:,:,i_base_operators_switch) + x(indices_cell{:}) * term_on_basis{i_base_operators_switch};
                end
            end

            for i_base_operators_switch=1:num_base_operators_switch
                temp_Jout_on_basis_switch=PartialTrace(Jout_1234 * base_operators_switch_global_permuted_adjoint{i_base_operators_switch}, ports_switch_with_symmetry, dims_switch_decomposed_global(perm_switch_decomposed));
                temp_Jout_on_basis_switch(abs(temp_Jout_on_basis_switch) < eps_E) = 0; % Replace "junk" with 0

                subject to
                    output_sum_on_basis_switch(:,:,i_base_operators_switch) == maxp * temp_Jout_on_basis_switch;
            end
        else
            expression output_sum_on_basis(num_combinations_base_matrices)

            for E_index_combination_index=1:size(E_index_combinations, 1)
                E_index_combination_index
                indices_cell=num2cell(E_index_combinations{E_index_combination_index,:});

                term_on_basis = terms_on_basis{E_index_combination_index};
                term_on_basis(abs(term_on_basis) < eps_E) = 0; % Replace "junk" with 0

                for i_output_basis=1:num_combinations_base_matrices
                    output_sum_on_basis(i_output_basis) = output_sum_on_basis(i_output_basis) + x(indices_cell{:}) * term_on_basis(i_output_basis);
                end
            end

            subject to
                output_sum_on_basis(1) / D_I == maxp
            for i_output_basis=2:num_combinations_base_matrices
                output_sum_on_basis(i_output_basis) == 0
            end
        end
        
        if ~is_Jin_identity
            disp('not is_Jin_identity')
            if isequal(dims, [4 2 2 2 2 4]) && ~is_switch
                for i=1:k_in
                    expression output_sum_on_basis_for_Uin(4^2,4^2,2^2)
        
                    for E_index_combination_index=1:size(E_index_combinations, 1)
                        indices_cell=num2cell(E_index_combinations{E_index_combination_index,:});
                        term_on_basis_for_Uin=terms_on_basis_for_Uin{i, E_index_combination_index};
                        for i_output_projectors = 1:num_combinations_for_Uin
                            output_sum_on_basis_for_Uin(:,:,i_output_projectors) = output_sum_on_basis_for_Uin(:,:,i_output_projectors) + x(indices_cell{:}) * term_on_basis_for_Uin{i_output_projectors};
                        end
                    end
        
                    Jout_1234=PermuteSystems(J_out_terms{i}, inserting_swap_perm, dims(inserting_swap_perm), 0, 1);
                    Jout_1234_permuted=PermuteSystems(Jout_1234, perm_decomposed, repelem(dims_decomposed, 2), 0, 1);
                    Jout_1234_A=PartialTrace(Jout_1234_permuted, [1 3], [2^2 4^2 2^2])  / (2*2);
    
                    subject to
                        output_sum_on_basis_for_Uin(:,:,1) / (2*2) == maxp * Jout_1234_A;
    
                    for i_output_projectors=2:num_combinations_for_Uin
                        subject to
                            output_sum_on_basis_for_Uin(:,:,i_output_projectors) == 0;
                    end
                end
            else
                for i=1:k_in
                    expression output_sum(D_R, D_R)
        
                    for E_index_combination_index=1:size(E_index_combinations, 1)
                        indices_cell=num2cell(E_index_combinations{E_index_combination_index,:});
                        term=terms{i, E_index_combination_index};
        
                        output_sum = output_sum + x(indices_cell{:}) * term;
                    end
        
                    Jout_1234=PermuteSystems(J_out_terms{i}, inserting_swap_perm, dims(inserting_swap_perm), 0, 1);
                    subject to
                        output_sum == maxp * Jout_1234;
                end
            end
        end
        
    cvx_end
end

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Start of Adjustable Parameters  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=1;              % Number of uses 'N'
        is_Jin_identity=true; % If true, target pure comb is supposed to be only identity comb. Then k is set to be 1 automatically
        k=1; % Number of candidates of target pure combs
        dimsList = {[4 2 2 4]}; % List of port dimensions of target pure combs to be SARed

        is_switch=false; % Set false unless you assume input is switch-like instead of pure combs; If true, then set dimsList as {[4 2 2 2 2 4]}
        
        % storage_kind ?
        % 1 for sequential, 2 for sequential with
        % inserting SWAP, 3 for parallel with inserting SWAP, 4 for
        % no-signalling (for indefinite causal order) 4 ((2 2) no-signalling (2 2)) 4
        % Inserting SWAP (2,3) means pure combs turn into unitary staircases in storage
        storage_kinds=1:3;
        % storage_kinds=2:3;
        % storage_kinds=2:4;
        % storage_kinds=4;
        
        % retrieval_kind ?
        % 1 for comb, 2 for channel, 7 for 2-slot indefinite causal order
        % Channel retrieval (2) means unitary staircases are retrieved
        retrieval_kinds = 1:2;
        % retrieval_kinds = 7:7;
        % retrieval_kinds = [1 2 7];
        % retrieval_kinds = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   End of Adjustable Parameters   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if is_Jin_identity
    k=1;
end

for index_dims=1:length(dimsList)
    dims = dimsList{index_dims};

    nports = length(dims);
    assert(mod(nports,2)==0);
    nU=nports/2;
    
    dims_I = dims(1:2:end);
    dims_O = dims(2:2:end);
    
    D  =prod(dims);
    D_I=prod(dims_I);
    D_O=prod(dims_O);
    assert(D_I == D_O)
    d=D_I;
    
    max_ent_times_d=IsotropicState(d,1)*d;
    
    Jout={};
    U_in={};
    for i=1:k
        if is_Jin_identity
            if is_switch
                U=UnitaryRepresentingQuantumSwitch();
            else
                U=speye(d);
            end
        else
            U=RandomPureComb(dims, true);
        end
        U_in{i}=U;
        
        Jout{i}=kron(speye(d),U) * max_ent_times_d * kron(speye(d),U');
    end
    
        
    for storage_kind=storage_kinds
        for retrieval_kind=retrieval_kinds
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%  Output setting starts %%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            current_time = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss-SSS');
            time_str = char(current_time);
            dims_str = sprintf('%d,', dims(1:end-1));
            dims_str = [dims_str, sprintf('%d', dims(end))];
            if is_switch
                file_path = ['./results/sar_result' '_dims=' dims_str '_N=' num2str(N) '_switch' '_storage_kind=' num2str(storage_kind) '_retrieval_kind=' num2str(retrieval_kind) '_' time_str '.txt'];
            else
                file_path = ['./results/sar_result' '_dims=' dims_str '_N=' num2str(N) '_storage_kind=' num2str(storage_kind) '_retrieval_kind=' num2str(retrieval_kind) '_' time_str '.txt'];
            end
            diary(file_path)
            disp('---------------------------------------------')
            disp('--------------run_psar_of_comb_form_unitary_with_GTO_basis_3 start-------------------')
            disp('---------------------------------------------')
            
            
            is_Jin_identity
            k
            is_switch
            dims
            N
            storage_kind
            retrieval_kind
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%  Output setting ends %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            tic;
            
            nports = length(dims);
            input_ports=1:2:nports;
            output_ports=2:2:nports;
            
            is_storage_2_slot_no_signal=false;
            if storage_kind == 1
                storage_ports_base_comb=[repmat(arrayfun(@(x) {x}, 1:nports), 1, N)];
            elseif storage_kind == 2
                storage_ports_base_comb=[repmat({input_ports,output_ports}, 1, N)];
            elseif storage_kind == 3
                storage_ports_base_comb=[{repmat(input_ports, 1, N)}, {repmat(output_ports, 1, N)}];
            elseif storage_kind == 4
                storage_ports_base_comb=[repmat(arrayfun(@(x) {x}, 1:nports), 1, N)];
                is_storage_2_slot_no_signal=true;
            end
            
            is_retrival_ico=false;
            if retrieval_kind == 1
                retrieval_ports_base_comb_before={[]};
                retrieval_ports_base_comb_after=[{[]}, arrayfun(@(x) {x}, 1:nports)];
            elseif retrieval_kind == 2
                retrieval_ports_base_comb_before={[]};
                retrieval_ports_base_comb_after=[{[]}, {input_ports}, {output_ports}];
            elseif retrieval_kind == 3
                retrieval_ports_base_comb_before={[]};
                retrieval_ports_base_comb_after=[{[2]}, {[3]}, {[]}, {[1]}, {[4]}];
            elseif retrieval_kind == 4
                retrieval_ports_base_comb_before={[]};
                retrieval_ports_base_comb_after=[{[2]}, {[1,3]}, {[4]}];
            elseif retrieval_kind == 5
                retrieval_ports_base_comb_before={[]};
                retrieval_ports_base_comb_after=[{[]}, arrayfun(@(x) {x}, nports:-1:1)];
            elseif retrieval_kind == 6
                retrieval_ports_base_comb_before={[]};
                retrieval_ports_base_comb_after=[{[]}, {output_ports}, {input_ports}];
            elseif retrieval_kind == 7
                retrieval_ports_base_comb_before={[]};
                retrieval_ports_base_comb_after={[]};
                is_retrival_ico=true;
            end
            PORTS_BASE_COMB=[retrieval_ports_base_comb_before, storage_ports_base_comb, retrieval_ports_base_comb_after];
            PORTS_BASE_COMB
            
            % Find the maximal success probability of PSAR
            maximal_success_probability=maxp_psar_pure_comb_with_GTO_basis(U_in,Jout,N,dims,PORTS_BASE_COMB,is_Jin_identity,is_retrival_ico,is_storage_2_slot_no_signal,is_switch);
            total_time_in_minutes=toc/60;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%  Output results  %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            format long
            maximal_success_probability
            total_time_in_minutes
    
        end
    end
end

diary off

function l2p = LeafToPaths(p, q, d)
    if     p == 0 && q == 0
        l2p = {[1]};
    elseif p == 1 && q == 0
        l2p = {[1]};
    elseif p == 2 && q == 0
        l2p = {[1], [2]};
        if d < p+q
            l2p(2) = [];
        end
    elseif p == 1 && q == 1
        l2p = {[1], [2]};
        if d < p+q
            l2p(2) = [];
        end
    elseif p == 3 && q == 0
        l2p = {[1], [2, 3; 4, 5], [6]};
        if d < p+q
            l2p(3) = [];
        end
    elseif p == 2 && q == 1
        l2p = {[1], [2, 3; 4, 5], [6]};
        if d < p+q
            l2p(3) = [];
        end
    elseif p == 4 && q == 0
        l2p = {[1], [2, 3, 4; 5, 6, 7; 8, 9, 10], [11, 12; 13, 14], [15, 16, 17; 18, 19, 20; 21, 22, 23], [24]};
        if d <= 3
            l2p(5) = [];
        end
        if d <= 2
            l2p(4) = [];
        end
    elseif p == 3 && q == 1
        l2p = {[1], [2, 3, 4; 5, 6, 7; 8, 9, 10], [11, 12; 13, 14], [15, 16, 17; 18, 19, 20; 21, 22, 23], [24]};
        if d <= 3
            l2p(5) = [];
        end
        if d <= 2
            l2p(4) = [];
        end
    else
        error('指定されたp, q, dの組み合わせに対応する行列は定義されていません。');
    end
end
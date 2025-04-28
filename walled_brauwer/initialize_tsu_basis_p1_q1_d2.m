p=1;
q=1;
d=2;

tmp_E_merged=[SparseArray(SparseArray((1/2),0,0,(-1/2)),SparseArray(0,1,0,0), ...
SparseArray(0,0,1,0),SparseArray((-1/2),0,0,(1/2))),SparseArray( ...
SparseArray((1/2),0,0,(1/2)),SparseArray(0,0,0,0),SparseArray(0,0, ...
0,0),SparseArray((1/2),0,0,(1/2)))];

p1_q1_d2_E = mat2cell(tmp_E_merged, d^(p+q), repmat(d^(p+q),1,DimesionOfGelfandTsetlinOperatorBasis(p,q,d)));
p=3;
q=0;
d=2;

tmp_E_merged=[SparseArray(SparseArray(1,0,0,0,0,0,0,0),SparseArray(0,(1/3),( ...
  1/3),0,(1/3),0,0,0),SparseArray(0,(1/3),(1/3),0,(1/3),0,0,0), ...
  SparseArray(0,0,0,(1/3),0,(1/3),(1/3),0),SparseArray(0,(1/3),(1/3) ...
  ,0,(1/3),0,0,0),SparseArray(0,0,0,(1/3),0,(1/3),(1/3),0), ...
  SparseArray(0,0,0,(1/3),0,(1/3),(1/3),0),SparseArray(0,0,0,0,0,0, ...
  0,1)),SparseArray(SparseArray(0,0,0,0,0,0,0,0),SparseArray(0,(2/3) ...
  ,(-1/3),0,(-1/3),0,0,0),SparseArray(0,(-1/3),(1/6),0,(1/6),0,0,0), ...
  SparseArray(0,0,0,(1/6),0,(1/6),(-1/3),0),SparseArray(0,(-1/3),( ...
  1/6),0,(1/6),0,0,0),SparseArray(0,0,0,(1/6),0,(1/6),(-1/3),0), ...
  SparseArray(0,0,0,(-1/3),0,(-1/3),(2/3),0),SparseArray(0,0,0,0,0, ...
  0,0,0)),SparseArray(SparseArray(0,0,0,0,0,0,0,0),SparseArray(0,0, ...
  3.^(-1/2),0,(-1).*3.^(-1/2),0,0,0),SparseArray(0,0,(-1/2).*3.^( ...
  -1/2),0,(1/2).*3.^(-1/2),0,0,0),SparseArray(0,0,0,(1/2).*3.^(-1/2) ...
  ,0,(-1/2).*3.^(-1/2),0,0),SparseArray(0,0,(-1/2).*3.^(-1/2),0,( ...
  1/2).*3.^(-1/2),0,0,0),SparseArray(0,0,0,(1/2).*3.^(-1/2),0,(-1/2) ...
  .*3.^(-1/2),0,0),SparseArray(0,0,0,(-1).*3.^(-1/2),0,3.^(-1/2),0, ...
  0),SparseArray(0,0,0,0,0,0,0,0)),SparseArray(SparseArray(0,0,0,0, ...
  0,0,0,0),SparseArray(0,0,0,0,0,0,0,0),SparseArray(0,3.^(-1/2),( ...
  -1/2).*3.^(-1/2),0,(-1/2).*3.^(-1/2),0,0,0),SparseArray(0,0,0,( ...
  1/2).*3.^(-1/2),0,(1/2).*3.^(-1/2),(-1).*3.^(-1/2),0),SparseArray( ...
  0,(-1).*3.^(-1/2),(1/2).*3.^(-1/2),0,(1/2).*3.^(-1/2),0,0,0), ...
  SparseArray(0,0,0,(-1/2).*3.^(-1/2),0,(-1/2).*3.^(-1/2),3.^(-1/2), ...
  0),SparseArray(0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0)), ...
  SparseArray(SparseArray(0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0, ...
  0,0),SparseArray(0,0,(1/2),0,(-1/2),0,0,0),SparseArray(0,0,0,(1/2) ...
  ,0,(-1/2),0,0),SparseArray(0,0,(-1/2),0,(1/2),0,0,0),SparseArray( ...
  0,0,0,(-1/2),0,(1/2),0,0),SparseArray(0,0,0,0,0,0,0,0), ...
  SparseArray(0,0,0,0,0,0,0,0))];

p3_q0_d2_E = mat2cell(tmp_E_merged, d^(p+q), repmat(d^(p+q),1,DimesionOfGelfandTsetlinOperatorBasis(p,q,d)));
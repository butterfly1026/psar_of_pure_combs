p=3;
q=1;
d=2;

tmp_E_merged=[SparseArray(SparseArray((1/4),0,0,(-1/4),0,(-1/4),0,0,0,(-1/4),0, ...
0,0,0,0,0),SparseArray(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,(1/6),0,(1/6),0,0,(-1/6),(1/6),0,0,(-1/6),0,(-1/6) ...
,0,0),SparseArray((-1/4),0,0,(1/4),0,(1/4),0,0,0,(1/4),0,0,0,0,0, ...
0),SparseArray(0,0,(1/6),0,(1/6),0,0,(-1/6),(1/6),0,0,(-1/6),0,( ...
-1/6),0,0),SparseArray((-1/4),0,0,(1/4),0,(1/4),0,0,0,(1/4),0,0,0, ...
0,0,0),SparseArray(0,0,0,0,0,0,(1/4),0,0,0,(1/4),0,(1/4),0,0,( ...
-1/4)),SparseArray(0,0,(-1/6),0,(-1/6),0,0,(1/6),(-1/6),0,0,(1/6), ...
0,(1/6),0,0),SparseArray(0,0,(1/6),0,(1/6),0,0,(-1/6),(1/6),0,0,( ...
-1/6),0,(-1/6),0,0),SparseArray((-1/4),0,0,(1/4),0,(1/4),0,0,0,( ...
1/4),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(1/4),0,0,0,(1/4),0,( ...
1/4),0,0,(-1/4)),SparseArray(0,0,(-1/6),0,(-1/6),0,0,(1/6),(-1/6), ...
0,0,(1/6),0,(1/6),0,0),SparseArray(0,0,0,0,0,0,(1/4),0,0,0,(1/4), ...
0,(1/4),0,0,(-1/4)),SparseArray(0,0,(-1/6),0,(-1/6),0,0,(1/6),( ...
-1/6),0,0,(1/6),0,(1/6),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0,1,0),SparseArray(0,0,0,0,0,0,(-1/4),0,0,0,(-1/4),0,(-1/4),0,0, ...
(1/4))),SparseArray(SparseArray((3/4),0,0,(1/4),0,(1/4),0,0,0,( ...
1/4),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,(1/6),0,(1/6),0,0,(1/6),(1/6),0,0,(1/6),0,(1/6),0, ...
0),SparseArray((1/4),0,0,(1/12),0,(1/12),0,0,0,(1/12),0,0,0,0,0,0) ...
,SparseArray(0,0,(1/6),0,(1/6),0,0,(1/6),(1/6),0,0,(1/6),0,(1/6), ...
0,0),SparseArray((1/4),0,0,(1/12),0,(1/12),0,0,0,(1/12),0,0,0,0,0, ...
0),SparseArray(0,0,0,0,0,0,(1/12),0,0,0,(1/12),0,(1/12),0,0,(1/4)) ...
,SparseArray(0,0,(1/6),0,(1/6),0,0,(1/6),(1/6),0,0,(1/6),0,(1/6), ...
0,0),SparseArray(0,0,(1/6),0,(1/6),0,0,(1/6),(1/6),0,0,(1/6),0,( ...
1/6),0,0),SparseArray((1/4),0,0,(1/12),0,(1/12),0,0,0,(1/12),0,0, ...
0,0,0,0),SparseArray(0,0,0,0,0,0,(1/12),0,0,0,(1/12),0,(1/12),0,0, ...
(1/4)),SparseArray(0,0,(1/6),0,(1/6),0,0,(1/6),(1/6),0,0,(1/6),0,( ...
1/6),0,0),SparseArray(0,0,0,0,0,0,(1/12),0,0,0,(1/12),0,(1/12),0, ...
0,(1/4)),SparseArray(0,0,(1/6),0,(1/6),0,0,(1/6),(1/6),0,0,(1/6), ...
0,(1/6),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,(1/4),0,0,0,(1/4),0,(1/4),0,0,(3/4))), ...
SparseArray(SparseArray(0,0,0,2.^(-1/2),0,(-1/2).*2.^(-1/2),0,0,0, ...
(-1/2).*2.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0, ...
0,0,0,0,0),SparseArray(0,0,(-1/3).*2.^(-1/2),0,(1/6).*2.^(-1/2),0, ...
0,(1/6).*2.^(-1/2),(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),0,(-1/3) ...
.*2.^(-1/2),0,0),SparseArray(0,0,0,(1/3).*2.^(-1/2),0,(-1/6).*2.^( ...
-1/2),0,0,0,(-1/6).*2.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,(-1/3) ...
.*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),(1/6).*2.^( ...
-1/2),0,0,(1/6).*2.^(-1/2),0,(-1/3).*2.^(-1/2),0,0),SparseArray(0, ...
0,0,(1/3).*2.^(-1/2),0,(-1/6).*2.^(-1/2),0,0,0,(-1/6).*2.^(-1/2), ...
0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(-1/6).*2.^(-1/2),0,0,0,( ...
-1/6).*2.^(-1/2),0,(1/3).*2.^(-1/2),0,0,0),SparseArray(0,0,(-1/3) ...
.*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),(1/6).*2.^( ...
-1/2),0,0,(1/6).*2.^(-1/2),0,(-1/3).*2.^(-1/2),0,0),SparseArray(0, ...
0,(-1/3).*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),(1/6) ...
.*2.^(-1/2),0,0,(1/6).*2.^(-1/2),0,(-1/3).*2.^(-1/2),0,0), ...
SparseArray(0,0,0,(1/3).*2.^(-1/2),0,(-1/6).*2.^(-1/2),0,0,0,( ...
-1/6).*2.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(-1/6).*2.^( ...
-1/2),0,0,0,(-1/6).*2.^(-1/2),0,(1/3).*2.^(-1/2),0,0,0), ...
SparseArray(0,0,(-1/3).*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0,(1/6).* ...
2.^(-1/2),(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),0,(-1/3).*2.^( ...
-1/2),0,0),SparseArray(0,0,0,0,0,0,(-1/6).*2.^(-1/2),0,0,0,(-1/6) ...
.*2.^(-1/2),0,(1/3).*2.^(-1/2),0,0,0),SparseArray(0,0,(-1/3).*2.^( ...
-1/2),0,(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),(1/6).*2.^(-1/2),0, ...
0,(1/6).*2.^(-1/2),0,(-1/3).*2.^(-1/2),0,0),SparseArray(0,0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(-1/2).*2.^(-1/2), ...
0,0,0,(-1/2).*2.^(-1/2),0,2.^(-1/2),0,0,0)),SparseArray( ...
SparseArray(0,0,0,0,0,(1/2).*(3/2).^(1/2),0,0,0,(-1/2).*(3/2).^( ...
1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,(-1/2).*6.^(-1/2),0,0,(1/2).*6.^(-1/2),(1/2).* ...
6.^(-1/2),0,0,(-1/2).*6.^(-1/2),0,0,0,0),SparseArray(0,0,0,0,0,( ...
1/2).*6.^(-1/2),0,0,0,(-1/2).*6.^(-1/2),0,0,0,0,0,0),SparseArray( ...
0,0,0,0,(-1/2).*6.^(-1/2),0,0,(1/2).*6.^(-1/2),(1/2).*6.^(-1/2),0, ...
0,(-1/2).*6.^(-1/2),0,0,0,0),SparseArray(0,0,0,0,0,(1/2).*6.^( ...
-1/2),0,0,0,(-1/2).*6.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0, ...
0,(-1/2).*6.^(-1/2),0,0,0,(1/2).*6.^(-1/2),0,0,0,0,0),SparseArray( ...
0,0,0,0,(-1/2).*6.^(-1/2),0,0,(1/2).*6.^(-1/2),(1/2).*6.^(-1/2),0, ...
0,(-1/2).*6.^(-1/2),0,0,0,0),SparseArray(0,0,0,0,(-1/2).*6.^(-1/2) ...
,0,0,(1/2).*6.^(-1/2),(1/2).*6.^(-1/2),0,0,(-1/2).*6.^(-1/2),0,0, ...
0,0),SparseArray(0,0,0,0,0,(1/2).*6.^(-1/2),0,0,0,(-1/2).*6.^( ...
-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(-1/2).*6.^(-1/2),0,0, ...
0,(1/2).*6.^(-1/2),0,0,0,0,0),SparseArray(0,0,0,0,(-1/2).*6.^( ...
-1/2),0,0,(1/2).*6.^(-1/2),(1/2).*6.^(-1/2),0,0,(-1/2).*6.^(-1/2), ...
0,0,0,0),SparseArray(0,0,0,0,0,0,(-1/2).*6.^(-1/2),0,0,0,(1/2).* ...
6.^(-1/2),0,0,0,0,0),SparseArray(0,0,0,0,(-1/2).*6.^(-1/2),0,0,( ...
1/2).*6.^(-1/2),(1/2).*6.^(-1/2),0,0,(-1/2).*6.^(-1/2),0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,(-1/2).*(3/2).^(1/2),0,0,0,(1/2).*(3/2).^(1/2),0,0,0,0,0)), ...
SparseArray(SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,( ...
-1/3).*2.^(-1/2),0,(-1/3).*2.^(-1/2),0,0,(-1/3).*2.^(-1/2),(-1/3) ...
.*2.^(-1/2),0,0,(-1/3).*2.^(-1/2),0,(-1/3).*2.^(-1/2),0,0), ...
SparseArray(2.^(-1/2),0,0,(1/3).*2.^(-1/2),0,(1/3).*2.^(-1/2),0,0, ...
0,(1/3).*2.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,(1/6).*2.^(-1/2), ...
0,(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),(1/6).*2.^(-1/2),0,0,(1/6) ...
.*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0),SparseArray((-1/2).*2.^(-1/2), ...
0,0,(-1/6).*2.^(-1/2),0,(-1/6).*2.^(-1/2),0,0,0,(-1/6).*2.^(-1/2), ...
0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(-1/6).*2.^(-1/2),0,0,0,( ...
-1/6).*2.^(-1/2),0,(-1/6).*2.^(-1/2),0,0,(-1/2).*2.^(-1/2)), ...
SparseArray(0,0,(1/6).*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0,(1/6).* ...
2.^(-1/2),(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),0,(1/6).*2.^(-1/2) ...
,0,0),SparseArray(0,0,(1/6).*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0,( ...
1/6).*2.^(-1/2),(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),0,(1/6).* ...
2.^(-1/2),0,0),SparseArray((-1/2).*2.^(-1/2),0,0,(-1/6).*2.^(-1/2) ...
,0,(-1/6).*2.^(-1/2),0,0,0,(-1/6).*2.^(-1/2),0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,(-1/6).*2.^(-1/2),0,0,0,(-1/6).*2.^(-1/2), ...
0,(-1/6).*2.^(-1/2),0,0,(-1/2).*2.^(-1/2)),SparseArray(0,0,(1/6).* ...
2.^(-1/2),0,(1/6).*2.^(-1/2),0,0,(1/6).*2.^(-1/2),(1/6).*2.^(-1/2) ...
,0,0,(1/6).*2.^(-1/2),0,(1/6).*2.^(-1/2),0,0),SparseArray(0,0,0,0, ...
0,0,(1/3).*2.^(-1/2),0,0,0,(1/3).*2.^(-1/2),0,(1/3).*2.^(-1/2),0, ...
0,2.^(-1/2)),SparseArray(0,0,(-1/3).*2.^(-1/2),0,(-1/3).*2.^(-1/2) ...
,0,0,(-1/3).*2.^(-1/2),(-1/3).*2.^(-1/2),0,0,(-1/3).*2.^(-1/2),0,( ...
-1/3).*2.^(-1/2),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) ...
,SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),SparseArray( ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,(1/3),0,(-1/6),0,0,(-1/6) ...
,(-1/6),0,0,(-1/6),0,(1/3),0,0),SparseArray(0,0,0,(2/3),0,(-1/3), ...
0,0,0,(-1/3),0,0,0,0,0,0),SparseArray(0,0,(-1/6),0,(1/12),0,0,( ...
1/12),(1/12),0,0,(1/12),0,(-1/6),0,0),SparseArray(0,0,0,(-1/3),0,( ...
1/6),0,0,0,(1/6),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(1/6),0,0,0, ...
(1/6),0,(-1/3),0,0,0),SparseArray(0,0,(-1/6),0,(1/12),0,0,(1/12),( ...
1/12),0,0,(1/12),0,(-1/6),0,0),SparseArray(0,0,(-1/6),0,(1/12),0, ...
0,(1/12),(1/12),0,0,(1/12),0,(-1/6),0,0),SparseArray(0,0,0,(-1/3), ...
0,(1/6),0,0,0,(1/6),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(1/6),0, ...
0,0,(1/6),0,(-1/3),0,0,0),SparseArray(0,0,(-1/6),0,(1/12),0,0,( ...
1/12),(1/12),0,0,(1/12),0,(-1/6),0,0),SparseArray(0,0,0,0,0,0,( ...
-1/3),0,0,0,(-1/3),0,(2/3),0,0,0),SparseArray(0,0,(1/3),0,(-1/6), ...
0,0,(-1/6),(-1/6),0,0,(-1/6),0,(1/3),0,0),SparseArray(0,0,0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)) ...
,SparseArray(SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,( ...
1/2).*3.^(-1/2),0,0,(-1/2).*3.^(-1/2),(-1/2).*3.^(-1/2),0,0,(1/2) ...
.*3.^(-1/2),0,0,0,0),SparseArray(0,0,0,0,0,3.^(-1/2),0,0,0,(-1).* ...
3.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,(-1/4).*3.^(-1/2),0,0,( ...
1/4).*3.^(-1/2),(1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,0,0,0), ...
SparseArray(0,0,0,0,0,(-1/2).*3.^(-1/2),0,0,0,(1/2).*3.^(-1/2),0, ...
0,0,0,0,0),SparseArray(0,0,0,0,0,0,(1/2).*3.^(-1/2),0,0,0,(-1/2).* ...
3.^(-1/2),0,0,0,0,0),SparseArray(0,0,0,0,(-1/4).*3.^(-1/2),0,0,( ...
1/4).*3.^(-1/2),(1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,0,0,0), ...
SparseArray(0,0,0,0,(-1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),(1/4).* ...
3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,0,0,0),SparseArray(0,0,0,0,0,( ...
-1/2).*3.^(-1/2),0,0,0,(1/2).*3.^(-1/2),0,0,0,0,0,0),SparseArray( ...
0,0,0,0,0,0,(1/2).*3.^(-1/2),0,0,0,(-1/2).*3.^(-1/2),0,0,0,0,0), ...
SparseArray(0,0,0,0,(-1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),(1/4).* ...
3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,0,0,0),SparseArray(0,0,0,0,0,0,( ...
-1).*3.^(-1/2),0,0,0,3.^(-1/2),0,0,0,0,0),SparseArray(0,0,0,0,( ...
1/2).*3.^(-1/2),0,0,(-1/2).*3.^(-1/2),(-1/2).*3.^(-1/2),0,0,(1/2) ...
.*3.^(-1/2),0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),SparseArray( ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0, ...
(-1/2).*6.^(-1/2),0,(-1/2).*6.^(-1/2),0,0,(-1/2).*6.^(-1/2),(-1/2) ...
.*6.^(-1/2),0,0,(-1/2).*6.^(-1/2),0,(-1/2).*6.^(-1/2),0,0), ...
SparseArray((1/2).*(3/2).^(1/2),0,0,(1/2).*6.^(-1/2),0,(1/2).*6.^( ...
-1/2),0,0,0,(1/2).*6.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0, ...
(-1/2).*6.^(-1/2),0,0,0,(-1/2).*6.^(-1/2),0,(-1/2).*6.^(-1/2),0,0, ...
(-1/2).*(3/2).^(1/2)),SparseArray(0,0,(1/2).*6.^(-1/2),0,(1/2).* ...
6.^(-1/2),0,0,(1/2).*6.^(-1/2),(1/2).*6.^(-1/2),0,0,(1/2).*6.^( ...
-1/2),0,(1/2).*6.^(-1/2),0,0),SparseArray(0,0,(1/2).*6.^(-1/2),0,( ...
1/2).*6.^(-1/2),0,0,(1/2).*6.^(-1/2),(1/2).*6.^(-1/2),0,0,(1/2).* ...
6.^(-1/2),0,(1/2).*6.^(-1/2),0,0),SparseArray((-1/2).*(3/2).^(1/2) ...
,0,0,(-1/2).*6.^(-1/2),0,(-1/2).*6.^(-1/2),0,0,0,(-1/2).*6.^(-1/2) ...
,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(1/2).*6.^(-1/2),0,0,0,(1/2) ...
.*6.^(-1/2),0,(1/2).*6.^(-1/2),0,0,(1/2).*(3/2).^(1/2)), ...
SparseArray(0,0,(-1/2).*6.^(-1/2),0,(-1/2).*6.^(-1/2),0,0,(-1/2).* ...
6.^(-1/2),(-1/2).*6.^(-1/2),0,0,(-1/2).*6.^(-1/2),0,(-1/2).*6.^( ...
-1/2),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0)),SparseArray(SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0),SparseArray(0,0,(1/2).*3.^(-1/2),0,(-1/4).*3.^(-1/2),0,0,( ...
-1/4).*3.^(-1/2),(-1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,(1/2).* ...
3.^(-1/2),0,0),SparseArray(0,0,0,3.^(-1/2),0,(-1/2).*3.^(-1/2),0, ...
0,0,(-1/2).*3.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,(1/2).* ...
3.^(-1/2),0,0,0,(1/2).*3.^(-1/2),0,(-1).*3.^(-1/2),0,0,0), ...
SparseArray(0,0,(-1/2).*3.^(-1/2),0,(1/4).*3.^(-1/2),0,0,(1/4).* ...
3.^(-1/2),(1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),0,(-1/2).*3.^( ...
-1/2),0,0),SparseArray(0,0,(-1/2).*3.^(-1/2),0,(1/4).*3.^(-1/2),0, ...
0,(1/4).*3.^(-1/2),(1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),0,(-1/2) ...
.*3.^(-1/2),0,0),SparseArray(0,0,0,(-1).*3.^(-1/2),0,(1/2).*3.^( ...
-1/2),0,0,0,(1/2).*3.^(-1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0, ...
(-1/2).*3.^(-1/2),0,0,0,(-1/2).*3.^(-1/2),0,3.^(-1/2),0,0,0), ...
SparseArray(0,0,(1/2).*3.^(-1/2),0,(-1/4).*3.^(-1/2),0,0,(-1/4).* ...
3.^(-1/2),(-1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,(1/2).*3.^( ...
-1/2),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0)),SparseArray(SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0),SparseArray(0,0,0,0,(1/4),0,0,(-1/4),(-1/4),0,0,(1/4),0,0,0, ...
0),SparseArray(0,0,0,0,0,(1/2),0,0,0,(-1/2),0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,(1/2),0,0,0,(-1/2),0,0,0,0,0),SparseArray( ...
0,0,0,0,(-1/4),0,0,(1/4),(1/4),0,0,(-1/4),0,0,0,0),SparseArray(0, ...
0,0,0,(-1/4),0,0,(1/4),(1/4),0,0,(-1/4),0,0,0,0),SparseArray(0,0, ...
0,0,0,(-1/2),0,0,0,(1/2),0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,( ...
-1/2),0,0,0,(1/2),0,0,0,0,0),SparseArray(0,0,0,0,(1/4),0,0,(-1/4), ...
(-1/4),0,0,(1/4),0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0,0,0)),SparseArray(SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) ...
,SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,( ...
1/3),0,(-1/6),0,0,(1/6),(-1/6),0,0,(1/6),0,(-1/3),0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,( ...
-1/6),0,(1/12),0,0,(-1/12),(1/12),0,0,(-1/12),0,(1/6),0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,(1/6),0,(-1/12),0,0,( ...
1/12),(-1/12),0,0,(1/12),0,(-1/6),0,0),SparseArray(0,0,(-1/6),0,( ...
1/12),0,0,(-1/12),(1/12),0,0,(-1/12),0,(1/6),0,0),SparseArray(0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0,0,0),SparseArray(0,0,(1/6),0,(-1/12),0,0,(1/12),(-1/12),0,0,( ...
1/12),0,(-1/6),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,(-1/3),0,(1/6),0,0,(-1/6),(1/6),0,0,(-1/6),0,(1/3) ...
,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),SparseArray(SparseArray(0,0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0),SparseArray(0,0,0,0,(1/2).*3.^(-1/2),0,0,(1/2).*3.^(-1/2),( ...
-1/2).*3.^(-1/2),0,0,(-1/2).*3.^(-1/2),0,0,0,0),SparseArray(0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,(-1/4).*3.^(-1/2), ...
0,0,(-1/4).*3.^(-1/2),(1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),0,0,0, ...
0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,(1/4).*3.^(-1/2),0, ...
0,(1/4).*3.^(-1/2),(-1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,0,0, ...
0),SparseArray(0,0,0,0,(-1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),( ...
1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0),SparseArray(0,0,0,0,(1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),( ...
-1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,0,0,0),SparseArray(0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,(-1/2).*3.^(-1/2), ...
0,0,(-1/2).*3.^(-1/2),(1/2).*3.^(-1/2),0,0,(1/2).*3.^(-1/2),0,0,0, ...
0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0)),SparseArray(SparseArray(0,0,0,0,0,0,0, ...
0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,(1/2).*3.^(-1/2),0,(-1/4) ...
.*3.^(-1/2),0,0,(1/4).*3.^(-1/2),(-1/4).*3.^(-1/2),0,0,(1/4).*3.^( ...
-1/2),0,(-1/2).*3.^(-1/2),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0, ...
0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,(1/2).*3.^(-1/2),0,(-1/4).*3.^(-1/2),0,0,(1/4).* ...
3.^(-1/2),(-1/4).*3.^(-1/2),0,0,(1/4).*3.^(-1/2),0,(-1/2).*3.^( ...
-1/2),0,0),SparseArray(0,0,(-1/2).*3.^(-1/2),0,(1/4).*3.^(-1/2),0, ...
0,(-1/4).*3.^(-1/2),(1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,(1/2) ...
.*3.^(-1/2),0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,( ...
-1/2).*3.^(-1/2),0,(1/4).*3.^(-1/2),0,0,(-1/4).*3.^(-1/2),(1/4).* ...
3.^(-1/2),0,0,(-1/4).*3.^(-1/2),0,(1/2).*3.^(-1/2),0,0), ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),SparseArray( ...
SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0, ...
0,0,(1/4),0,0,(1/4),(-1/4),0,0,(-1/4),0,0,0,0),SparseArray(0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0,0),SparseArray(0,0,0,0,(1/4),0,0,(1/4),(-1/4),0,0,(-1/4),0,0, ...
0,0),SparseArray(0,0,0,0,(-1/4),0,0,(-1/4),(1/4),0,0,(1/4),0,0,0, ...
0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,(-1/4),0,0,(-1/4),( ...
1/4),0,0,(1/4),0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0, ...
0,0,0,0,0,0,0,0,0,0,0,0,0),SparseArray(0,0,0,0,0,0,0,0,0,0,0,0,0, ...
0,0,0))];

p3_q1_d2_E = mat2cell(tmp_E_merged, d^(p+q), repmat(d^(p+q),1,DimesionOfGelfandTsetlinOperatorBasis(p,q,d)));
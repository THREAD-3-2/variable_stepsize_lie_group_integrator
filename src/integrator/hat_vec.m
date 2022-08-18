function mat=hat_vec(a)
% skew-symmetric matrix with axial vector a
% see eq(28) in Betsch[06]
mat=[    0,     -a(3),   a(2);
    a(3),          0,  -a(1);
    -a(2),       a(1),      0];
end
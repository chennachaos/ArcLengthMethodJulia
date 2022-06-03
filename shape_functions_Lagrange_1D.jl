function shape_functions_Lagrange_1D(NodeNums, XX, p, xi)
### Shape Function Routine for a 1D Lagrange polynomials


N, dN_dxi, d2N_dxi2 = shape_functions_Lagrange_base(p, xi);

nen = p+1;

Jx=0.0;
Jy=0.0;
x=0.0;
for kk=1:nen
    Jx = Jx + dN_dxi[kk] * XX[NodeNums[kk],1];
    Jy = Jy + dN_dxi[kk] * XX[NodeNums[kk],2];
    x = x + N[kk] * XX[NodeNums[kk],1];
end
J = sqrt(Jx*Jx+Jy*Jy);

dN_dx   = (1/J) * dN_dxi;
d2N_dx2 = (1/J^2) * d2N_dxi2;


return  N,dN_dx,d2N_dx2,J,x
end
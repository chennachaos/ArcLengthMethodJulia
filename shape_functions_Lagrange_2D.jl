function shape_functions_Lagrange_2D(nodeNums, XX, degree, xi, eta)
### Shape Function Routine for a 1D Lagrange polynomials

nen  = degree+1;
nlbf = nen*nen;

N, dN_dxi, dN_deta  = shape_functions_Lagrange_Quad(degree, xi, eta);

#display(N);
#println("\n\n");
#display(dN_dxi);
#println("\n\n");
#display(dN_deta);
#println("\n\n");

# Gradient of mapping from parameter space to physical space
Bmat=zeros(2,2);

ii = 1;
for ii=1:nlbf
    xc = XX[nodeNums[ii],1];
    yc = XX[nodeNums[ii],2];

    Bmat[1,1] = Bmat[1,1] + xc * dN_dxi[ii] ;
    Bmat[2,1] = Bmat[2,1] + xc * dN_deta[ii] ;

    Bmat[1,2] = Bmat[1,2] + yc * dN_dxi[ii] ;
    Bmat[2,2] = Bmat[2,2] + yc * dN_deta[ii] ;
end

#display(XX);
#println("");
#display(Bmat);


Jac  = det(Bmat);

Binv = inv(Bmat);

# Compute derivatives of basis functions w.r.t physical coordinates

dN_dx = zeros(nlbf,1);
dN_dy = zeros(nlbf,1);

for ii = 1:nlbf
    dN_dx[ii] = dN_dxi[ii] * Binv[1,1] + dN_deta[ii] * Binv[1,2];
    dN_dy[ii] = dN_dxi[ii] * Binv[2,1] + dN_deta[ii] * Binv[2,2];
end

return  N, dN_dx, dN_dy, Jac

end
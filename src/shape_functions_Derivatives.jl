function shape_functions_Derivatives_1D(NodeNums, XX, p, xi)
### Shape Function Routine for a 1D Lagrange polynomials


N, dN_dxi, d2N_dxi2 = shape_functions_Lagrange_1D(p, xi);

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





function shape_functions_Derivatives_2D(elshape, degree, nlbf, xi, eta, xNode, yNode)
    ### Shape Function Derivatives for 2D elements
    
    if(elshape == "TRIA")
        N, dN_dxi, dN_deta  = shape_functions_Lagrange_Tria(degree, xi, eta);
    else
        N, dN_dxi, dN_deta  = shape_functions_Lagrange_Quad(degree, xi, eta);
    end
    
    
    
    # Gradient of mapping from parameter space to physical space
    Bmat=zeros(2,2);
    
    ii = 1;
    for ii=1:nlbf
        xc = xNode[ii];
        yc = yNode[ii];
    
        Bmat[1,1] = Bmat[1,1] + xc * dN_dxi[ii] ;
        Bmat[2,1] = Bmat[2,1] + xc * dN_deta[ii] ;
    
        Bmat[1,2] = Bmat[1,2] + yc * dN_dxi[ii] ;
        Bmat[2,2] = Bmat[2,2] + yc * dN_deta[ii] ;
    end
    
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
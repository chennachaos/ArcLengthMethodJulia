function shape_functions_Lagrange_1D(p, xi)
    ### Shape Function Routine for a 1D Lagrange polynomials

    nen = p+1;

    N        = zeros(nen,1); 
    dN_dxi   = zeros(nen,1);
    d2N_dxi2 = zeros(nen,1);

    if(p == 1)
        N[1] = 0.5*(1.0 - xi);
        N[2] = 0.5*(1.0 + xi);

        dN_dxi[1] = -0.5;
        dN_dxi[2] =  0.5;

        d2N_dxi2[1] = 0.0;
        d2N_dxi2[2] = 0.0;
    elseif(p == 2)
        val = xi*xi;

        N[1] = -0.5 * (xi - val);
        N[2] = 1.0 - val;
        N[3] = 0.5 *(xi + val);

        val = 2.0*xi;

        dN_dxi[1] = -0.5*(1.0 - val);
        dN_dxi[2] = -val;
        dN_dxi[3] =  0.5*(1.0 + val);

        d2N_dxi2[1] =  1.0;
        d2N_dxi2[2] = -2.0;
        d2N_dxi2[3] =  1.0;
    elseif(p == 3)
        fact1 =  9.0/16.0;
        fact2 = 27.0/16.0;
        val1 = xi*xi;

        N[1] = -fact1 * (1 - xi)   * (1/9 - val1);
        N[2] =  fact2 * (1 - val1) * (1/3 - xi);
        N[3] =  fact2 * (1 - val1) * (1/3 + xi);
        N[4] = -fact1 * (1 + xi)   * ( 1/9 - val1);

        val2 = 3.0*val1;

        dN_dxi[1] = -fact1*(-1/9 - 2.0*xi   +  val2);
        dN_dxi[2] =  fact2*(-1   - 2.0/3*xi +  val2);
        dN_dxi[3] =  fact2*(1    - 2.0/3*xi -  val2);
        dN_dxi[4] = -fact1*(1/9  - 2.0*xi   -  val2);

        val2 = 6.0*xi;

        d2N_dxi2[1] = -fact1 * (-2   + val2);
        d2N_dxi2[2] =  fact2 * (-2/3 + val2);
        d2N_dxi2[3] =  fact2 * (-2/3 - val2);
        d2N_dxi2[4] = -fact1 * (-2   - val2);
    elseif(p == 4)
        fact1 = 2.0/3.0;
        fact2 = 8.0/3.0;
        val1 = xi*xi;
        val2 = val1*xi;
        val3 = val2*xi;

        N[1] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
        N[2] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
        N[3] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
        N[4] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
        N[5] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

        val4 = 4.0*val2;

        dN_dxi[1] =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
        dN_dxi[2] = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
        dN_dxi[3] =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
        dN_dxi[4] =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
        dN_dxi[5] = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);

        val4 = 12.0*val1;

        d2N_dxi2[1] =  fact1 * (-0.5  -  6.0*xi  +  val4);
        d2N_dxi2[2] = -fact2 * (-2.0  -  3.0*xi  +  val4);
        d2N_dxi2[3] =    4.0 * (-2.5  -  0.0     +  val4);
        d2N_dxi2[4] =  fact2 * ( 2.0  -  3.0*xi  -  val4);
        d2N_dxi2[5] = -fact1 * ( 0.5  -  6.0*xi  -  val4);
    else
        println("no basis functions defined for this degree");
    end    
    
    return  N, dN_dxi, d2N_dxi2
end





function  shape_functions_Lagrange_Quad(degree, xi1, xi2)

    nlbf = (degree+1)*(degree+1);

    N       = zeros(nlbf, 1);
    dN_dxi1 = zeros(nlbf, 1);
    dN_dxi2 = zeros(nlbf, 1);

    if (degree == 0)
            N[1]       = 1.0;

            dN_dxi1[1] = 0.0;
            dN_dxi2[1] = 0.0;

    elseif (degree == 1)
            """
             v2          u1,v2  u2,v2
             v1          u1,v1  u2,v1
               u1 u2
            """

            v1 = 0.5*(1.0 - xi2);
            v2 = 0.5*(1.0 + xi2);

            dv1 = -0.5;
            dv2 =  0.5;

            u1 = 0.5*(1.0 - xi1);
            u2 = 0.5*(1.0 + xi1);

            du1 = -0.5;
            du2 =  0.5;

            N[1] = u1*v1;
            N[2] = u2*v1;
            N[3] = u2*v2;
            N[4] = u1*v2;

            dN_dxi1[1] = du1*v1;
            dN_dxi1[2] = du2*v1;
            dN_dxi1[3] = du2*v2;
            dN_dxi1[4] = du1*v2;

            dN_dxi2[1] = u1*dv1;
            dN_dxi2[2] = u2*dv1;
            dN_dxi2[3] = u2*dv2;
            dN_dxi2[4] = u1*dv2;

    elseif (degree == 2)
            """
             v2          u1,v2  u3,v2  u2,v2
             v3          u1,v3  u3,v3  u2,v3
             v1          u1,v1  u3,v1  u2,v1
               u1 u3 u2
            """

            v1 = 0.5*(xi2*xi2 - xi2);
            v2 = 0.5*(xi2*xi2 + xi2);
            v3 = 1.0 - xi2*xi2;

            dv1 = 0.5*(2.0*xi2 - 1.0);
            dv2 = 0.5*(2.0*xi2 + 1.0);
            dv3 = -2.0*xi2;

            u1 = 0.5*(xi1*xi1 - xi1);
            u2 = 0.5*(xi1*xi1 + xi1);
            u3 = 1.0 - xi1*xi1;

            du1 = 0.5*(2.0*xi1 - 1.0);
            du2 = 0.5*(2.0*xi1 + 1.0);
            du3 = -2.0*xi1;

            N[1] = u1*v1;
            N[2] = u2*v1;
            N[3] = u2*v2;
            N[4] = u1*v2;
            N[5] = u3*v1;
            N[6] = u2*v3;
            N[7] = u3*v2;
            N[8] = u1*v3;
            N[9] = u3*v3;

            dN_dxi1[1] = du1*v1;
            dN_dxi1[2] = du2*v1;
            dN_dxi1[3] = du2*v2;
            dN_dxi1[4] = du1*v2;
            dN_dxi1[5] = du3*v1;
            dN_dxi1[6] = du2*v3;
            dN_dxi1[7] = du3*v2;
            dN_dxi1[8] = du1*v3;
            dN_dxi1[9] = du3*v3;

            dN_dxi2[1] = u1*dv1;
            dN_dxi2[2] = u2*dv1;
            dN_dxi2[3] = u2*dv2;
            dN_dxi2[4] = u1*dv2;
            dN_dxi2[5] = u3*dv1;
            dN_dxi2[6] = u2*dv3;
            dN_dxi2[7] = u3*dv2;
            dN_dxi2[8] = u1*dv3;
            dN_dxi2[9] = u3*dv3;

    else
            printf("no basis functions defined for this degree = %5d \n", p);

    end

    return  N, dN_dxi1, dN_dxi2;
end



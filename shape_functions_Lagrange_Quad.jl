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
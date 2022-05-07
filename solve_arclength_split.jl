function  solve_arclength_split(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, ds)

    psi = 1.0;

    FextReduced = Fext[assy4r];
    FtF = Fext'*Fext;

    #display(ds)
    #display(Du)
    #display(Dl)
    #display(FtF)


    if (loadStep > 1)
        A = (Du'*Du)[1] + (psi*Dl*Dl)*FtF[1] - ds*ds;

        a = 2.0*Du[assy4r]';
        b = 2.0*psi*Dl*FtF[1];
    else
        A = 0.0;
        a = 0.0*Du[assy4r]';
        b = 1.0;
    end
    #println("A")
    #display(A)
    #println("a")
    #display(a)
    #println("b")
    #display(b)

    ### Applying Boundary Conditions

    R = Rglobal[assy4r];

    rNorm = norm(R,2);
    rNorm = sqrt(rNorm*rNorm + A*A);

    println(" Iter : ", iter, " rNorm : ", rNorm);
    du = R*0.0;
    dl = 0.0;
    converged = false;

    if (rNorm < 1.0e-6)
       converged = true;

       return  converged, du, dl, du
    end

    #K1 = Kglobal[assy4r,assy4r];
    #LU = lu(sparse(K1));
    #LU = lu(K1);

    ### solve the matrix system
    #duu = LU.L\(LU.P*FextReduced);
    #du1 = LU.U\duu;
    ##println("du1 = ")
    ##display(du1)

    #duu = LU.L\(LU.P*R);
    #du2 = LU.U\duu;
    ##println("du2 = ")
    ##display(du2)

    K1 = sparse(Kglobal[assy4r,assy4r]);

    du1 = K1\FextReduced;
    du2 = K1\R;

    du2 = -du2; # this is because the Residual is added to the RHS

    dl = ((a*du2)[1] - A)/(b+(a*du1)[1]);

    #println("dl = ")
    #display(dl)

    du = -du2 + dl*du1;

return  converged, du, dl, du1
end
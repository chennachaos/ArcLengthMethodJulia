

module  NewtonRaphson

using LinearAlgebra
using SparseArrays
using DelimitedFiles
using Plots

include("./Utilities.jl")
using .Utilities
include("./Elements.jl")
using .Elements



export main_newtonraphson


function main_newtonraphson()



#######################

fname = "input_plate_square_Q4-nelem2.txt"
#fname = "input_plate_square_Q4-nelem20.txt"
#fname = "input_plate_square_Q9-nelem10.txt"

#fname = "input_Sze_Ex1_Q4-nelem10x1.txt"
#fname = "input_Sze_Ex3_Q4-nelem10x80.txt"

#fname = "input_Truss_2D_3members_model1.txt";
#fname = "input_Truss_3D_2members.txt";
#fname = "input_Truss_3D_12members.txt";
#fname = "input_LeeFrame-nelem10.txt";
#fname = "input_LeeFrame-nelem20.txt";
#fname = "input_arch-215deg.txt";
#fname = "input_Arch_semicircle-nelem50-sym.txt";
#fname = "input_Arch_semicircle-nelem50-unsym.txt";
#fname = "input-beamEndMoment-nelem10.txt";
#fname = "lattice2.txt"

ndim, ndof, nnode, nelem, nodecoords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist = processfile(fname)

#display(LM)

npElem = size(elemConn)[2] - 2;

disp      = zeros(neq,1);

#display(disp)

dispPrev  = deepcopy(disp);
dispPrev2 = deepcopy(disp);
dispPrev3 = deepcopy(disp);
dispPrev4 = deepcopy(disp);

#dispPrev *= 0.0;
#dispPrev2 *= 0.0;


Kglobal = zeros(neq,neq);
Rglobal = zeros(neq,1);


bf = [0.0; 0.0];

loadincrmax     = deepcopy(loadincr);
loadincrmin     = deepcopy(loadincr);
loadincrmin     = loadincrmin*16.0;
loadfactor      = deepcopy(loadincr);
loadfactorPrev2 = 0.0;
loadfactorPrev  = 0.0;

converged = false;
convergedPrev = false;

loadStepConverged = 0;
output = disp[outputlist]';
llist  = 0.0;

dispFull = [disp];

for  loadStep=1:maxloadSteps
    println("load step = ", loadStep);
    println("load incr = ", loadincr);

    loadfactor = loadfactorPrev + loadincr;

    if (loadStep > 1)
      DsFactor1 = 1.0;
      #DsFactor1  = loadfactor/loadfactorPrev;
      #display(DsFactor1)
      disp       = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2;
    end

    #println("")
    #display(dispPrev2)
    #println("")
    #display(dispPrev)
    #println("aaaaaaaaa")
    #display(disp)

    convergedPrev = converged;
    converged = false;

    for iter = 1:2
        Kglobal *= 0.0;
        Rglobal *= 0.0;

        if (ndim == 2)
          if (ndof == 2) # Truss element
            for e = 1:nelem
                Klocal, Flocal = Truss_2D_model1(elemData, elemConn, e, nodecoords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
            end
          else # Beam element
            for e = 1:nelem
                Klocal, Flocal = GeomExactBeam_2D(elemData, elemConn, e, nodecoords, disp, bf);
                #display(Klocal)
                #display(Flocal)

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
            end
          end
        else
          if (ndof == 3) # Truss element
            for e = 1:nelem
                #Klocal, Flocal = Truss_3D_model2(elemData, elemConn, e, nodecoords, disp, bf);
                nodeNums = elemConn[e,3:end];

                #Klocal, Flocal = Mindlinplate_Linear(elemData, nodeNums, e, nodecoords, disp, bf);
                Klocal, Flocal = Mindlinplate_NonLinear_Model1(elemData, nodeNums, e, nodecoords, disp, bf);

                #display(Flocal)

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
            end
          end
        end

        Rglobal = Rglobal + loadfactor*Fext;

        #display(Kglobal)
        #display(Rglobal)

        R = Rglobal[assy4r];

        #display(R)

        rNorm = norm(R,2);

        println(" Iter : ", iter, " rNorm : ", rNorm);
    
        if (rNorm < 1.0e-5)
           converged = true;
           break;
        end

        #display(Kglobal[assy4r,assy4r]);
        #println("\n\n");

        K1 = sparse(Kglobal[assy4r,assy4r]);

        du = K1\R;

        #display(du);

        disp[assy4r] = disp[assy4r] + du;
    end

    if (converged)
      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;

      dispPrev2 = dispPrev;
      dispPrev  = disp;

      loadincr = minimum([maximum([2.0*loadincr, loadincrmin]), loadincrmax]);

      loadStepConverged += 1;

      output = [output; disp[outputlist]'];
      llist = [llist; loadfactor];

      writeoutputvtk(ndim, nelem, nnode, npElem, ndof, nodecoords, elemConn, disp);

      #pp=plot(nodecoords[:,1], nodecoords[:,2], line=(:black,0.9,5,:dot))
      #for e=1:nelem
      #  n1 = elemConn[e,3];
      #  n2 = elemConn[e,4];
      #  xx = [nodecoords[n1,1]+disp[ndof*(n1-1)+1], nodecoords[n2,1]+disp[ndof*(n2-1)+1]];
      #  yy = [nodecoords[n1,2]+disp[ndof*(n1-1)+2], nodecoords[n2,2]+disp[ndof*(n2-1)+2]];
      #  pp = plot!(xx, yy, line=(:black,0.9,1,:solid))
      #end
      #display(pp)
      #current()

    else
      loadincr *= 0.5;
    end
    #current()

    #println("===========================")
    #display(dispPrev2)
    #println("")
    #display(dispPrev)
    #println("")
    #display(disp)
    #println("===========================")

#    waitforbuttonpress
end
#current()
#println()
#display(output)
#println()
#display(llist)
#println()


fileID = open("solution.dat","w");

#write(fileID, (llist, output[:,1], output[:,2]) );

for ii=1:size(llist,1)
  writedlm(fileID, [llist[ii] output[ii,1] output[ii,2] ]);
  #writedlm(fileID,"%12.8f \t %12.8f \t %12.8f", llist[ii], output[ii,1], output[ii,2]);
end

close(fileID)


return
end


end


#main()

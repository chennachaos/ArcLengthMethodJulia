

#module  arclength

using LinearAlgebra
using SparseArrays
using DelimitedFiles

include("./src/Utilities.jl")
using .Utilities
include("./src/Elements.jl")
using .Elements



#export main_arclength


function main_arclength(dirname, fname)


inpfile = dirname * "/" * fname;
ndim, ndof, nnode, nelem, nodecoords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist = processfile(inpfile);

npElem = size(elemConn)[2] - 2;

disp      = zeros(neq,1);


dispPrev  = deepcopy(disp);
dispPrev2 = deepcopy(disp);
dispPrev3 = deepcopy(disp);
dispPrev4 = deepcopy(disp);


Kglobal = zeros(neq,neq);
Rglobal = zeros(neq,1);


bf = [0.0; 0.0];

Ds     = deepcopy(loadincr);
DsPrev = deepcopy(Ds);
DsMax  = deepcopy(Ds);
DsMin  = deepcopy(Ds);

loadfactor      = loadincr;
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

    #println("Ds = ", Ds)
    #println("DsPrev = ", DsPrev)

    if (loadStep > 1)
      DsFactor1  = Ds/DsPrev
      disp       = (1.0+DsFactor1)*dispPrev - DsFactor1*dispPrev2;
      loadfactor = (1.0+DsFactor1)*loadfactorPrev - DsFactor1*loadfactorPrev2;
    end
    #println("load step = ", loadStep);

    Du = disp - dispPrev;
    Dl = loadfactor - loadfactorPrev;
    #display(Du)
    #println(" \n\n\n")
    #display(Dl)

    convergedPrev = converged;
    converged = false;

    for iter = 1:10
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
                Klocal, Flocal = Truss_3D_model2(elemData, elemConn, e, nodecoords, disp, bf);

                Kglobal = Assembly_Matrix(Kglobal,Klocal,LM,e);
                Rglobal = Assembly_Vector(Rglobal,Flocal,LM,e);
            end
          end
        end

        Rglobal = Rglobal + loadfactor*Fext;

        #display(Kglobal)
        #display(Rglobal)

        #[converged, du, dl] = solve_arclength(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, Ds);
        converged, du, dl = solve_arclength_split(loadStep, neq, iter, Kglobal, Rglobal, dof_force, Fext, assy4r, Du, Dl, Ds);

        #display(du)

        if (converged)
          break;
        end

        disp[assy4r] = disp[assy4r] + du;
        loadfactor = loadfactor + dl;

        Du[assy4r] = Du[assy4r] + du;
        Dl = Dl + dl;
    end

    if (converged)

      # write VTK file
      linestr  = split(fname, ".");
      vtkfilename = dirname * "/" * linestr[1] * "-" * string(loadStep) * ".vtk";
      writeoutputvtk(vtkfilename, ndim, nelem, nnode, npElem, ndof, nodecoords, elemConn, disp);

      if (loadStep == 1)
         Ds = (sqrt(Du'*Du + loadfactor*loadfactor*Fext'*Fext))[1];

         DsMax = Ds;
         DsMin = Ds/1024.0;
      end

      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;
      dispPrev2 = dispPrev;
      dispPrev  = disp;

      DsPrev = Ds;
      if (convergedPrev)
        Ds = minimum([maximum([2.0*Ds, DsMin]), DsMax]);
      end

      dispFull = [dispFull; disp];
      output = [output; disp[outputlist]'];
      llist = [llist; loadfactor];

      loadStepConverged += 1;
    else
      if (convergedPrev)
        Ds = maximum([Ds*0.5, DsMin]);
      else
        Ds = maximum([Ds*0.25, DsMin]);
      end
    end

end

linestr  = split(fname, ".");
slnfilename = dirname * "/" * linestr[1] * "-" * "solution.dat";
fileID = open(slnfilename,"w");

for ii=1:size(llist,1)
  writedlm(fileID, [llist[ii] output[ii,1] output[ii,2] ]);
  #writedlm(fileID,"%12.8f \t %12.8f \t %12.8f", llist[ii], output[ii,1], output[ii,2]);
end

close(fileID)


return
end



main_arclength(ARGS[1], ARGS[2]);





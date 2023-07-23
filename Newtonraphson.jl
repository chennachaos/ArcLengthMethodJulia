

#module  NewtonRaphson

using LinearAlgebra
using SparseArrays
using DelimitedFiles

include("./src/Utilities.jl")
using .Utilities
include("./src/Elements.jl")
using .Elements




function main_newtonraphson(dirname, fname)


inpfile = dirname * "/" * fname;
eltype, ndim, ndof, nnode, nelem, nodecoords, elemConn, elemData, NodeDoFArray, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist = processfile(inpfile);


npElem = size(elemConn)[2] - 2;

disp      = zeros(neq,1);

#display(disp)

dispPrev  = deepcopy(disp);
dispPrev2 = deepcopy(disp);
dispPrev3 = deepcopy(disp);
dispPrev4 = deepcopy(disp);


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

    convergedPrev = converged;
    converged = false;

    # iteration loop
    for iter = 1:10
        Kglobal *= 0.0;
        Rglobal *= 0.0;

        # loop over elements
        for e = 1:nelem

            nodeNums = elemConn[e,3:end];
            forAssyVec = NodeDoFArray[e,:];

            if(eltype == "Beam_EulerBernoulli_2D")
              #
              Klocal, Flocal = Beam_EulerBernoulli_2D(elemData, nodeNums, nodecoords, disp, bf);
              #
            elseif(eltype == "Beam_Timoshenko_2D")
              #
              Klocal, Flocal = Beam_Timoshenko_2D(elemData, nodeNums, nodecoords, disp, bf);
              #
            elseif(eltype == "Beam_GeomExact_2D")
              #
              Klocal, Flocal = Beam_GeomExact_2D(elemData, nodeNums, nodecoords, disp, bf);
              #
            elseif(eltype == "Plate_Mindlin_Linear")
              #
              Klocal, Flocal = Plate_Mindlin_Linear(elemData, nodeNums, nodecoords, disp, bf);
              #
            elseif(eltype == "Plate_Mindlin_NonLinear_Model1")
              #
              Klocal, Flocal = Plate_Mindlin_NonLinear_Model1(elemData, nodeNums, nodecoords, disp, bf);
              #
            elseif(eltype == "Shell_Flat_Linear")
              #
              Klocal, Flocal = Shell_Flat_Linear(elemData, nodeNums, nodecoords, disp, bf);
              #
            elseif(eltype == "Shell_Flat_Linear_Rotation")
              #
              Klocal, Flocal = Shell_Flat_Linear_Rotation(elemData, nodeNums, nodecoords, disp, bf);
              #
            end

            Kglobal = Assembly_Matrix(Kglobal, Klocal, forAssyVec);
            Rglobal = Assembly_Vector(Rglobal, Flocal, forAssyVec);
        end

        Rglobal = Rglobal + loadfactor*Fext;

        R = Rglobal[assy4r];

        rNorm = norm(R,2);

        println(" Iter : ", iter, " rNorm : ", rNorm);
    
        if (rNorm < 1.0e-6)
           converged = true;
           break;
        end

        #fact = 1000.0;
        #println("fact = ", fact);
        #for ii=1:neq
        #
        #  if (abs(Kglobal[ii,ii]) < 1.0e-5)
        #    sum=0.0;
        #    for jj=1:neq
        #      sum += abs(Kglobal[ii,jj]);
        #    end
        #    
        #    if(sum < 1.0e-5)
        #      Kglobal[ii,ii] = fact;
        #    end
        #  end
        #end
        
        #display(Kglobal[assy4r,assy4r]);
        #println("\n\n");

        K1 = sparse(Kglobal[assy4r,assy4r]);

        du = K1\R;

        #display(du);

        disp[assy4r] = disp[assy4r] + du;
    end

    #display(disp);

    if (converged)

      # write VTK file
      linestr  = split(fname, ".");
      vtkfilename = dirname * "/" * linestr[1] * "-" * string(loadStep) * ".vtk";
      writeoutputvtk(vtkfilename, ndim, nelem, nnode, npElem, ndof, nodecoords, elemConn, disp);
  
      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;

      dispPrev2 = dispPrev;
      dispPrev  = disp;

      loadincr = minimum([maximum([2.0*loadincr, loadincrmin]), loadincrmax]);

      loadStepConverged += 1;

      output = [output; disp[outputlist]'];
      llist = [llist; loadfactor];

    else
      loadincr *= 0.5;
    end
end


#linestr  = split(fname, ".");
#slnfilename = dirname * "/" * linestr[1] * "-" * "solution.dat";
#fileID = open(slnfilename,"w");
#
#for ii=1:size(llist,1)
#  writedlm(fileID, [llist[ii] output[ii,:]' ]);
#end
#
#close(fileID)

return
end





main_newtonraphson(ARGS[1],ARGS[2]);






using LinearAlgebra

function Shell_Flat_Linear_Rotation(elmDat, nodeNums, elnum, XX, soln, bf)

af = 1.0;

nlbf   = size(nodeNums)[1];
ndof   = 6;
nsize  = nlbf*ndof;

bforce = elmDat[1];
rho    = elmDat[2];
h      = elmDat[3];
E      = elmDat[4];
nu     = elmDat[5];
#kappa  = elmDat[6];
kappa  = 5.0/6.0;

#println("bforce ", bforce);

G = E/2.0/(1.0+nu);
I = h*h*h/12.0;

# Dm - elasticity matrix for the membrane part
Dm = zeros(3,3);
Dm[1,1] = 1.0;   Dm[1,2] = nu;
Dm[2,1] = nu;    Dm[2,2] = 1.0;
Dm[3,3] = (1.0-nu)*0.5;

Dm = (E*h/(1.0-nu*nu))*Dm;

# Db - elasticity matrix for the bending part
Db = zeros(3,3);
Db = (h*h/12.0)*Dm;


# Ds - elasticity matrix for the shear part
Ds = zeros(2,2);
Ds[1,1] = G;
Ds[2,2] = G;
Ds = (h*kappa)*Ds;

# Combined D matrix for membrane and bending parts
Dmat = zeros(6,6);

Dmat[1:3,1:3] = Dm;
Dmat[4:6,4:6] = Db;


#printf(" %14.12f \t %14.12f \t %14.12f \n ", Db[0][0], Db[0][1], Db[0][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n ", Db[1][0], Db[1][1], Db[1][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n\n\n ", Db[2][0], Db[2][1], Db[2][2]);
#cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;

coordsGlobal = zeros(3,nlbf);

xNode = zeros(nlbf,1);
yNode = zeros(nlbf,1);
zNode = zeros(nlbf,1);

solnElemGlobal = zeros(nsize,1);
solnElemLocal  = zeros(nsize,1);

for ii=1:nlbf
  coordsGlobal[1,ii] = XX[nodeNums[ii],1];
  coordsGlobal[2,ii] = XX[nodeNums[ii],2];
  coordsGlobal[3,ii] = XX[nodeNums[ii],3];

  r1 = ndof*(ii-1)+1;
  r2 = r1+1;
  r3 = r1+2;
  r4 = r1+3;
  r5 = r1+4;
  r6 = r1+5;

  jj = ndof*(nodeNums[ii]-1);

  solnElemGlobal[r1] = soln[jj+1];
  solnElemGlobal[r2] = soln[jj+2];
  solnElemGlobal[r3] = soln[jj+3];
  solnElemGlobal[r4] = soln[jj+4];
  solnElemGlobal[r5] = soln[jj+5];
  solnElemGlobal[r6] = soln[jj+6];
end

#printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
#printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);


# Coordinate Transformation - global to local
#

RotMat, coordsLocal, solnElemLocal = compute_transformation_matrix(coordsGlobal, solnElemGlobal);

#if (elnum in [31, 32, 1023, 1024])
#if (elnum in [8, 64, 4, 60])
if (elnum in [100000])
    println("elnum = \n", elnum);
    display(coordsGlobal);
    println("\n\n");

    display(RotMat);
    println("\n\n");

    display(coordsLocal);
    println("\n\n");
end

xNode = coordsLocal[1,:];
yNode = coordsLocal[2,:];


Klocal=zeros(nsize,nsize); # Local stiffness matrix
Flocal=zeros(nsize,1);   # Local load vector

############################################
#
#   part 1. -- contribution due to bending


if(nlbf == 6)
  elShape = "TRIA";
  degree  = 2;
  nGP = 3;
  gps1, gps2, gws = getGaussPointsTriangle(nGP);
elseif (nlbf == 4)
  elShape = "QUAD";
  degree  = 1;
  nGP = 4;
  gps1, gps2, gws = getGaussPointsQuad(nGP);
elseif (nlbf == 9)
  elShape = "QUAD";
  degree  = 2;
  nGP = 9;
  gps1, gps2, gws = getGaussPointsQuad(nGP);
end


param =[0.0, 0.0];
Bmat = zeros(6,ndof*nlbf);

for gp=1:nGP
    param[1] = gps1[gp];
    param[2] = gps2[gp];

    N, dN_dx, dN_dy, Jac = shape_functions_Derivatives_2D(elShape, degree, nlbf, param[1], param[2], xNode, yNode);

    dvol = gws[gp] * Jac;

    grad    = zeros(3,3);
    phis    = [0.0, 0.0, 0.0];
    dphis   = zeros(3,3);
    strains = zeros(6,1);
    
    for ii=1:nlbf
      bb1 = dN_dx[ii];
      bb2 = dN_dy[ii];
      bb3 = N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;
      r4 = r1+3;
      r5 = r1+4;
      r6 = r1+5;

      grad[1,1] = grad[1,1]   + solnElemLocal[r1] * bb1;
      grad[1,2] = grad[1,2]   + solnElemLocal[r1] * bb2;
      grad[2,1] = grad[2,1]   + solnElemLocal[r2] * bb1;
      grad[2,2] = grad[2,2]   + solnElemLocal[r2] * bb2;
      #grad[3,1] = grad[3,1]   + solnElemLocal[r3] * bb1;
      #grad[3,2] = grad[3,2]   + solnElemLocal[r3] * bb2;

      #phis[1]   = phis[1]     + solnElemLocal[r4] * bb3;
      #phis[2]   = phis[2]     + solnElemLocal[r5] * bb3;
      phis[3]   = phis[3]     + solnElemLocal[r6] * bb3;

      dphis[1,1] = dphis[1,1] + solnElemLocal[r4] * bb1;
      dphis[1,2] = dphis[1,2] + solnElemLocal[r4] * bb2;

      dphis[2,1] = dphis[2,1] + solnElemLocal[r5] * bb1;
      dphis[2,2] = dphis[2,2] + solnElemLocal[r5] * bb2;

      dphis[3,1] = dphis[3,1] + solnElemLocal[r6] * bb1;
      dphis[3,2] = dphis[3,2] + solnElemLocal[r6] * bb2;
    end

    for ii=1:nlbf
        bb1 = dN_dx[ii];
        bb2 = dN_dy[ii];
        bb3 = N[ii];

        r1 = ndof*(ii-1)+1;
        r2 = r1+1;
        r3 = r1+2;
        r4 = r1+3;
        r5 = r1+4;
        r6 = r1+5;

        # Bm - membrane part
        Bmat[1,r1] = bb1;
                             Bmat[2,r2] = bb2;
        Bmat[3,r1] = bb2;    Bmat[3,r2] = bb1;

        # Bb - bending part
        Bmat[4,r4] = bb1;
                             Bmat[5,r5] = bb2;
        Bmat[6,r4] = bb2;    Bmat[6,r5] = bb1;
    end

    # Bm - membrane part
    strains[1] = grad[1,1];
    strains[2] = grad[2,2];
    strains[3] = grad[1,2] + grad[2,1];

    # Bb - bending part
    strains[4] = dphis[1,1];
    strains[5] = dphis[2,2];
    strains[6] = dphis[1,2] + dphis[2,1];

    sigf = Dmat*strains;

    Klocal = Klocal + (transpose(Bmat)*Dmat)*(dvol*Bmat);
    Flocal = Flocal - (transpose(Bmat)*(dvol*sigf));

    temp = 0.0 - phis[3];
    #for ii=1:nlbf
    #  bb3 = (E*h*h*h)*dvol*N[ii];
    #  
    #  r6 = ndof*(ii-1)+6;
    #  
    #  Flocal[r6  ] = Flocal[r6  ] + (bb3*temp) ;
    #  
    #  for jj=1:nlbf
    #    c6 = ndof*(jj-1)+6;
    #  
    #    Klocal[r6, c6] = Klocal[r6, c6] + (bb3*N[jj]) ;
    #  end
    #end

    #for ii=1:nlbf
    #  bb3 = dvol*N[ii];
    #  
    #  r1 = ndof*(ii-1)+1;
    #  r2 = r1+1;
    #  r3 = r1+2;
    #  
    #  #Flocal[r1  ] = Flocal[r1  ] + (bb3*bforce) ;
    #  #Flocal[r3  ] = Flocal[r3  ] + (bb3*bforce) ;
    #end
end #gp

#printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);


############################################
#
#   part 2. -- contribution due to shear

# Gauss points for reduced integration
if(nlbf == 6)
  nGP = 2;
  gps1, gps2, gws = getGaussPointsTriangle(nGP);
elseif (nlbf == 4)
  nGP = 1;
  gps1, gps2, gws = getGaussPointsQuad(nGP);
elseif (nlbf == 9)
  nGP = 4;
  gps1, gps2, gws = getGaussPointsQuad(nGP);
end


Bmat = zeros(2,ndof*nlbf);

for gp=1:nGP
    param[1] = gps1[gp];
    param[2] = gps2[gp];

    N, dN_dx, dN_dy, Jac = shape_functions_Derivatives_2D(elShape, degree, nlbf, param[1], param[2], xNode, yNode);

    dvol = gws[gp] * Jac;

    grad    = zeros(3,3);
    phis    = [0.0, 0.0, 0.0];
    strains = zeros(2,1);
    
    for ii=1:nlbf
      bb1 = dN_dx[ii];
      bb2 = dN_dy[ii];
      bb3 = N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;
      r4 = r1+3;
      r5 = r1+4;
      r6 = r1+5;

      grad[3,1] = grad[3,1] + solnElemLocal[r3] * bb1;
      grad[3,2] = grad[3,2] + solnElemLocal[r3] * bb2;

      phis[1] = phis[1] + solnElemLocal[r4] * bb3;
      phis[2] = phis[2] + solnElemLocal[r5] * bb3;
      #phis[3] = phis[3] + solnElemLocal[r6] * bb3;

      # Bs - shear pat
      Bmat[1,r3] = bb1;    Bmat[1,r4] = -bb3;
      Bmat[2,r3] = bb2;    Bmat[2,r5] = -bb3;
    end

    # Bs - shear part
    strains[1] = -phis[1] + grad[3,1];
    strains[2] = -phis[2] + grad[3,2];

    sigf = Ds*strains;

    Klocal = Klocal + (transpose(Bmat)*Ds)*(dvol*Bmat);
    Flocal = Flocal - (transpose(Bmat)*(dvol*sigf));
end #gp

#fact = abs(maximum(Diagonal(Klocal)));
#fact = 10000.0;
#fact = E*h*h*h*100;
#println("fact = ", fact);
#for ii=1:nlbf
#  r1 = ndof*(ii-1)+6;
#
#  Klocal[r1, r1] += fact;
#
#  Flocal[r1]     += fact*(0.0-solnElemLocal[r1]);
#end


RotMatFull = zeros(nsize, nsize);

#for i=1:nlbf
#    # displacements
#    ind1 = ndof*(i-1)+1;
#    ind2 = ind1 + 2;
#    RotMatFull[ind1:ind2, ind1:ind2] = RotMat;
#    
#    # rotations
#    ind1 = ind1 + 3;
#    ind2 = ind1 + 2;
#    RotMatFull[ind1:ind2, ind1:ind2] = RotMat;
#end

for i=1:2*nlbf
  ind1 = 3*(i-1)+1;
  ind2 = ind1 + 2;

  RotMatFull[ind1:ind2, ind1:ind2] = RotMat;
end


#if (elnum in [1, 10000])
#  display(Klocal);
#end

# Transform stiffness matrix and force vector to global coordinate system
#
Klocal = RotMatFull'*Klocal*RotMatFull;
Flocal = RotMatFull'*Flocal;

#if (elnum in [1, 10000])
#  display(Klocal);
#end

# Global body force, e.g. gravity

if(nlbf == 6)
  elShape = "TRIA";
  degree  = 2;
  nGP = 3;
  gps1, gps2, gws = getGaussPointsTriangle(nGP);
elseif (nlbf == 4)
  elShape = "QUAD";
  degree  = 1;
  nGP = 4;
  gps1, gps2, gws = getGaussPointsQuad(nGP);
elseif (nlbf == 9)
  elShape = "QUAD";
  degree  = 2;
  nGP = 9;
  gps1, gps2, gws = getGaussPointsQuad(nGP);
end


bforce = -90.0;

param =[0.0, 0.0];

for gp=1:nGP
    param[1] = gps1[gp];
    param[2] = gps2[gp];

    N, dN_dx, dN_dy, Jac = shape_functions_Derivatives_2D(elShape, degree, nlbf, param[1], param[2], xNode, yNode);

    dvol = gws[gp] * Jac;

    for ii=1:nlbf
      bb3 = dvol*N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;

      #Flocal[r1  ] = Flocal[r1  ] + (bb3*bforce) ;
      #Flocal[r2  ] = Flocal[r2  ] + (bb3*bforce) ;
      Flocal[r3  ] = Flocal[r3  ] + (bb3*bforce) ;
    end
end #gp

return Klocal, Flocal
end

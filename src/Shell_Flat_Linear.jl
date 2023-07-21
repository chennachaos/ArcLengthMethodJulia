using LinearAlgebra

function Shell_Flat_Linear(elmDat, nodeNums, e, XX, soln, bf)

af = 1.0;

nlbf   = size(nodeNums)[1];
degree = 1;
if (nlbf == 9)
    degree = 2;
end
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

xNode = zeros(nlbf,1);
yNode = zeros(nlbf,1);
zNode = zeros(nlbf,1);

solnElem = zeros(nsize,1);

for ii=1:nlbf
  xNode[ii] = XX[nodeNums[ii],1];
  yNode[ii] = XX[nodeNums[ii],2];
  zNode[ii] = XX[nodeNums[ii],3];

  r1 = ndof*(ii-1)+1;
  r2 = r1+1;
  r3 = r1+2;
  r4 = r1+3;
  r5 = r1+4;
  r6 = r1+5;

  jj = ndof*(nodeNums[ii]-1);

  solnElem[r1] = soln[jj+1];
  solnElem[r2] = soln[jj+2];
  solnElem[r3] = soln[jj+3];
  solnElem[r4] = soln[jj+4];
  solnElem[r5] = soln[jj+5];
  solnElem[r6] = soln[jj+6];
end


#printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
#printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

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

      grad[1,1] = grad[1,1] + solnElem[r1] * bb1;
      grad[1,2] = grad[1,2] + solnElem[r1] * bb2;
      grad[2,1] = grad[2,1] + solnElem[r2] * bb1;
      grad[2,2] = grad[2,2] + solnElem[r2] * bb2;
      grad[3,1] = grad[3,1] + solnElem[r3] * bb1;
      grad[3,2] = grad[3,2] + solnElem[r3] * bb2;

      phis[1] = phis[1] + solnElem[r4] * bb3;
      phis[2] = phis[2] + solnElem[r5] * bb3;
      phis[3] = phis[3] + solnElem[r6] * bb3;

      dphis[1,1] = dphis[1,1] + solnElem[r4] * bb1;
      dphis[1,2] = dphis[1,2] + solnElem[r4] * bb2;

      dphis[2,1] = dphis[2,1] + solnElem[r5] * bb1;
      dphis[2,2] = dphis[2,2] + solnElem[r5] * bb2;

      dphis[3,1] = dphis[3,1] + solnElem[r6] * bb1;
      dphis[3,2] = dphis[3,2] + solnElem[r6] * bb2;
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


    for ii=1:nlbf
      bb3 = dvol*N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;

      Flocal[r3  ] = Flocal[r3  ] + (bb3*bforce) ;
    end
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

      grad[1,1] = grad[1,1] + solnElem[r1] * bb1;
      grad[1,2] = grad[1,2] + solnElem[r1] * bb2;
      grad[2,1] = grad[2,1] + solnElem[r2] * bb1;
      grad[2,2] = grad[2,2] + solnElem[r2] * bb2;
      grad[3,1] = grad[3,1] + solnElem[r3] * bb1;
      grad[3,2] = grad[3,2] + solnElem[r3] * bb2;

      phis[1] = phis[1] + solnElem[r4] * bb3;
      phis[2] = phis[2] + solnElem[r5] * bb3;
      phis[3] = phis[3] + solnElem[r6] * bb3;
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


fact = abs(maximum(Diagonal(Klocal)));
#println("fact = ", fact);
for ii=1:nlbf
  r1 = ndof*(ii-1)+6;

  Klocal[r1, r1] += fact;

  Flocal[r1]     += fact*(0.0-solnElem[r1]);
end


return Klocal, Flocal
end
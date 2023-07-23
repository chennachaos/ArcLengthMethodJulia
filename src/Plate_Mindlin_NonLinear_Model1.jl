function Plate_Mindlin_NonLinear_Model1(elmDat, nodeNums, e, XX, soln, bf)

af = 1.0;

nlbf   = size(nodeNums)[1];
degree = 1;
if (nlbf == 9)
    degree = 2;
end
ndof   = 5;
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

Dm = zeros(3,3);
Dm[1,1] = 1.0;   Dm[1,2] = nu;
Dm[2,1] = nu;    Dm[2,2] = 1.0;
Dm[3,3] = (1.0-nu)*0.5;

Dm = (E*h/(1.0-nu*nu))*Dm;

Db = zeros(3,3);
Db = (h*h/12.0)*Dm;



Ds = zeros(2,2);
Ds[1,1] = G;
Ds[2,2] = G;
Ds = (h*kappa)*Ds;

Dmat = zeros(8,8);

Dmat[1:3,1:3] = Dm;
Dmat[4:5,4:5] = Ds;
Dmat[6:8,6:8] = Db;



#printf(" %14.12f \t %14.12f \t %14.12f \n ", Dm[0][0], Dm[0][1], Dm[0][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n ", Dm[1][0], Dm[1][1], Dm[1][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n\n\n ", Dm[2][0], Dm[2][1], Dm[2][2]);
#cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;

xNode = zeros(nlbf,1);
yNode = zeros(nlbf,1);

solnElem = zeros(nsize,1);

for ii=1:nlbf
  xNode[ii] = XX[nodeNums[ii],1];
  yNode[ii] = XX[nodeNums[ii],2];

  r1 = ndof*(ii-1)+1;
  r2 = r1+1;
  r3 = r1+2;
  r4 = r1+3;
  r5 = r1+4;

  jj = ndof*(nodeNums[ii]-1);

  solnElem[r1] = soln[jj+1];
  solnElem[r2] = soln[jj+2];
  solnElem[r3] = soln[jj+3];
  solnElem[r4] = soln[jj+4];
  solnElem[r5] = soln[jj+5];
end

#printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
#printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

Klocal=zeros(nsize,nsize); # Local stiffness matrix
Flocal=zeros(nsize,1);   # Local load vector


############################################
#
#   part 1. -- contribution due to bending

nGP = 2;
if (degree == 2)
  nGP = 3;
end

gpvec, gwvec = get_Gauss_points(nGP);

param =[0.0, 0.0];

for gp2=1:nGP
    JacMult = gwvec[gp2];

    param[2] = gpvec[gp2];

for gp1=1:nGP
    param[1] = gpvec[gp1];

    N, dN_dx, dN_dy, Jac = shape_functions_Lagrange_2D(nodeNums, XX, degree, param[1], param[2]);

    dvol = gwvec[gp1] * JacMult * Jac;

    grad = zeros(3,3);
    phis = [0.0 0.0];
    dphis = zeros(2,2);
    strains = zeros(8,1);
    
    for ii=1:nlbf
      bb1 = dN_dx[ii];
      bb2 = dN_dy[ii];
      bb3 = N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;
      r4 = r1+3;
      r5 = r1+4;

      grad[1,1] = grad[1,1] + solnElem[r1] * bb1;
      grad[1,2] = grad[1,2] + solnElem[r1] * bb2;
      grad[2,1] = grad[2,1] + solnElem[r2] * bb1;
      grad[2,2] = grad[2,2] + solnElem[r2] * bb2;
      grad[3,1] = grad[3,1] + solnElem[r3] * bb1;
      grad[3,2] = grad[3,2] + solnElem[r3] * bb2;

      phis[1] = phis[1] + solnElem[r4] * bb3;
      phis[2] = phis[2] + solnElem[r5] * bb3;

      dphis[1,1] = dphis[1,1] + solnElem[r4] * bb1;
      dphis[1,2] = dphis[1,2] + solnElem[r4] * bb2;

      dphis[2,1] = dphis[2,1] + solnElem[r5] * bb1;
      dphis[2,2] = dphis[2,2] + solnElem[r5] * bb2;
    end

    Bmat = zeros(8,ndof*nlbf);
    Ga   = zeros(2,ndof*nlbf);
    for ii=1:nlbf
        bb1 = dN_dx[ii];
        bb2 = dN_dy[ii];
        bb3 = N[ii];

        r1 = ndof*(ii-1)+1;
        r2 = r1+1;
        r3 = r1+2;
        r4 = r1+3;
        r5 = r1+4;

        # Bp
        Bmat[1,r1] = bb1;    Bmat[1,r3] = grad[3,1]*bb1;
        Bmat[2,r2] = bb2;    Bmat[2,r3] = grad[3,2]*bb2;
        Bmat[3,r1] = bb2;    Bmat[3,r2] = bb1;           Bmat[3,r3] = grad[3,1]*bb2+grad[3,2]*bb1;

        # Bs
        #Bmat[4,r3] = bb1;    Bmat[4,r4] = bb3;
        #Bmat[5,r3] = bb2;    Bmat[5,r5] = bb3;

        # Bb
        Bmat[6,r4] = bb1;
        Bmat[7,r5] = bb2;
        Bmat[8,r4] = bb2;    Bmat[8,r5] = bb1;


        Ga[1,r3] = bb1;
        Ga[2,r3] = bb2;
    end


    strains[1] = grad[1,1] + 0.5*grad[3,1]*grad[3,1];
    strains[2] = grad[2,2] + 0.5*grad[3,2]*grad[3,2];
    strains[3] = grad[1,2] + grad[2,1] + grad[3,1]*grad[3,2];

    #strains[4] = phis[1] + grad[3,1];
    #strains[5] = phis[2] + grad[3,2];

    strains[6] = dphis[1,1];
    strains[7] = dphis[2,2];
    strains[8] = dphis[1,2] + dphis[2,1];

    sigf = Dmat*strains;

    Tmat = zeros(2,2);
    Tmat[1,1] = sigf[1];    Tmat[1,2] = sigf[3];
    Tmat[2,1] = sigf[3];    Tmat[2,2] = sigf[2];

    Klocal = Klocal + (transpose(Bmat)*Dmat)*(dvol*Bmat);
    Klocal = Klocal + (transpose(Ga)*Tmat)*(dvol*Ga);
    Flocal = Flocal - (transpose(Bmat)*(dvol*sigf));


    for ii=1:nlbf
      bb3 = dvol*N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;

      Flocal[r3  ] = Flocal[r3  ] + (bb3*bforce) ;
    end
end #gp1
end #gp2

#printMatrix(Klocal);  printf("\n\n\n"); printVector(Flocal);


############################################
#
#   part 2. -- contribution due to shear

nGP = 1;
if (degree == 2)
  nGP = 2;
end

gpvec, gwvec = get_Gauss_points(nGP);


for gp2=1:nGP
  JacMult = gwvec[gp2];

  param[2] = gpvec[gp2];

for gp1=1:nGP
  param[1] = gpvec[gp1];

  N, dN_dx, dN_dy, Jac = shape_functions_Lagrange_2D(nodeNums, XX, degree, param[1], param[2]);

  dvol = gwvec[gp1] * JacMult * Jac;

  grad = zeros(3,3);
  phis = [0.0 0.0];
  dphis = zeros(2,2);
  strains = zeros(8,1);
  
  for ii=1:nlbf
    bb1 = dN_dx[ii];
    bb2 = dN_dy[ii];
    bb3 = N[ii];

    r1 = ndof*(ii-1)+1;
    r2 = r1+1;
    r3 = r1+2;
    r4 = r1+3;
    r5 = r1+4;

    grad[1,1] = grad[1,1] + solnElem[r1] * bb1;
    grad[1,2] = grad[1,2] + solnElem[r1] * bb2;
    grad[2,1] = grad[2,1] + solnElem[r2] * bb1;
    grad[2,2] = grad[2,2] + solnElem[r2] * bb2;
    grad[3,1] = grad[3,1] + solnElem[r3] * bb1;
    grad[3,2] = grad[3,2] + solnElem[r3] * bb2;

    phis[1] = phis[1] + solnElem[r4] * bb3;
    phis[2] = phis[2] + solnElem[r5] * bb3;

    dphis[1,1] = dphis[1,1] + solnElem[r4] * bb1;
    dphis[1,2] = dphis[1,2] + solnElem[r4] * bb2;

    dphis[2,1] = dphis[2,1] + solnElem[r5] * bb1;
    dphis[2,2] = dphis[2,2] + solnElem[r5] * bb2;
  end

  Bmat = zeros(8,ndof*nlbf);
  Ga   = zeros(2,ndof*nlbf);
  for ii=1:nlbf
      bb1 = dN_dx[ii];
      bb2 = dN_dy[ii];
      bb3 = N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;
      r4 = r1+3;
      r5 = r1+4;

      # Bp
      #Bmat[1,r1] = bb1;    Bmat[1,r3] = grad[3,1]*bb1;
      #Bmat[2,r2] = bb2;    Bmat[2,r3] = grad[3,2]*bb2;
      #Bmat[3,r1] = bb2;    Bmat[3,r2] = bb1;           Bmat[3,r3] = grad[3,1]*bb2+grad[3,2]*bb1;

      # Bs
      Bmat[4,r3] = bb1;    Bmat[4,r4] = bb3;
      Bmat[5,r3] = bb2;    Bmat[5,r5] = bb3;

      # Bb
      #Bmat[6,r4] = bb1;
      #Bmat[7,r5] = bb2;
      #Bmat[8,r4] = bb2;    Bmat[8,r5] = bb1;


      #Ga[1,r3] = bb1;
      #Ga[2,r3] = bb2;
  end


  #strains[1] = grad[1,1] + 0.5*grad[3,1]*grad[3,1];
  #strains[2] = grad[2,2] + 0.5*grad[3,2]*grad[3,2];
  #strains[3] = grad[1,2] + grad[2,1] + grad[3,1]*grad[3,2];

  strains[4] = phis[1] + grad[3,1];
  strains[5] = phis[2] + grad[3,2];

  #strains[6] = dphis[1,1];
  #strains[7] = dphis[2,2];
  #strains[8] = dphis[1,2] + dphis[2,1];

  sigf = Dmat*strains;

  Tmat = zeros(2,2);
  Tmat[1,1] = sigf[1];    Tmat[1,2] = sigf[3];
  Tmat[2,1] = sigf[3];    Tmat[2,2] = sigf[2];

  Klocal = Klocal + (transpose(Bmat)*Dmat)*(dvol*Bmat);
  #Klocal = Klocal + (transpose(Ga)*Tmat)*(dvol*Ga);
  Flocal = Flocal - (transpose(Bmat)*(dvol*sigf));

end #gp1
end #gp2

return Klocal, Flocal
end
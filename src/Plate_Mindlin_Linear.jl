function Plate_Mindlin_Linear(elmDat, nodeNums, XX, soln, bf)

af = 1.0;

nlbf   = size(nodeNums)[1];
ndof   = 3;
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

Db = zeros(3,3);
Db[1,1] = 1.0;   Db[1,2] = nu;
Db[2,1] = nu;    Db[2,2] = 1.0;
Db[3,3] = (1.0-nu)*0.5;

Db = (I*E/(1.0-nu*nu))*Db;

Ds = zeros(2,2);
Ds[1,1] = G;
Ds[2,2] = G;
Ds = (h*kappa)*Ds;

Bb = zeros(3,3*nlbf);
Bs = zeros(2,3*nlbf);

#printf(" %14.12f \t %14.12f \t %14.12f \n ", Db[0][0], Db[0][1], Db[0][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n ", Db[1][0], Db[1][1], Db[1][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n\n\n ", Db[2][0], Db[2][1], Db[2][2]);
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

  jj = ndof*(nodeNums[ii]-1);

  solnElem[r1] = soln[jj+1];
  solnElem[r2] = soln[jj+2];
  solnElem[r3] = soln[jj+3];
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

for gp=1:nGP
    param[1] = gps1[gp];
    param[2] = gps2[gp];

    N, dN_dx, dN_dy, Jac = shape_functions_Derivatives_2D(elShape, degree, nlbf, param[1], param[2], xNode, yNode);

    dvol = gws[gp] * Jac;

    for ii=1:nlbf
        bb1 = dN_dx[ii];
        bb2 = dN_dy[ii];
        bb3 = N[ii];

        r1 = 3*(ii-1)+1;
        r2 = r1+1;
        r3 = r1+2;

        Bb[1,r2] = bb1;
        Bb[2,r3] = bb2;
        Bb[3,r2] = bb2;    Bb[3,r3] = bb1;  
    end

    sigf = Db*(Bb*solnElem);

    Klocal = Klocal + (transpose(Bb)*Db)*(dvol*Bb);
    Flocal = Flocal - (transpose(Bb)*(dvol*sigf));


    for ii=1:nlbf
      bb3 = dvol*N[ii];

      r1 = ndof*(ii-1)+1;
      r2 = r1+1;
      r3 = r1+2;

      Flocal[r1  ] = Flocal[r1  ] + (bb3*bforce) ;
      #Flocal[r2] = Flocal[r2] + (bb3*bforce) ;
      #Flocal[r3] = Flocal[r3] + (bb3*bforce) ;
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


for gp=1:nGP
    param[1] = gps1[gp];
    param[2] = gps2[gp];

    N, dN_dx, dN_dy, Jac = shape_functions_Derivatives_2D(elShape, degree, nlbf, param[1], param[2], xNode, yNode);

    dvol = gws[gp] * Jac;

    for ii=1:nlbf
        bb1 = dN_dx[ii];
        bb2 = dN_dy[ii];
        bb3 = N[ii];

        r1 = 3*(ii-1)+1;
        r2 = r1+1;
        r3 = r1+2;

        Bs[1,r1] = bb1;    Bs[1,r2] = bb3;
        Bs[2,r1] = bb2;    Bs[2,r3] = bb3;
    end

    sigc = Ds*(Bs*solnElem);

    Klocal = Klocal + (transpose(Bs)*Ds)*(dvol*Bs);
    Flocal = Flocal - (transpose(Bs)*(dvol*sigc));
end #gp

return Klocal, Flocal
end
function Mindlinplate_Linear(elmDat, nodeNums, e, XX, soln, bf)

af = 1.0;

nlbf   = size(nodeNums)[1];
degree = 1;
if (nlbf == 9)
    degree = 2;
end
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

Df = zeros(3,3);
Df[1,1] = 1.0;   Df[1,2] = nu;
Df[2,1] = nu;    Df[2,2] = 1.0;
Df[3,3] = (1.0-nu)*0.5;

Df = (I*E/(1.0-nu*nu))*Df;

Dc = zeros(2,2);
Dc[1,1] = G;
Dc[2,2] = G;
Dc = (h*kappa)*Dc;

Bf = zeros(3,3*nlbf);
Bc = zeros(2,3*nlbf);

#printf(" %14.12f \t %14.12f \t %14.12f \n ", Df[0][0], Df[0][1], Df[0][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n ", Df[1][0], Df[1][1], Df[1][2]);
#printf(" %14.12f \t %14.12f \t %14.12f \n\n\n ", Df[2][0], Df[2][1], Df[2][2]);
#cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;

xNode = zeros(nlbf,1);
yNode = zeros(nlbf,1);

solnElem = zeros(nsize,1);

for ii=1:nlbf
  xNode[ii] = XX[nodeNums[ii],1];
  yNode[ii] = XX[nodeNums[ii],2];

  TI   = 3*(ii-1)+1;
  TIp1 = TI+1;
  TIp2 = TI+2;

  jj = ndof*(nodeNums[ii]-1);

  solnElem[TI]   = soln[jj+1];
  solnElem[TIp1] = soln[jj+2];
  solnElem[TIp2] = soln[jj+3];
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

    for ii=1:nlbf
        bb1 = dN_dx[ii];
        bb2 = dN_dy[ii];
        bb3 = N[ii];

        TI   = 3*(ii-1)+1;
        TIp1 = TI+1;
        TIp2 = TI+2;

        Bf[1,TIp1] = bb1;
        Bf[2,TIp2] = bb2;
        Bf[3,TIp1] = bb2;    Bf[3,TIp2] = bb1;  
    end

    sigf = Df*(Bf*solnElem);

    Klocal = Klocal + (transpose(Bf)*Df)*(dvol*Bf);
    Flocal = Flocal - (transpose(Bf)*(dvol*sigf));


    for ii=1:nlbf
      bb3 = dvol*N[ii];

      TI   = 3*(ii-1)+1;
      TIp1 = TI+1;
      TIp2 = TI+2;

      Flocal[TI  ] = Flocal[TI  ] + (bb3*bforce) ;
      #Flocal[TIp1] = Flocal[TIp1] + (bb3*bforce) ;
      #Flocal[TIp2] = Flocal[TIp2] + (bb3*bforce) ;
    end
end #gp1
end #gp2

#printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);


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

    for ii=1:nlbf
        bb1 = dN_dx[ii];
        bb2 = dN_dy[ii];
        bb3 = N[ii];

        TI   = 3*(ii-1)+1;
        TIp1 = TI+1;
        TIp2 = TI+2;

        Bc[1,TI] = bb1;    Bc[1,TIp1] = bb3;
        Bc[2,TI] = bb2;    Bc[2,TIp2] = bb3;
    end

    sigc = Dc*(Bc*solnElem);

    Klocal = Klocal + (transpose(Bc)*Dc)*(dvol*Bc);
    Flocal = Flocal - (transpose(Bc)*(dvol*sigc));
end #gp1
end #gp2

return Klocal, Flocal
end
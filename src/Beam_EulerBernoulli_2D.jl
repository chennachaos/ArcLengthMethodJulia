function Beam_EulerBernoulli_2D(elmDat, nodeNums, nodecoords, soln, bf)

nlocal = 4;
npElem = 2;
ndof   = 3;
nsize  = 6;

rho  = elmDat[2];
A    = elmDat[3];
I    = elmDat[4];
E    = elmDat[5];
nu   = elmDat[6];
kappa= elmDat[7];

G  = E/2.0/(1.0+nu);
EA = E*A;
EI = E*I;
GA = G*A*kappa;


Klocal=zeros(nsize,nsize); # Local stiffness matrix
Flocal=zeros(nsize,1);   # Local load vector

node1 = nodeNums[1]
node2 = nodeNums[2]

x0 = [0.0, 0.0];
y0 = [0.0, 0.0];

x0[1] = nodecoords[node1,1];
y0[1] = nodecoords[node1,2];
x0[2] = nodecoords[node2,1];
y0[2] = nodecoords[node2,2];

dx = x0[2] - x0[1];
dy = y0[2] - y0[1];
h  = sqrt(dx*dx+dy*dy);

cth0 = dx/h;
sth0 = dy/h;
    
RotMat=zeros(6,6);

RotMat[1,1] =  cth0; RotMat[1,2] = -sth0;
RotMat[2,1] =  sth0; RotMat[2,2] =  cth0;
RotMat[3,3] =  1.0;
RotMat[4,4] =  cth0; RotMat[4,5] = -sth0;
RotMat[5,4] =  sth0; RotMat[5,5] =  cth0;
RotMat[6,6] =  1.0;


uxn=[0.0,0.0];
uzn=[0.0,0.0];
btn=[0.0,0.0];


dispC = zeros(nsize,1);

dispC[1] = soln[ndof*(node1-1)+1];
dispC[2] = soln[ndof*(node1-1)+2];
dispC[3] = soln[ndof*(node1-1)+3];

dispC[4] = soln[ndof*(node2-1)+1];
dispC[5] = soln[ndof*(node2-1)+2];
dispC[6] = soln[ndof*(node2-1)+3];

#dummy = RotMat'*[uxn[1]; uzn[1]; btn[1]; uxn[2]; uzn[2]; btn[2]];

fact = EA/h;

Klocal[1,1] =  fact;
Klocal[1,4] = -fact;
Klocal[4,1] = -fact;
Klocal[4,4] =  fact;

fact = 12.0*EI/h/h/h;

Klocal[2,2] =  fact;
Klocal[2,5] = -fact;
Klocal[5,2] = -fact;
Klocal[5,5] =  fact;

fact = 6.0*EI/h/h;

Klocal[2,3] =  fact;
Klocal[3,2] =  fact;
Klocal[2,6] =  fact;
Klocal[6,2] =  fact;

Klocal[3,5] = -fact;
Klocal[5,3] = -fact;
Klocal[5,6] = -fact;
Klocal[6,5] = -fact;

fact = EI/h;

Klocal[3,3] =  4.0*fact;
Klocal[3,6] =  2.0*fact;
Klocal[6,3] =  2.0*fact;
Klocal[6,6] =  4.0*fact;


Flocal[1] = Flocal[1] + 0.5*h*bf[1];
Flocal[4] = Flocal[4] + 0.5*h*bf[1];

Flocal[2] = Flocal[2] + 0.5*h*bf[2];
Flocal[5] = Flocal[5] + 0.5*h*bf[2];

Flocal[3] = Flocal[3] + h*h*bf[2]/12.0;
Flocal[6] = Flocal[6] - h*h*bf[2]/12.0;

#Flocal = RotMat*Flocal;
#Klocal = RotMat*Klocal*RotMat';

Flocal  = Flocal - Klocal*dispC;

return Klocal, Flocal
end
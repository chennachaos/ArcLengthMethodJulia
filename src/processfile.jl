function  processfile(fname)
# nodecoords: global coordinates of the nodes, x, y, and z
# elemConn: element connectivities
# nelem: total number of elements
# nnode: total number of nodes
# nperelem: number of nodes per element for all elements
# ndof: number of degrees of per node

fid=open(fname,"r");

# ndim

line    = readline(fid);
println("line = ", line)
linestr = split(line, ",");
ndim    = parse(Int32, linestr[2])
println("ndim = ", ndim)

# ndof

line    = readline(fid);
println("line = ", line)
linestr = split(line, ",");
ndof    = parse(Int32, linestr[2])
println("ndof = ", ndof)

# nodes

line    = readline(fid);
println("line = ", line)
linestr = split(line, ",");
nnode   = parse(Int32, linestr[2])


nodecoords = zeros(nnode,ndim);
for i=1:nnode
    line = readline(fid);
    linestr = split(line, ",");

    nodecoords[i,1] = parse(Float64, linestr[2]);
    nodecoords[i,2] = parse(Float64, linestr[3]);
    if(ndim == 3)
      nodecoords[i,3] = parse(Float64, linestr[4]);
    end
end    


# element data

line      = readline(fid);
linestr   = split(line, ",");
nelemData = parse(Int32, linestr[2])

elemData = zeros(nelemData, 10);
for i=1:nelemData
    line = readline(fid);
    linestr = split(line, ",");

    for j=1:10
      elemData[i,j] = parse(Float64, linestr[j+1]);
    end
end

# elements

line    = readline(fid);
linestr = split(line, ",");
nelem   = parse(Int32, linestr[2])
npElem   = parse(Int32, linestr[3])

nsize = npElem*ndof;
neq = nnode*ndof;

elemConn = zeros(Int32, nelem, npElem+2);
for i=1:nelem
    line    = readline(fid);
    linestr = split(line, ",");

    elemConn[i,1] = parse(Int32, linestr[2]);
    elemConn[i,2] = parse(Int32, linestr[3]);

    for j=1:npElem
      elemConn[i,2+j] = parse(Int32, linestr[3+j]);
    end
end

# Dirichlet boundary conditions

line    = readline(fid);
linestr = split(line, ",");
nDBC    = parse(Int32, linestr[2])

#dbclist = zeros(nDBC, 3);
dbcnodes = zeros(Int64, nDBC, 1);
for i=1:nDBC
    line    = readline(fid);
    linestr = split(line, ",");

    n1 = parse(Int32, linestr[1]);
    n2 = parse(Int32, linestr[2]);
#    dbclist(i,3) = double(str2num(linestr{1,3}));
    dbcnodes[i] = (n1-1)*ndof+n2;
end

#display("text/plain", dbcnodes)
assy4r = setdiff(collect(1:neq), dbcnodes);

# Force boundary conditions

line    = readline(fid);
linestr = split(line, ",");
nFBC    = parse(Int32, linestr[2])
fbclist = zeros(nFBC, 3);

dof_force = zeros(Int32, nFBC, 1);
Fext = zeros(neq, 1);
for i=1:nFBC
    line = readline(fid);
    linestr = split(line, ",");

    n1 = parse(Int32, linestr[1]);
    n2 = parse(Int32, linestr[2]);
    ind = (n1-1)*ndof + n2;

    dof_force[i] = ind;
    Fext[ind] = parse(Float64, linestr[3]);
end

# for output

line    = readline(fid);
linestr = split(line, ",");
nOutput = parse(Int32, linestr[2])
outputlist = zeros(Int32, nOutput, 1);

for i=1:nOutput
    line    = readline(fid);
    linestr = split(line, ",");

    n1 = parse(Int32, linestr[1]);
    n2 = parse(Int32, linestr[2]);

    outputlist[i,1] = (n1-1)*ndof+n2;
end


# Solver parameters

line    = readline(fid);
linestr = split(line, ",");

line    = readline(fid)
linestr = split(line, ",");
maxloadSteps = parse(Int32, linestr[1]);

line     = readline(fid)
linestr  = split(line, ",");
loadincr = parse(Float64, linestr[1]);

close(fid);

# data structures

LM  = zeros(Int32, nelem, nsize);

for e=1:nelem
    count = 1;
    for jj=1:npElem
        ind = ndof*(elemConn[e,jj+2]-1)+1;
        for kk=1:ndof
            LM[e,count] = ind;
            ind = ind + 1;
            count = count + 1;
        end
    end
    count = count - 1;
end


return  ndim, ndof, nnode, nelem, nodecoords, elemConn, elemData, LM, neq, assy4r, dof_force, Fext, maxloadSteps, loadincr, outputlist
end
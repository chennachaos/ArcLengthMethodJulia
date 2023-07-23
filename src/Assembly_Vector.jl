function Assembly_Vector(global_vector,local_vector,forAssyVec)
nen=size(forAssyVec,1);
for aa=1:nen
    mm=forAssyVec[aa];
    if mm!=0
        global_vector[mm,1]=global_vector[mm,1]+local_vector[aa,1];
    end
end

return global_vector
end
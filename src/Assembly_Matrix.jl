function Assembly_Matrix(global_matrix,local_matrix,forAssyVec)
# Assemble global_matrix
nen=size(forAssyVec,1);
for aa=1:nen
    mm=forAssyVec[aa];
    if mm!=0
        for bb=1:nen
            nn=forAssyVec[bb];
            if nn!=0
                global_matrix[mm,nn]=global_matrix[mm,nn]+local_matrix[aa,bb];
            end
        end
    end
end

return global_matrix
end
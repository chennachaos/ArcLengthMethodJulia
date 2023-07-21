

#
# ylocal
# ^
# |
# 4------3
# |      |
# |      |
# 1------2--->xlocal
#
#

function  compute_transformation_matrix(coordsGlobal, solnElemGlobal)

    # base vectors in global coordinate system
    EX = [1.0, 0.0, 0.0];
    EY = [0.0, 1.0, 0.0];
    EZ = [0.0, 0.0, 1.0];

    # vectors in local coordinate system
    # line 1-2
    v21 = [coordsGlobal[1,2]-coordsGlobal[1,1],  coordsGlobal[2,2]-coordsGlobal[2,1],   coordsGlobal[3,2]-coordsGlobal[3,1]];
    # line 1-4
    v41 = [coordsGlobal[1,4]-coordsGlobal[1,1],  coordsGlobal[2,4]-coordsGlobal[2,1],   coordsGlobal[3,4]-coordsGlobal[3,1]];
    
    v21 = v21/norm(v21,2);
    v41 = v41/norm(v41,2);
    
    exl = v21;
    ezl = cross(exl, v41);
    eyl = cross(ezl, exl);

    # transformation matrix from global to local
    RotMat = zeros(3,3);
    RotMat[1,1] = dot(exl, EX);    RotMat[1,2] = dot(exl, EY);    RotMat[1,3] = dot(exl, EZ);
    RotMat[2,1] = dot(eyl, EX);    RotMat[2,2] = dot(eyl, EY);    RotMat[2,3] = dot(eyl, EZ);
    RotMat[3,1] = dot(ezl, EX);    RotMat[3,2] = dot(ezl, EY);    RotMat[3,3] = dot(ezl, EZ);


    ## midpoing of line 1-2
    #Xi = 0.5*[coordsGlobal[1,1]+coordsGlobal[1,2],  coordsGlobal[2,1]+coordsGlobal[2,2],   coordsGlobal[3,1]+coordsGlobal[3,2]];
    ## midpoing of line 2-3
    #Xj = 0.5*[coordsGlobal[1,2]+coordsGlobal[1,3],  coordsGlobal[2,2]+coordsGlobal[2,3],   coordsGlobal[3,2]+coordsGlobal[3,3]];
    ## midpoing of line 3-4
    #Xk = 0.5*[coordsGlobal[1,3]+coordsGlobal[1,4],  coordsGlobal[2,3]+coordsGlobal[2,4],   coordsGlobal[3,3]+coordsGlobal[3,4]];
    ## midpoing of line 4-1
    #Xl = 0.5*[coordsGlobal[1,4]+coordsGlobal[1,1],  coordsGlobal[2,4]+coordsGlobal[2,1],   coordsGlobal[3,4]+coordsGlobal[3,1]];
    #
    ## line l-j
    #Vlj = [Xj[1]-Xl[1], Xj[2]-Xl[2], Xj[3]-Xl[3]];
    ## line i-k
    #Vik = [Xk[1]-Xi[1], Xk[2]-Xi[2], Xk[3]-Xi[3]];
    #
    #Vlj = Vlj/norm(Vlj,2);
    #Vik = Vik/norm(Vik,2);
    #
    #exl = Vlj;
    #ezl = cross(exl, Vik);
    #eyl = cross(ezl, exl);
    #
    ## Transformation matrix
    #RotMat = zeros(3,3);
    #
    #RotMat = [exl eyl ezl]';

    IsFlat = false;
    if (RotMat[1,1] == 1.0 && RotMat[2,2] == 1.0 && RotMat[3,3] == 1.0)
        IsFlat = true;
    end
    #println("\n IsFlat = ", IsFlat);

    if (IsFlat)
        coordsLocal   = 1.0*coordsGlobal;
        solnElemLocal = 1.0*solnElemGlobal;
    else
        # transform nodal coordinates into local system
        coordsLocal = 0.0*coordsGlobal;
        nnode = size(coordsGlobal,2);

        for i=1:nnode
            coordsLocal[:,i] = RotMat*(coordsGlobal[:,i]-coordsGlobal[:,1]);
        end
    
        # transform nodal solution into local system
        solnElemLocal = 0.0*solnElemGlobal;
        for i=1:nnode
            ind1 = 6*(i-1)+1;
            ind2 = ind1 + 2;
        
            # displacements
            solnElemLocal[ind1:ind2] = RotMat*solnElemGlobal[ind1:ind2];

            # rotations
            ind1 = ind1 + 3;
            ind2 = ind1 + 2;
            solnElemLocal[ind1:ind2] = RotMat*solnElemGlobal[ind1:ind2];
        end
    end

    return RotMat, coordsLocal, solnElemLocal;
end
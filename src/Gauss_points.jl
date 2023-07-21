#module Quadrature

function getGaussPoints1D(nGP)

#######################
#######################
if(nGP == 1)
    gps=0.0;
    gws=2.0;
elseif(nGP == 2)
    ### 2 Point quadrature rule
    gps=[-0.577350269,0.577350269];
    gws=[1.0,1.0];
elseif(nGP == 3)
    ### 3 Point quadrature rule
    gps=[-sqrt(3/5),0,sqrt(3/5)];
    gws=[5/9,8/9,5/9];
elseif(nGP == 4)
    ### 4 Point quadrature rule
    gps = [-0.861136311594953, -0.339981043584856, 0.339981043584856, 0.861136311594953];
    gws = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454];
else
    ### 5 Point quadrature rule
    gps=[-0.9061798459,-0.5384693101,0,0.5384693101,0.9061798459];
    gws=[0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851];
end

return  gps, gws
end


function getGaussPointsTriangle(ngp)

    r1d3 = 1.0/3.0;
    fact = 0.5

    if(ngp == 1) # 1 Point quadrature rule
        fact = 0.5;

        gps1 = [r1d3];
        gps2 = [r1d3];
        gws  = [fact];
    elseif (ngp == 3) # 3 Point quadrature rule

        fact = 0.5*r1d3;

        gps1 = [1.0/6.0, 4.0/6.0, 1.0/6.0];
        gps2 = [1.0/6.0, 1.0/6.0, 4.0/6.0];
        gws  = [fact, fact, fact];
    elseif  (ngp == 6)  # 6 Point quadrature rule
        gps1 = [0.10810301816807022736, 0.44594849091596488632, 0.44594849091596488632, 0.81684757298045851308, 0.09157621350977074346, 0.09157621350977074346];
        gps2 = [0.44594849091596488632, 0.10810301816807022736, 0.44594849091596488632, 0.09157621350977074346, 0.81684757298045851308, 0.09157621350977074346];

        fact1 = 0.5*0.22338158967801146570;
        fact2 = 0.5*0.10995174365532186764;
        gws  = [fact1, fact1, fact1, fact2, fact2, fact2]
    end

    return gps1, gps2, gws;
end




function getGaussPointsQuad(ngp)

    if(ngp == 1) # 1 Point quadrature rule
        gps1 = [0.0];
        gps2 = [0.0];
        gws  = [4.0];
    elseif (ngp == 4) # 4 Point quadrature rule
        gps1 = [-0.577350269189626, 0.577350269189626, 0.577350269189626, -0.577350269189626];
        gps2 = [-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626];

        gws  = [1.0, 1.0, 1.0, 1.0];
    elseif  (ngp == 9)  # 9 Point quadrature rule
        gps1 = [-0.774596669241483,  0.774596669241483, 0.774596669241483, -0.774596669241483,  0.0,               0.774596669241483, 0.0,              -0.774596669241483, 0.0];
        gps2 = [-0.774596669241483, -0.774596669241483, 0.774596669241483,  0.774596669241483, -0.774596669241483, 0.0,               0.774596669241483, 0.0,               0.0];

        fact1 = 25.0/81.0;
        fact2 = 40.0/81.0;
        gws  = [fact1, fact1, fact1, fact1, fact2, fact2, fact2, fact2, 64.0/81.0]
    end

    return gps1, gps2, gws;
end


#end

using DelimitedFiles

function writeoutputvtk(filename,ndim,nElem,nNode,npElem,ndof,coords,elemNodeConn,soln)

    #Open the file and write to it
    fileID = open(filename,"w");

      # Directives
      println(fileID, "# vtk DataFile Version 4.0");
      println(fileID, "PoissonTwoD example");

      # ASCII or Binary
      println(fileID, "ASCII");

      # Type of dataset : Structured/Unstructured/Rectilinear/Polydata...
      println(fileID, "DATASET UNSTRUCTURED_GRID");
      
      # Coordinates of the points (nodes)
      println(fileID, "POINTS ", nNode, " float");

      if (ndim == 2)
        for ii=1:nNode
          writedlm(fileID, [coords[ii,1] coords[ii,2] 0.0]);
        end
      else
        for ii=1:nNode
          writedlm(fileID, [coords[ii,1]  coords[ii,2]  coords[ii,3]]);
        end
      end

      #WRITE(*,*) "aaaaaaaaaaaaaa"

      # Element<->Nodes connectivity
      # In VTK terminology, Cell<->Points
      # <number of nodes> <n1> <n2> <n3> ...
      # Starting index is 0
      ind = nElem*(npElem+1);
      println(fileID, "CELLS ", nElem, " ", ind);

      #if (ndim == 2)
        if (npElem == 2)
          # 2-noded line element
          n1 = 3;
          for ee=1:nElem
            writedlm(fileID, [npElem  elemNodeConn[ee,3]-1  elemNodeConn[ee,4]-1]);
          end
        elseif (npElem == 4)
          # Linear Quadrilateral element
          n1 = 9;
          for ee=1:nElem
            writedlm(fileID, [npElem  elemNodeConn[ee,3]-1  elemNodeConn[ee,4]-1  elemNodeConn[ee,5]-1  elemNodeConn[ee,6]-1]);
          end
        else #    VTK_BIQUADRATIC_QUAD = 28
          n1 = 28;
          for ee=1:nElem
            writedlm(fileID, [npElem  elemNodeConn[ee,3]-1  elemNodeConn[ee,4]-1  elemNodeConn[ee,5]-1  elemNodeConn[ee,6]-1 elemNodeConn[ee,7]-1 elemNodeConn[ee,8]-1 elemNodeConn[ee,9]-1 elemNodeConn[ee,10]-1 elemNodeConn[ee,11]-1]);
          end
        end
      #end

      # Cell type, as per VTK
      println(fileID, "CELL_TYPES ", nElem);
      for ee=1:nElem
        println(fileID, n1);
      end

      #WRITE(*,*) "bbbbbbbbbbbbbbbbb"
      # Cell data. Processor ID
      #println(fileID, "CELL_DATA ", nElem);
      #write(fileID, "SCALARS procid int 1\n");
      #write(fileID, "LOOKUP_TABLE default\n");
      #for ee=1:nElem
      #  println(fileID, 1);
      #end

      #WRITE(*,*) "ccccccccccccccccccc"

      # Point data
      println(fileID, "POINT_DATA ", nNode);
      if( (npElem == 2) && (ndof == 3) )
        #println(fileID, "SCALARS solution float 1");
        #println(fileID, "LOOKUP_TABLE default");

        #for ii=1:nNode
        #  println(fileID, soln[ndof*(ii-1)+1]);
        #end
        println(fileID, "VECTORS solution float");

        for ii=1:nNode
          writedlm(fileID, [soln[ndof*(ii-1)+1] soln[ndof*(ii-1)+2] 0.0]);
        end
      elseif (ndof == 3)
        #println(fileID, "SCALARS solution float 1");
        #println(fileID, "LOOKUP_TABLE default");

        #for ii=1:nNode
        #  println(fileID, soln[ndof*(ii-1)+1]);
        #end
        println(fileID, "VECTORS solution float");

        for ii=1:nNode
          writedlm(fileID, [0.0 0.0 soln[ndof*(ii-1)+1]]);
        end
      elseif (ndof == 5)
        println(fileID, "SCALARS solution float 1");
        println(fileID, "LOOKUP_TABLE default");

        for ii=1:nNode
          println(fileID, soln[ndof*(ii-1)+3]);
        end
      elseif (ndof == 6)
        println(fileID, "VECTORS solution float");

        for ii=1:nNode
          writedlm(fileID, [soln[ndof*(ii-1)+1] soln[ndof*(ii-1)+2] soln[ndof*(ii-1)+3]]);
        end
      end

      # close the file
      close(fileID);

    return;
end
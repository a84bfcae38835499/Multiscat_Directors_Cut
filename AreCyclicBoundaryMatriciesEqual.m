function [result] = AreCyclicBoundaryMatriciesEqual(A,B,Ndefect)
  if(~isequal(size(A),size(B)))
    disp(size(A))
    disp(size(B))
    error("Matricies not of equal size!")
  end
  Nsuper = length(A);
  %operations to check:
  % - rot90 (and 180,270)
  % - fliplr and flipud
  % - transpose
  %oh god wait some of these are the same. Like, a transpose is the same as
  %flipud followed by rot90. how can we generalise transformations and make
  %sure they're not repeated??
  %Hypothesis: Doing each rot, checking that, then checking its transpose,
  %then undoing the transpose and rotating again, gives the right result.
  %Then I just do that for each possible collumn/row shift!
  result = false;
  for mshift = 0:Nsuper-1
    for nshift = 0:Nsuper-1
      %Keep a constant and compare to a shifted version of B
      Bs = [B(:,2+mshift:end)  B(:,1:1+mshift)];
      Bs = [Bs(2+nshift:end,:);Bs(1:1+nshift,:)];
      indicies1 = ones(Ndefect,2,'int64')*137;
      indicies2 = ones(Ndefect,2,'int64')*137;
      counter1 = 1;
      counter2 = 1;
      for m = 0:int64(Nsuper-1)
        for n = 0:int64(Nsuper-1)
          if(A(m+1,n+1) == true)
            indicies1(counter1,1) = m;
            indicies1(counter1,2) = n;
            counter1 = counter1 + 1;
          end
          if(Bs(m+1,n+1) == true)
            indicies2(counter2,1) = m;
            indicies2(counter2,2) = n;
            counter2 = counter2 + 1;
          end
        end
      end
      i1 = sort(indicies1,1);
      i2 = sort(indicies2,1);
      %disp("i1 = ")
      %disp(i1)
      %disp("i2 = ")
      %disp(i2)
      if(isequal(i1,i2))
        disp("isequal!")
        result = true;
        return;
      end
    end
  end
end
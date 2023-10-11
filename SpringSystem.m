
[u, w, e] = spring(3, [1 1 1], [1 1 1], 1)
function [u,w,e] = spring(nmasses, constants, masses, fixed)
  g = 9.81;  

  %% One free, one fixed
  if fixed == 1
      %% Initialize the displacement matrix of size 1 x the number of masses
      u = zeros(1, nmasses);
      %% Initialize the elongation matrix of size 1 x the number of masses
      e = zeros(1, nmasses);
      %% Initiliaze the A matrix of size number of masses x number of masses
      A = zeros (nmasses, nmasses);
      %% Solving f = mg
      f = masses*g;
      for i = 1:nmasses
          %% Setting the e = u(x+1)-u(x) for the 1s
          A(i,i) = 1;
          %% Setting the e = u(x+1)-u(x) for the -1s
          if i > 1
              A(i, i-1) = -1;
          end
      end    
      %% Creating the C matrix of spring constants in the diagonal
      C = diag(constants);
      %% Creating the K matrix
      K = A'*C*A;
      %% Calling the SVD function for K  
      [U, S, Vt, c, I] = svdd(K);
      fprintf("Singular Values: \n")
      disp(diag(S))
      fprintf("\nEigenvalues: \n")
      disp((diag(S)).^2)
      fprintf("\nCondition Number:\n")
      disp(c)
      %% Solve for displacements, u = K inverse * f
      u = inv(K)*f';
      %% Back-calculate elongation (e = A*u) and internal stresses
      e = A*u;
      w = C*e;

  end

  %% Zero free, meaning both are fixed
  if fixed == 2
      u = zeros(1, nmasses);
      e = zeros(1, nmasses);
      A = zeros (nmasses+1, nmasses);
      f = masses*g;
      for i = 1:nmasses
          if i == 1
              A(i,i) = 1;
          else
              A(i, i-1) = -1;
              A(i, i) = 1;
          end
      end
      A(nmasses+1, nmasses) = -1;
      C = diag(constants);
      K = A'*C*A;
      [U, S, Vt, c, I] = svdd(K);
      fprintf("Singular Values: \n")
      disp(diag(S))
      fprintf("\nEigenvalues: \n")
      disp((diag(S)).^2)
      fprintf("\nCondition Number:\n" )
      disp(c)
      u = inv(K)*f';
      e = A*u;
      w = C*e;
  end

end

function [U, S, Vt, c, I] = svdd(A)
    AAT = A * A';
    [v1,d1] = eig(AAT);
    ATA = A' *A;
    [v2,d2] = eig(ATA);
    d1 = diag(d1);
    d2 = diag(d2);
    [d1, order1] = sort(d1, 'descend');
    [d2, order2] = sort(d2, 'descend');
    S = sqrt(diag(d2));
    svals = diag(S);
    sort(S, 'descend');
    U = v1(:, order1);
    V = v2(:, order2);
    Vt = V';
    c = max(svals) / min(svals);
    if any(svals == 0)
        error ("Matrix has no inverse.");
    else
       %% Calculating the inverse of A given the SVD
       I = V*inv(S)*U';
    end
end
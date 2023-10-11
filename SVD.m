%% SVD Decomposition

%% Invertible Matrix
A = [1 2 3;
     0 1 4;
     5 6 0];
%% Calling the SVD
%% Returns U, S, and V transpose of the SVD 
%% Returns c which is the condition number
%% Returns I which is A inverse
[U, S, Vt, c, I] = svdd(A)

disp("Built-in SVD and Inverse")
[U1, S1, V1] = svd(A);
U1
S1
VT = V1'
AI = inv(A)


function [U, S, Vt, c, I] = svdd(A)
    %% Calculating A * A transpose
    AAT = A * A';

    %% Finding v1 whose columns are the corresponding right eigenvectors
    %% Finding d1, a diagonal matrix of eigenvalues
    %% Both are of AAT
    [v1,d1] = eig(AAT);

    %% Calculating A tranpose * A
    ATA = A' *A;

    %% Same as v1 and d1 from above, except they're for ATA
    [v2,d2] = eig(ATA);

    %% Getting the eigenvalues from the diagonal eigenvalue matrices
    d1 = diag(d1);
    d2 = diag(d2);

    %% Sorting the eigenvalues in descending order and keeping track of the order
    [d1, order1] = sort(d1, 'descend');
    [d2, order2] = sort(d2, 'descend');

    %% Taking the square root of the values in d2 to get a matrix of singular values
    S = sqrt(diag(d2));

    %% Getting just the singular values
    svals = diag(S);

    %% Sorting S in descending order and keeping track of the order
    sort(S, 'descend');

    %% Sorting U and V based on the order the eigenvalues were sorted with
    U = v1(:, order1);
    V = v2(:, order2);

    %% Finding the transpose of V
    Vt = V';

    %% Calculating the condition number from the singular values
    c = max(svals) / min(svals);
    
    %% Checking if any singular values are equal to 0 because that would mean
    %% there is no inversw
    %% Returns an error if there is no inverse
    if any(svals == 0)
        error ("Matrix has no inverse.");
    else
       %% Calculating the inverse of A given the SVD
       I = V*inv(S)*U';
    end

end

        

    






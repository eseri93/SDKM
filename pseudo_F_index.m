% Pseudo-F %

function pFdk = pseudo_F_index(X, U, V, K, Q)
    % Inputs:
    % X - (n x J) data matrix
    % U - (n x K) binary membership matrix for units
    % V - (J x Q) binary membership matrix for variables
    % K - number of clusters for units
    % Q - number of clusters for variables

    % Size of the data matrix
    [n, J] = size(X);

    % Compute projection matrices H_U and H_V
    % H_U projects the data onto the subspace defined by the clustering of units
    H_U = U / (U' * U) * U';  % (n x n) matrix

    % H_V projects the data onto the subspace defined by the clustering of variables
    H_V = V / (V' * V) * V';  % (J x J) matrix

    % Create a matrix of ones to compute the grand mean properly
    one_n = ones(n, 1);   % (n x 1) vector of ones
    one_J = ones(J, 1);   % (J x 1) vector of ones

    % Compute the grand mean matrix
    grand_mean_matrix = (1 / (n * J)) * (one_n * one_n') * X * (one_J * one_J');

    % Compute between-cluster sum of squares (SS_between)
    % SS_between compares the clustered data projected by H_U and H_V to the grand mean
    X_clustered = H_U * X * H_V;
    SS_between = norm(X_clustered - grand_mean_matrix, 'fro')^2;
    
    % Compute within-cluster sum of squares (SS_within)
    % SS_within compares the actual data to the clustered data (projected data)
    SS_within = norm(X - X_clustered, 'fro')^2;
    
    % Degrees of freedom
    df_between = (K * Q) - 1;
    df_within = (n * J) - (K * Q);
    
    % Compute pseudo-F index
    pFdk = (SS_between / df_between) / (SS_within / df_within);
end


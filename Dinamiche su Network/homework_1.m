W = [0,1,1,1,1,1,0,0,0,0,0,0,0,0; ...
     1,0,1,1,1,1,0,0,0,0,0,0,0,0; ...
     1,1,0,1,1,1,0,0,0,0,0,0,0,0; ...
     1,1,1,0,1,1,0,0,0,0,0,0,0,0; ...
     1,1,1,1,0,1,0,0,0,0,0,0,0,0; ...
     1,1,1,1,1,0,1,0,0,0,0,0,0,0; ...
     0,0,0,0,0,1,0,1,0,0,0,0,0,0; ...
     0,0,0,0,0,0,1,0,1,0,0,0,0,0; ...
     0,0,0,0,0,0,0,1,0,1,1,1,1,1; ...
     0,0,0,0,0,0,0,0,1,0,0,0,0,0; ...
     0,0,0,0,0,0,0,0,1,0,0,0,0,0; ...
     0,0,0,0,0,0,0,0,1,0,0,0,0,0; ...
     0,0,0,0,0,0,0,0,1,0,0,0,0,0; ...
     0,0,0,0,0,0,0,0,1,0,0,0,0,0];
G = graph(W);
plot(G);

% DEGREE CENTRALITY
% degreeCentrality = centrality(G, 'degree', 'Importance', G.Edges.Weight);
degrees = zeros(length(W), 1);
for i = 1:length(W)
    degrees(i) = sum(W(:, i));
end
degreeCentrality = degrees;

% EIGENVECTOR CENTRALITY
% eigenvectorCentrality = centrality(G, 'eigenvector', 'Importance', G.Edges.Weight);
[~, eigenvalues] = eig(W);
lambdaW = max(eigenvalues, [], 'all');
[eigenvectors, eigenvalues] = eig((1/lambdaW) * transpose(W));
eigenvalues = diag(eigenvalues);
eigenvector = eigenvectors(:, eigenvalues == max(eigenvalues, [], 'all')); % prendo l'autovettore di autovalore 1
eigenvectorCentrality = (1/sum(eigenvector) * eigenvector); % normalizzo rispetto alla somma

% INVARIANT DISTRIBUTION CENTRALITY
D = diag(degrees);
P = D^(-1) * W;
[eigenvectorsP, eigenvaluesP] = eig(P);
eigenvaluesP = diag(eigenvaluesP);
eigenvectorP = eigenvectorsP(:, eigenvaluesP == max(eigenvaluesP)); % prendo l'autovettore di autovalore 1
invariantDistributionCentrality = (1/sum(eigenvectorP) * eigenvectorP); % normalizzo rispetto alla somma

% KATZ CENTRALITY
beta = 0.15;
mu = ones(length(W), 1);

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
% Check: centrality(G, 'degree', 'Importance', G.Edges.Weight)
degrees = zeros(length(W), 1);
for i = 1:length(W)
    degrees(i) = sum(W(:, i));
end
degreeCentrality = degrees;

% EIGENVECTOR CENTRALITY
% Check: centrality(G, 'eigenvector', 'Importance', G.Edges.Weight);
[~, eigenvalues] = eig(W);
lambdaW = max(eigenvalues, [], 'all');
[eigenvectors, eigenvalues] = eig((1 / lambdaW) * transpose(W));
eigenvalues = diag(eigenvalues);
eigenvector = eigenvectors(:, eigenvalues == max(eigenvalues, [], 'all')); % eigenvector with lambda=1
eigenvectorCentrality = ((1 / sum(eigenvector)) * eigenvector); % normalize: 1 / sum(eigenvec elements)

% INVARIANT DISTRIBUTION CENTRALITY
D = diag(degrees);
P = D^(-1) * W;
[eigenvectorsP, eigenvaluesP] = eig(transpose(P));
eigenvaluesP = diag(eigenvaluesP);
eigenvectorP = eigenvectorsP(:, eigenvaluesP == max(eigenvaluesP, [], 'all')); % eigenvector with lambda=1
invariantDistributionCentrality = ((1 / sum(eigenvectorP)) * eigenvectorP); % normalize: 1 / sum(eigenvec elements)

% KATZ CENTRALITY
% Check: (eye(14) - (((1 - beta) / lambdaW) * transpose(W)) )^(-1) * (beta * mu)
beta = 0.15;
mu = ones(length(W), 1);
katzCentrality = zeros(length(W), 1);
while true
    oldCentrality = katzCentrality;
    katzCentrality = ((( (1 - beta) / lambdaW) * transpose(W)) * katzCentrality) + (beta * mu);
    if katzCentrality == oldCentrality
        break;
    end
end

% PAGERANK


function [muCluster] = computeClusterSpectra(A, S_abstract, numCluster)

%% Eigendecomposition-based routine
[U,Sigma,V] = svd(A,'econ');
Lambda = diag(Sigma).^2;
C = U;
R = Sigma*V';

C_reduce = C(:,1:S_abstract);
R_reduce = R(1:S_abstract,:);

%% Clustering 

% Clustering
idxClusters = kmeans(R_reduce, kseeds(R_reduce, numCluster))';

% Find means
clusterCentre = zeros(size(A,1),numCluster);
for i = 1:numCluster 
    clusterCentre(:,i) = mean(A(:,(idxClusters==i)),2);
end
clusterCentre = clusterCentre;

%% Compare with ground truth

% Get rid of noise cluster
[~, noise_index] = min(mean(clusterCentre));

clusterCentreReduced = clusterCentre;
clusterCentreReduced(:, noise_index) = [];

muCluster = clusterCentreReduced;

end
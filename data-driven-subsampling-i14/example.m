%%%% Example
clear, clc

% Load a dataset as an example
load("i14-331373_xanes_autoprocessing0.mat")

% Reshape
nE = length(energy);
nY = length(y);
nX = length(x);
A = reshape(data, nE,[]); 

% Set subsampling ratio
p = 0.06;

% Load spectrum dictionary
load("spectrumDictionaryFe.mat")
spectrumDictionary = table2array(spectraFeHemMag);
spectrumDictionary = spectrumDictionary(:,2:end);

%% Let's try RISS first
[energiesToSampleRISS, spatialRowsToSampleRISS, energyIndicesRISS] = RISS(p, spectrumDictionary, energy, nX, nY);

% A little bit of work to transform this to a sparsity pattern
pixellabels = reshape(1:(nX*nY), nY, nX);
OmegaRISS = zeros(size(A)); OmegaRISS = logical(OmegaRISS);
OmegaRISS(energyIndicesRISS, :) = 1;

for ii = 1:nE
    if ~ismember(ii, energyIndicesRISS)
        pixel_indices = pixellabels(spatialRowsToSampleRISS(ii,2:end)', :);
        OmegaRISS(ii, reshape(pixel_indices, 1, [])) = 1;
    end
end
A_sparseRISS = A .* OmegaRISS;
A_sparseRISStensor = reshape(A_sparseRISS, [nE, nY, nX]);

% One can now do loopedASD on this sparse A

% We show the sparsity pattern of two consecutive frames
figure(1)
subplot(4,2,1), imshow(reshape(A_sparseRISStensor(100, :, :), nY, nX))
title("Sparsity pattern frame 100 from RISS")
subplot(4,2,2), imshow(reshape(A_sparseRISStensor(101, :, :), nY, nX))
title("Sparsity pattern frame 101 from RISS")

%% Let's try CUR now
[energiesToSampleCURISS, spatialRowsToSampleCURISS, energyIndicesCURISS] = CURISS(p, spectrumDictionary, energy, nX, nY);

% A little bit of work to transform this to a sparsity pattern
pixellabels = reshape(1:(nX*nY), nY, nX);
OmegaCURISS = zeros(size(A)); OmegaCURISS = logical(OmegaCURISS);
OmegaCURISS(energyIndicesCURISS, :) = 1;
pixel_indicesCUR = pixellabels(spatialRowsToSampleCURISS', :);
OmegaCURISS(:, reshape(pixel_indicesCUR, 1, [])) = 1;

A_sparseCURISS = A .* OmegaCURISS;
A_sparseCURISStensor = reshape(A_sparseCURISS, [nE, nY, nX]);

% We show the sparsity pattern of two consecutive frames
figure(1)
subplot(4,2,3), imshow(reshape(A_sparseCURISStensor(100, :, :), nY, nX))
title("Sparsity pattern frame 100 from CURISS")
subplot(4,2,4), imshow(reshape(A_sparseCURISStensor(101, :, :), nY, nX))
title("Sparsity pattern frame 101 from CURISS")

%% Reconstruct CURISS
% Get C, U, and R
C = A_sparseCURISS(:, pixel_indicesCUR);
R = A_sparseCURISS(energyIndicesCURISS, :);
U = A_sparseCURISS(energyIndicesCURISS, pixel_indicesCUR);

% Invert U in a stable manner
[Qu,Ru] = qr(U, 'econ');
A_completed = (C/Ru) * (Qu' * R);

figure(1)
subplot(4,2,5), imagesc(reshape(A_completed(100, :, :), nY, nX), [0 max(A_completed,[], 'all')]), axis off,
title("Frame 100 completed using CURISS")
subplot(4,2,6), imagesc(reshape(A(100, :, :), nY, nX), [0 max(A_completed,[], 'all')]), axis off
title("Frame 100 from full data")

%% Let's try ACURISS now

% Set parameters for ACURISS
p0 = 0.02;
maxIte = 30;
tolCompletionVar = 1e-3;
tolSpectralVar = 1e-3;
% Clustering for spectra
S_abstract = 3;
numCluster = 4;
% PostEdgeEnergyIndex
post_edge_energy = 7.2*1e3;
[~, postEdgeEnergyIndex] = min(abs(energy - post_edge_energy));
% You need codes for the cluster analysis
addpath("../clustering_codes/")

% Run ACURISS
[energiesToSampleACURISS,spatialRowsToSampleACURISS,energyIndicesACURISS] = ...
ACURISS(A,p0,spectrumDictionary,energy,nX,nY,maxIte, ...
            tolCompletionVar,tolSpectralVar,S_abstract,...
    numCluster,postEdgeEnergyIndex);

% A little bit of work to transform this to a sparsity pattern
pixellabels = reshape(1:(nX*nY), nY, nX);
OmegaACURISS = zeros(size(A)); OmegaACURISS = logical(OmegaACURISS);
OmegaACURISS(energyIndicesACURISS, :) = 1;
pixel_indicesCUR_A = pixellabels(spatialRowsToSampleACURISS', :);
OmegaACURISS(:, reshape(pixel_indicesCUR_A, 1, [])) = 1;

A_sparseACURISS = A .* OmegaACURISS;
A_sparseACURISStensor = reshape(A_sparseACURISS, [nE, nY, nX]);

%% Reconstruct ACURISS
% Get C, U, and R
C = A_sparseACURISS(:, pixel_indicesCUR_A);
R = A_sparseACURISS(energyIndicesACURISS, :);
U = A_sparseACURISS(energyIndicesACURISS, pixel_indicesCUR_A);

% Invert U in a stable manner
[Qu,Ru] = qr(U, 'econ');
A_completed = (C/Ru) * (Qu' * R);

figure(1)
subplot(4,2,7), imagesc(reshape(A_completed(100, :, :), nY, nX), [0 max(A_completed,[], 'all')]), axis off,
title("Frame 100 completed using ACURISS")
subplot(4,2,8), imagesc(reshape(A(100, :, :), nY, nX), [0 max(A_completed,[], 'all')]), axis off
title("Frame 100 from full data")

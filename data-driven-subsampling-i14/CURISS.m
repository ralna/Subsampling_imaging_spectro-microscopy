function [energiesToSample, spatialRowsToSample, energyIndices] = CURISS(p, spectrumDictionary, energies, nX, nY)
% Input:
% p                   = subsampling ratio upper bound, a scalar between 0 
%                       and 1
% spectrumDictionary  = a matrix representing a dictionary of possible
%                       spectra present in the sample. Should correspond to 
%                       the energies in the energies input. If no
%                       spectrumDictionary is available, input a vector of
%                       ones the size of 
% energies            = a vector of the possible photon energies
% nX, nY              = the number of pixels in the X and Y directions

% Output:
% energiesToSample    = the list of energies of which a full scan should be
%                       taken
% spatialRowsToSample = the list of row indices to sample at all other
%                       energies

%% Decide on the number of spatial rows and the number of energies for full scans
% Take approximately the same number of rows as of full scans
nE = length(energies);
totalmeasurements = nE * nX * nY;
numberOfEnergies = ceil(nY * nE * p/(nY + nE));
fullscans = numberOfEnergies * nX * nY;
numberOfRows = floor((p * totalmeasurements - fullscans) / ((nE - numberOfEnergies) * nX));

% Make sure there's at least one spatial row
if numberOfRows < 1
    numberOfEnergies = floor(nY * nE * p/(nY + nE));
    fullscans = numberOfEnergies * nX * nY;
    numberOfRows = floor((p * totalmeasurements - fullscans) / ((nE - numberOfEnergies) * nX));
end

% Throw an error if the subsampling ratio is too low
if numberOfRows < 1 || numberOfEnergies < 1
    error("Subsampling ratio too small")
end


%% Identify important energies
% Perform a leverage score based investigation into important energies -
% once an energy has been identified, the surrounding energies are taken
% out of contention

% Orthogonal decomposition of the spectrum
[Q,~] = qr(spectrumDictionary, 'econ');
% Compute leverage scores
lisq = vecnorm(Q, 2, 2).^2;
% Randomly sample energies according to the leverage scores
energyIndices = datasample(1:nE, numberOfEnergies, 'Weights', lisq, 'Replace', false);
% Return the energies where full scans should be taken
energiesToSample = energies(energyIndices);
    
%% Load data from full scans for some images
% Full scans should be taken at the sampled energies and the XRF scans
% should be placed next to eachother - we load scans from existing data
% here as an example
load("i14-331373_xanes_autoprocessing0.mat")
A = reshape(data, length(energies),[]); 
row_wise_data = [];
for indx = energyIndices
    row_wise_data = [row_wise_data reshape(A(indx, :), nY, nX)];
end

%% Low-rank decomposition using ARP
[U,~,~] = svd(row_wise_data, "econ");
Uk = U(:, 1:numberOfRows);
spatialRowsToSample = sort(ARP(Uk));     

end
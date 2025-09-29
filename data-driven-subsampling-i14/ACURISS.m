function [energiesToSample,spatialRowsToSample,energyIndices] = ...
    ACURISS(Afull,p0,spectrumDictionary,energies,nX,nY,maxIte, ...
            tolCompletionVar,tolSpectralVar,S_abstract,...
    numCluster,postEdgeEnergyIndex)
% Input:
% Afull               = the whole data-set. This is never accessed entirely,
%                       but it is needed in input in order to extract from 
%                       it energies and spatial rows, and to define 
%                       different unfolding of it.
% p0                  = intial subsampling ratio upper bound, a scalar 
%                       between 0 and 1
% spectrumDictionary  = a matrix representing a dictionary of possible
%                       spectra present in the sample. Should correspond to 
%                       the energies in the energies input. If no
%                       spectrumDictionary is available, input a vector of
%                       ones the size of 
% energies            = a vector of the possible photon energies
% nX, nY              = the number of pixels in the X and Y directions
% maxIte              = the number of maximum adaptive cycles. Used in case
%                       where the tolerances for the stopping criterias are
%                       not satisfied 
% tolCompletionVar    = the tolerance for the stopping criteria defined by
%                       the completion variation (if not insert the criteria
%                       is ignored)
% tolSpectralVar      = the tolerance for the stopping criteria defined by
%                       the spectral variation (if not insert the criteria
%                       is ignored)
% S_abstract          =
% numCluster          = number of clusters in the cluster analysis (needed
%                       for computing the Spectral Variantion)
% postEdgeEnergyIndex =
%
%
%
% Output:
% energiesToSample    = the list of energies of which a full scan should be
%                       taken
% spatialRowsToSample = the list of row indices to sample at all other
%                       energies
%
%
% Needs the functions: ARP, computeClusterSpectra, compareUpToColPerm


%% A first iteration of CURISS

% ------
% Data aware Sampling 
% (Afull, nE, nX, nY, p0, muTrue) --> [lisq, energy_indices, srow_indices]
    
    % Decide on the number of spatial rows and the number of energies for full scans
        % Take approximately the same number of rows as of full scans
        nE = length(energies);
        totalmeasurements = nE * nX * nY;
        numberOfEnergies = ceil(nY * nE * p0/(nY + nE));
        fullscans = numberOfEnergies * nX * nY;
        numberOfRows = floor((p0 * totalmeasurements - fullscans) / ((nE - numberOfEnergies) * nX));
        
        if numberOfRows < 1
            numberOfEnergies = floor(nY * nE * p0/(nY + nE));
            fullscans = numberOfEnergies * nX * nY;
            numberOfRows = floor((p0 * totalmeasurements - fullscans) / ((nE - numberOfEnergies) * nX));
        end
    
    
        % Throw an error if the subsampling ratio is too low
        if numberOfRows < 1 || numberOfEnergies < 1
            error("Subsampling ratio too small")
        end

    % Identify important energies
        % Perform a leverage score based investigation into important energies -
        % once an energy has been identified, the surrounding energies are taken
        % out of contention

        % Orthogonal decomposition of the spectrum
        [Q,~] = qr(spectrumDictionary, 'econ');
        % Compute leverage scores
        lisq = vecnorm(Q, 2, 2).^2;
        % Randomly sample energies according to the leverage scores
        energyIndices = datasample(1:nE, numberOfEnergies, 'Weights', lisq, 'Replace', false);
        % Set probabilities of measured energies to 0
        lisq(energyIndices) = 0;

    % Change unfolding

        row_wise_data = [];
        for indx = energyIndices
            row_wise_data = [row_wise_data reshape(Afull(indx, :), nY, nX)];
        end
        row_wise_data = row_wise_data / max(row_wise_data, [], 'all');
       
    % Low-rank decomposition using ARP

        [U,~,~] = svd(row_wise_data, "econ");
        Uk = U(:, 1:numberOfRows); 
    
        pixellabels = reshape(1:(nX*nY), nY, nX);
    
        spatialRowsToSample = sort(ARP(Uk));
        pixel_indices = pixellabels(spatialRowsToSample', :);

% ------
% Completion by CUR 
% (Afull, energy_indices, pixel_indices) --> [A_completed,C,U,R,resmat,
%                                                               resmatF]

    C = Afull(:, pixel_indices); 
    R = Afull(energyIndices, :);
    Up = Afull(energyIndices, pixel_indices);

    [QU,RU] = qr(Up, "econ");
    A_completed = (C/RU) * (QU' * R);
    normFirstComp = norm(A_completed,'fro');

    stock_A_completed(:,:,1) = A_completed;


% ------
% Custer Analysis 
% (A_completed, S_abstract, numCluster) --> [old_muClusterNormalised, 
%                                         resmatSpecNorm2, resmatSpecNormF] 

    muCluster = computeClusterSpectra(A_completed, S_abstract, numCluster);
    postEdgeMean = mean(muCluster(postEdgeEnergyIndex:end, :));
    muClusterNormalised = muCluster./postEdgeMean;
   
    old_muClusterNormalised = muClusterNormalised;
%% Adaptive cycles

% Check if tolerance for stopping criteria has been inputed

    if (nargin < 8 || isempty(tolCompletionVar))
    tolCompletionVar = Inf;
    end

    if (nargin < 9 || isempty(tolSpectralVar))
        tolSpectralVar = Inf;
    end

% Initialize variations and corresponding tests 
   
    spectralVar = [];
    completionVar = [];
    test_spectralVar = 1;
    test_completionVar = 1;
    numIte = 0;

% ------
% Check critiria 

    while (test_spectralVar > tolSpectralVar && ... % Check spectral variation
           test_completionVar > tolCompletionVar && ... % Check completion variation
           numIte < maxIte) % Check number of iteration already done
        
        % Count number of iterations performed
        numIte = numIte + 1;

% ------
% Find new indices    
                % Add one new energy scan
                % (nE, lisq, energy_indices, Afull) --> 
                % [energy_indices, lisq, R, Up]

                new_energy_idx = datasample(1:nE,1,'Weights', lisq); % Note: lisq is 0 for already selected energies
                energyIndices = [energyIndices new_energy_idx];
                energyIndices = sort(energyIndices);

                lisq(new_energy_idx) = 0;

                % Update C, U

                R = Afull(energyIndices, :);
                Up = Afull(energyIndices, pixel_indices);

%% Computing stopping criterias
% Completion
        
        % Form CUR
        [QU,RU] = qr(Up, "econ");
        A_completed = (C/RU) * (QU' * R);

        stock_A_completed(:,:,numIte +1) = A_completed;

        % ------- Completion variation
        completionVar_current = norm(A_completed - stock_A_completed(:,:,numIte),'fro');
        
        completionVar = [completionVar completionVar_current];

        test_completionVar = completionVar_current;


% Spectral

        % Cluster Analysis
        muCluster = computeClusterSpectra(A_completed, S_abstract, numCluster);
        postEdgeMean = mean(muCluster(postEdgeEnergyIndex:end, :));
        muClusterNormalised = muCluster./postEdgeMean;

        % ------- Spectral variation
        [~, spectralVar_current] =compareUpToColPerm(muClusterNormalised, old_muClusterNormalised,1);        
        
        spectralVar = [spectralVar spectralVar_current];
    
        test_spectralVar = spectralVar_current;

        old_muClusterNormalised = muClusterNormalised;
    
    end
    
    % Return the energies where full scans should be taken
    energiesToSample = energies(energyIndices);
end



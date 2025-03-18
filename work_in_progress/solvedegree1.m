%% SOLVEDEGREE1
% Recomputes the GRACE degree-1 surface mass density coefficients based on Sun et al. (2016) with the specified GIA model.
%
% Inputs
%   Pcenter - Data centre
%       - 'CSR': Center for Space Research
%       - 'GFZ': GeoForschungsZentrum Potsdam
%       - 'JPL': Jet Propulsion Laboratory
%       The default data center is 'CSR'.
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06'.
%       The default release level is 'RL06'.
%   L - The bandwidth of the GRACE product
%       The default L is 60.
%   Lsle - The bandwidth used to solve the sea level equation
%       The default Lsle is 96.
%   GiaModel - GIA model
%       The default GIA model is 'ice6gd'.
%   OceanDomain - Ocean domain
%       The input format should be whatever KERNELCP_NEW accepts, e.g. a GeoDomain object.
%       The default ocean domain is the global ocean (ALLOCEANS) with a 0.5° buffer.
%   ReferenceEpoch - Reference epoch
%       When specified, the coefficients will be removed by the mean of the specified time range.
%       Possible input formats include:
%       - A DATETIME object as the centre of the epoch
%       - A 2-element DATETIME array as the start and end of the epoch
%       - A 2-element cell array as the centre (DATETIME) and half-length (DURATION) of the epoch
%       The defualt reference epoch is centred at Jan 1, 2008, with a 5-year half-length (similar to TN-13).
%   IncludeC20 - Whether to also recompute C20
%       The default option is false, as C20 is usually replaced with SLR values.
%   ReplaceWithGAD - Whether to remove GAD instead of GAC
%       The default option is true (see Sun et al., 2016).
%   Method - Method to solve the sea level equation
%       - 'fingerprint': Solve the sea level equation using the fingerprint method (see SOLVESLE).
%       - 'uniform': Assume a uniform barystatic load over the ocean domain.
%       The default method is 'fingerprint'.
%   ForceNew - Whether to force new generation of a save file
%       The default option is false.
%   BeQuiet - Whether to suppress output
%       The default option is soft true.
%
% See also
%   SOLVESLE, GRACE2PLMT
%
% Authored by
%   2025/03/18, williameclee@arizona.edu (@williameclee)

function varargout = solvedegree1(varargin)
    %% Initialisation
    % Parse inputs
    [pcenter, rlevel, L, Lsle, giaModel, oceanDomain, ...
         refEpochLim, includeC20, replaceWGad, method, forceNew, beQuiet] = ...
        parseinputs(varargin);

    [dataPath, ~, dataExists] = getoutputfile(pcenter, rlevel, replaceWGad, giaModel, includeC20, L, Lsle, oceanDomain, method);

    [graceSpht, ~, dates] = grace2plmt_new(pcenter, rlevel, 60, ...
        "OutputFormat", 'traditional', "BeQuiet", beQuiet);
    graceSpht = graceSpht(:, 1:4, :);

    if isnumeric(dates)
        dates = datetime(dates, "ConvertFrom", 'datenum');
    end

    if dataExists && ~forceNew
        load(dataPath, 'coeffs')
    else
        % Add back GAC and remove GAD instead (Sun et al., 2016)
        if replaceWGad
            [~, gacSpht] = aod1b2plmt(pcenter, rlevel, 'GAC', L, ...
                "OutputFormat", 'traditional', "BeQuiet", beQuiet);
            [~, gadSpht] = aod1b2plmt(pcenter, rlevel, 'GAD', L, ...
                "OutputFormat", 'traditional', "BeQuiet", beQuiet);
            graceSpht(:, 3:4, :) = graceSpht(:, 3:4, :) + gacSpht(:, 3:4, :) - gadSpht(:, 3:4, :);
        end

        % Remove GIA signal for l >= 2 (Sun et al., 2016)
        giaSpht = gia2plmt(dates, giaModel, "L", L, ...
            "OutputFormat", 'traditional', "BeQuiet", beQuiet);
        graceSpht(3:end, 3:4, :) = graceSpht(3:end, 3:4, :) - giaSpht(3:end, 3:4, :);

        % Remove time mean, similar as in TN13
        if ~isempty(refEpochLim)
            isMeanEpoch = dates >= refEpochLim(1) & dates <= refEpochLim(2);
            graceSpht(:, 3:4, :) = graceSpht(:, 3:4, :) - mean(graceSpht(:, 3:4, isMeanEpoch), 3);
        end

        if includeC20
            coeffs = nan([length(dates), 4]);
            coeffsId = [2, 3, addmup(L) + 3, 4]';
        else
            coeffs = nan([length(dates), 3]);
            coeffsId = [2, 3, addmup(L) + 3]';
        end

        oceanKernel = kernelcp_new(L, oceanDomain, "BeQuiet", beQuiet);
        landKernel = eye(size(oceanKernel)) - oceanKernel;
        coeffsKernel = oceanKernel(2:1 + length(coeffsId), 2:1 + length(coeffsId));
        oceanKernelSle = kernelcp_new(Lsle, oceanDomain, "BeQuiet", true);
        [~, ~, ~, ~, ~, oceanFunSph] = geoboxcap(Lsle, oceanDomain, "BeQuiet", beQuiet);

        for iDate = 1:length(dates)
            fprintf('Processing date %3d of %3d\n', iDate, length(dates));
            coeffs(iDate, :) = ...
                solvedegree1_iter(graceSpht(:, 3:4, iDate), oceanDomain, L, Lsle, ...
                coeffsKernel, coeffsId, oceanKernel, landKernel, oceanKernelSle, oceanFunSph, ...
                method);
        end

        save(dataPath, 'coeffs')

        if beQuiet <= 1
            fprintf('%s saved %s\n', upper(mfilename), dataPath)
        end

    end

    if nargout >= 2
        varargout = {dates, coeffs};
        return
    end

    graceSphtAux = reshape(graceSpht, [size(graceSpht, 1) * size(graceSpht, 2), size(graceSpht, 3)]);
    graceSphtAux(addmup(L) * 2 + coeffsId, :) = coeffs';
    graceSpht = reshape(graceSphtAux, size(graceSpht));

    if ~isempty(refEpochLim)
        isMeanEpoch = dates >= refEpochLim(1) & dates <= refEpochLim(2);
        graceSpht(:, 3:4, isMeanEpoch) = ...
            graceSpht(:, 3:4, isMeanEpoch) + mean(graceSpht(:, 3:4, isMeanEpoch), 3);
    end

    varargout = {dates, coeffs, graceSpht};

end

%% Subfunctions
function [coeffs, graceSph] = ...
        solvedegree1_iter(graceSph, oceanDomain, L, Lsle, coeffsKernel, coeffsId, ...
        oceanKernel, landKernel, oceanKernelSle, oceanFunSph, method)
    maxIter = 5;

    graceNoUnknownSph = graceSph;
    graceNoUnknownSph(coeffsId) = 0;
    graceNoUnknownLclSph = localise(graceNoUnknownSph, oceanDomain, L, "K", oceanKernel);
    oceanCoeffsNoUnkown = graceNoUnknownLclSph(coeffsId);

    for iIter = 1:maxIter
        landSph = localise(graceSph, oceanDomain, L, "K", landKernel);

        switch method
            case 'fingerprint'
                [~, oceanSph] = solvesle(landSph, Lsle, "Ocean", oceanDomain, ...
                    "OceanKernel", oceanKernelSle, "OceanFUnction", oceanFunSph, ...
                    "RotationFeedback", false);
                oceanSph = oceanSph(1:addmup(L), :); % Truncate to original L
                oceanCoeffs = oceanSph(coeffsId);
            case 'uniform'
                barystaticLoad = -landSph(1, 1) / oceanDomain.SphArea;
                oceanCoeffs = barystaticLoad * oceanKernel(2:length(coeffsId) + 1, 1);
        end

        coeffs = coeffsKernel \ (oceanCoeffs - oceanCoeffsNoUnkown);
        graceSph(coeffsId) = coeffs;

    end

end

function parseinputs(inputs)
    % Set default parameters
    defaultPcenter = 'CSR';
    defaultRlevel = 'RL06';
    defaultL = 60;
    defaultLsle = 96;
    defaultGiaModel = 'ice6gd';
    defaultOceanDomain = GeoDomain('alloceans', "Buffer", 0.5);
    defaultRefEpoch = {datetime(2008, 01, 01), years(5)}; % similar to TN-13
    defaultIncludeC20 = false;
    refaultReplaceWGad = true;
    defaultMethod = 'fingerprint';
    % Construct input parser
    ip = inputParser;
    addOptional(ip, 'Pcenter', defaultPcenter, ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', defaultRlevel, ...
        @(x) ischar(validatestring(x, {'RL04', 'RL05', 'RL06'})));
    addOptional(ip, 'L', defaultL, ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addOptional(ip, 'Lsle', defaultLsle, ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addOptional(ip, 'GiaModel', defaultGiaModel, ...
        @(x) ischar(x));
    addOptional(ip, 'OceanDomain', defaultOceanDomain, ...
        @(x) isa(x, 'GeoDomain'));
    addOptional(ip, 'ReferenceEpoch', defaultRefEpoch);
    addOptional(ip, 'IncludeC20', defaultIncludeC20, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'ReplaceWithGAD', refaultReplaceWGad, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'Method', defaultMethod, ...
        @(x) ischar(validatestring(x, {'fingerprint', 'uniform'})));
    addParameter(ip, 'ForceNew', false, ...
        @(x) isnumeric(x) || islogical(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) isnumeric(x) || islogical(x));
    parse(ip, inputs{:});

    pcenter = ip.Results.Pcenter;
    rlevel = ip.Results.Rlevel;
    L = ip.Results.L;
    Lsle = ip.Results.Lsle;
    giaModel = ip.Results.GiaModel;
    oceanDomain = ip.Results.OceanDomain;
    refEpoch = ip.Results.ReferenceEpoch;
    includeC20 = logical(ip.Results.IncludeC20);
    replaceWGad = logical(ip.Results.ReplaceWithGAD);
    method = ip.Results.Method;
    forceNew = ip.Results.ForceNew;
    beQuiet = ip.Results.BeQuiet * 2;

    % Obtain reference epoch
    refEpochHLength = years(5);

    if isempty(refEpoch) || (isscalar(refEpoch) && (isnan(refEpoch) || isnat(refEpoch)))
        refEpochLim = [];
    elseif isscalar(refEpoch)

        if isnumeric(refEpoch)
            refEpoch = datetime(refEpoch, "ConvertFrom", 'datenum');
        end

        refEpochCentre = refEpoch;

        refEpochLim = refEpochCentre + refEpochHLength * [-1, 1];
    elseif isdatetime(refEpoch) && length(refEpoch) == 2
        refEpochLim = refEpoch;
    elseif iscell(refEpoch) && length(refEpoch) == 2
        refEpochCentre = refEpoch{1};
        refEpochHLength = refEpoch{2};

        if isnumeric(refEpochCentre)
            refEpochCentre = datetime(refEpochCentre, "ConvertFrom", 'datenum');
        end

        if isnumeric(refEpochHLength)
            refEpochHLength = years(refEpochHLength);
        end

        refEpochLim = refEpochCentre + refEpochHLength * [-1, 1];
    end

    varargout = {pcenter, rlevel, L, Lsle, giaModel, oceanDomain, refEpochLim, includeC20, replaceWGad, method, forceNew, beQuiet};
end

function [outputPath, outputFile, outputExists] = ...
        getoutputfile(pcenter, rlevel, replaceWGad, giaModel, ...
        includeC20, L, Lsle, oceanDomain, method)
    outputFolder = fullfile(getenv('GRACEDATA'), 'Degree1');

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
        fprintf('%s created folder %s\n', upper(mfilename), outputFolder);
    end

    if replaceWGad
        replaceWGadFlag = '_RGAD';
    else
        replaceWGadFlag = '';
    end

    if includeC20
        includeC20Flag = '-WC20';
    else
        includeC20Flag = '';
    end

    switch method
        case 'fingerprint'
            methodFlag = '';
        case 'uniform'
            methodFlag = '-uniform';
        otherwise
            methodFlag = sprintf('-%s', method);
    end

    outputFile = sprintf('%s_%s%s-%s%s-L%d_Lsle%d-%s%s.mat', ...
        pcenter, rlevel, replaceWGadFlag, giaModel, includeC20Flag, L, Lsle, oceanDomain.Id, methodFlag);

    outputPath = fullfile(outputFolder, outputFile);

    outputExists = exist(outputPath, 'file');
end

%% SOLVEDEGREE1
% Recomputes the GRACE degree-1 surface mass density coefficients based on
% Sun et al. (2016) with the specified GIA model.
%
% Syntax
%   coeffs = SOLVEDEGREE1(pcenter, rlevel, Ldata, Lsle, GIAModel, oceanDomain)
%   coeffs = SOLVEDEGREE1(__, ...
%       "ReplaceWithGAD", rwGad, "IncludeC20", includeC20, "TimeRange", timelim)
%   coeffs = SOLVEDEGREE1(__, ...
%       "ForceNew", forceNew, "SaveData", saveData, "BeQuiet", beQuiet)
%   [dates, coeffs, graceSpht] = SOLVEDEGREE1(__)
%
% Input arguments
%   Pcenter - Data centre
%       - 'CSR': Center for Space Research
%       - 'GFZ': GeoForschungsZentrum Potsdam
%       - 'JPL': Jet Propulsion Laboratory
%       The default data center is 'CSR'.
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06'.
%       The default release level is 'RL06'.
%   Ldata - The bandwidth of the GRACE product
%       The default Ldata is 60.
%   Lsle - The bandwidth used to solve the sea level equation
%       The default Lsle is 96.
%   GiaModel - GIA model
%       The default GIA model is 'ice6gd'.
%   OceanDomain - Ocean domain
%       The input format should be whatever KERNELCP_NEW accepts, e.g. a
%       GeoDomain object.
%       The default ocean domain is the global ocean (ALLOCEANS) with a
%       0.5° buffer.
%   ReferenceEpoch - Reference epoch
%       When specified, the coefficients will be removed by the mean of the
%       specified time range.
%       Possible input formats include (but not necessarily limited to):
%       - A DATETIME/DATENUM object as the centre of the epoch
%       - A 2-element DATETIME/DATENUM array as the start and end of the
%           epoch
%       - A 2-element cell array as the centre (DATETIME/DATENUM) and
%           half-length (DURATION/DATENUM in years) of the epoch
%       The defualt reference epoch is centred at Jan 1, 2008, with a
%           5-year half-length (similar to TN-13).
%   IncludeC20 - Whether to also recompute C20
%       The default option is false, as C20 is usually replaced with SLR
%       values.
%   ReplaceWithGAD - Whether to remove GAD instead of GAC
%       The default option is true (see Sun et al., 2016).
%   Method - Method to solve the sea level equation
%       - 'fingerprint': Solve the sea level equation using the fingerprint
%           method (see SOLVESLE).
%       - 'uniform': Assume a uniform barystatic load over the ocean
%           domain.
%       The default method is 'fingerprint'.
%   ForceNew - Whether to force new generation of a save file
%       The default option is false.
%   BeQuiet - Whether to suppress output
%       The default option is soft true.
%
% Output arguments
%   dates - dates of the coefficients
%       In DATETIME format.
%   coeffs - C10, C11, S11 (and maybe C20) surface mass density
%       coefficients
%       Same unit as GRACE2PLMT with SD.
%   graceSpht - GRACE spherical harmonics coefficients with the recomputed
%       coefficients
%
% Example
%   coeffs = SOLVEDEGREE1('CSR', 'RL06', 60, 96, ...
%       'ice6gd', GeoDomain('alloceans', "Buffer", 2))
%   coeffs = SOLVEDEGREE1('CSR', 'RL06', 45, 180, ...
%       'ice6gd', GeoDomain('alloceans', "Buffer", 2), ...
%       "IncludeC20", true, "ReplaceWithGAD", true)
%
% See also
%   SOLVESLE, GRACE2PLMT
%
% Authored by
%   2025/03/18, williameclee@arizona.edu (@williameclee)
%
% Last modified by
%   2025/05/13, williameclee@arizona.edu (@williameclee)

function varargout = solvedegree1(varargin)
    %% Initialisation
    % Parse inputs
    [pcenter, rlevel, Ldata, Lsle, giaModel, oceanDomain, ...
         refEpochLim, includeC20, replaceWGad, method, timelim, ...
         forceNew, saveData, beQuiet] = ...
        parseinputs(varargin);

    [dataPath, ~, dataExists] = ...
        getoutputfile(pcenter, rlevel, replaceWGad, giaModel, ...
        includeC20, Ldata, Lsle, oceanDomain, method);

    [graceSpht, ~, dates] = grace2plmt_new(pcenter, rlevel, 60, ...
        "OutputFormat", 'traditional', "BeQuiet", beQuiet);
    graceSpht = graceSpht(1:addmup(Ldata), 1:4, :);

    if Ldata < Lsle
        [order, degree] = addmon(Lsle);
        graceSpht(1:addmup(Lsle), 1, :) = ...
            repmat(degree, [1, 1, size(graceSpht, 3)]);
        graceSpht(1:addmup(Lsle), 2, :) = ...
            repmat(order, [1, 1, size(graceSpht, 3)]);
    end

    if isnumeric(dates)
        dates = datetime(dates, "ConvertFrom", 'datenum');
    end

    if includeC20
        coeffsId = [2, 3, addmup(Lsle) + 3, 4]';
    else
        coeffsId = [2, 3, addmup(Lsle) + 3]';
    end

    if dataExists && ~forceNew
        load(dataPath, 'coeffs')

        if beQuiet <= 1
            fprintf('%s loaded %s\n', upper(mfilename), dataPath)
        end

    else
        coeffs = computeDeg1s(pcenter, rlevel, graceSpht, replaceWGad, Lsle, dates, giaModel, refEpochLim, includeC20, oceanDomain, coeffsId, method, beQuiet);

        if saveData
            save(dataPath, 'coeffs')

            if beQuiet <= 1
                fprintf('%s saved %s\n', upper(mfilename), dataPath)
            end

        end

    end

    if ~isempty(timelim)
        isTimeRange = dates >= timelim(1) & dates <= timelim(2);
        dates = dates(isTimeRange);
        coeffs = coeffs(isTimeRange, :);
        graceSpht = graceSpht(:, :, isTimeRange);
    end

    if nargout == 1
        varargout = {coeffs};
        return
    elseif nargout == 2
        varargout = {dates, coeffs};
        return
    end

    graceSpht(2, 3, :) = coeffs(:, 1);
    graceSpht(3, 3, :) = coeffs(:, 2);
    graceSpht(3, 4, :) = coeffs(:, 3);

    if includeC20
        graceSpht(4, 3, :) = coeffs(:, 4);
    end

    if ~isempty(refEpochLim)
        isMeanEpoch = dates >= refEpochLim(1) & dates <= refEpochLim(2);
        graceSpht(:, 3:4, isMeanEpoch) = ...
            graceSpht(:, 3:4, isMeanEpoch) + mean(graceSpht(:, 3:4, isMeanEpoch), 3);
    end

    varargout = {dates, coeffs, graceSpht};

end

%% Subfunctions
function varargout = parseinputs(inputs)
    % Set default parameters
    defaultPcenter = 'CSR';
    defaultRlevel = 'RL06';
    defaultL = 45;
    defaultLsle = 96;
    defaultGiaModel = 'ice6gd';
    defaultOceanDomain = GeoDomain('alloceans', "Buffer", 1.5); % Sun et al. (2016)
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
    addOptional(ip, 'Ldata', defaultL, ...
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
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || isdatetime(x) || isnumeric(x));
    addParameter(ip, 'ForceNew', false, ...
        @(x) isnumeric(x) || islogical(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) isnumeric(x) || islogical(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) isnumeric(x) || islogical(x));
    parse(ip, inputs{:});

    pcenter = ip.Results.Pcenter;
    rlevel = ip.Results.Rlevel;
    Ldata = ip.Results.Ldata;
    Lsle = ip.Results.Lsle;
    giaModel = ip.Results.GiaModel;
    oceanDomain = ip.Results.OceanDomain;
    refEpoch = ip.Results.ReferenceEpoch;
    includeC20 = logical(ip.Results.IncludeC20);
    replaceWGad = logical(ip.Results.ReplaceWithGAD);
    method = ip.Results.Method;
    timelim = ip.Results.TimeRange;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = ip.Results.BeQuiet * 2;

    % Obtain reference epoch
    refEpochHLength = years(5);

    if isempty(refEpoch) || ...
            (isscalar(refEpoch) && (isnan(refEpoch) || isnat(refEpoch)))
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

    varargout = ...
        {pcenter, rlevel, Ldata, Lsle, giaModel, oceanDomain, refEpochLim, ...
         includeC20, replaceWGad, method, timelim, ...
         forceNew, saveData, beQuiet};
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
        pcenter, rlevel, replaceWGadFlag, giaModel, ...
        includeC20Flag, L, Lsle, oceanDomain.Id, methodFlag);

    outputPath = fullfile(outputFolder, outputFile);

    outputExists = exist(outputPath, 'file');
end

function coeffs = computeDeg1s(pcenter, rlevel, graceSpht, replaceWGad, Lsle, dates, giaModel, refEpochLim, includeC20, oceanDomain, coeffsId, method, beQuiet)
    % Add back GAC and remove GAD instead (Sun et al., 2016)
    if replaceWGad
        [~, gacSpht] = aod1b2plmt(pcenter, rlevel, 'GAC', Lsle, ...
            "OutputFormat", 'traditional', "BeQuiet", beQuiet);
        [~, gadSpht] = aod1b2plmt(pcenter, rlevel, 'GAD', Lsle, ...
            "OutputFormat", 'traditional', "BeQuiet", beQuiet);
        graceSpht(:, 3:4, :) = graceSpht(:, 3:4, :) ...
            + gacSpht(:, 3:4, :) - gadSpht(:, 3:4, :);
    end

    % Remove GIA signal for l >= 2 (Sun et al., 2016)
    giaSpht = gia2plmt(dates, giaModel, "L", Lsle, ...
        "OutputFormat", 'traditional', "BeQuiet", beQuiet);
    graceSpht(3:end, 3:4, :) = graceSpht(3:end, 3:4, :) - giaSpht(3:end, 3:4, :);

    % Remove time mean, similar as in TN13
    if ~isempty(refEpochLim)
        isMeanEpoch = dates >= refEpochLim(1) & dates <= refEpochLim(2);
        graceSpht(:, 3:4, :) = ...
            graceSpht(:, 3:4, :) - mean(graceSpht(:, 3:4, isMeanEpoch), 3);
    end

    % Whether to also reestimate C20
    if includeC20
        coeffs = nan([length(dates), 4]);
    else
        coeffs = nan([length(dates), 3]);
    end

    % Precompute kernels
    oceanKernelSle = kernelcp_new(Lsle, oceanDomain, "BeQuiet", true);
    landKernelSle = eye(size(oceanKernelSle)) - oceanKernelSle;
    coeffsKernel = ...
        oceanKernelSle(2:1 + length(coeffsId), 2:1 + length(coeffsId));
    [~, ~, ~, ~, ~, oceanFunSph] = ...
        geoboxcap(Lsle, oceanDomain, "BeQuiet", beQuiet);
    kernelOrder = kernelorder(Lsle);

    % Decently fast, no need to use parfor
    for iDate = 1:length(dates)
        coeffs(iDate, :) = ...
            solvedegree1_iter(graceSpht(:, 3:4, iDate), oceanDomain, ...
            Lsle, coeffsKernel, coeffsId, ...
            oceanKernelSle, landKernelSle, oceanFunSph, kernelOrder, ...
            method, beQuiet);

        if beQuiet <= 1 && mod(iDate, 10) == 0
            fprintf('%s processed %s (%3d/%3d)\n', ...
                upper(mfilename), ...
                datetime(dates(iDate), "Format", 'uuuu/MM/dd'), ...
                iDate, length(dates));
        end

    end

    % Restore GAC/GAD
    if replaceWGad
        coeffs(:, 1) = ...
            coeffs(:, 1) - squeeze(gacSpht(2, 3, :) - gadSpht(2, 3, :));
        coeffs(:, 2) = ...
            coeffs(:, 2) - squeeze(gacSpht(3, 3, :) - gadSpht(3, 3, :));
        coeffs(:, 3) = ...
            coeffs(:, 3) - squeeze(gacSpht(3, 4, :) - gadSpht(3, 4, :));

        if includeC20
            coeffs(:, 4) = ...
                coeffs(:, 4) - squeeze(gacSpht(4, 3, :) - gadSpht(4, 3, :));
        end

    end

    % Restore GIA signal for l >= 2
    if includeC20
        coeffs(:, 4) = coeffs(:, 4) + ...
            squeeze(giaSpht(4, 3, :));
    end

end

function [coeffs, graceSph] = ...
        solvedegree1_iter(graceSph, oceanDomain, Lsle, coeffsKernel, coeffsId, ...
        oceanKernelSle, landKernelSle, oceanFunSph, kernelOrder, method, beQuiet)
    maxIter = 5;

    graceNoUnknownSph = graceSph;
    graceNoUnknownSph(coeffsId) = 0;
    graceNoUnknownLclSph = ...
        localise(graceNoUnknownSph, oceanDomain, Lsle, "K", oceanKernelSle);
    oceanCoeffsNoUnkown = graceNoUnknownLclSph(coeffsId);

    for iIter = 1:maxIter
        landSph = localise(graceSph, oceanDomain, Lsle, "K", landKernelSle);

        switch method
            case 'fingerprint'
                [~, oceanSph] = ...
                    solvesle(landSph, Lsle, "Ocean", oceanDomain, ...
                    "OceanKernel", oceanKernelSle, "OceanFunction", oceanFunSph, ...
                    "KernelOrder", kernelOrder, ...
                    "RotationFeedback", true, "BeQuiet", beQuiet);
                oceanSph = oceanSph(1:addmup(Lsle), :); % Truncate to original L
                oceanCoeffs = oceanSph(coeffsId);
            case 'uniform'
                barystaticLoad = -landSph(1, 1) / oceanFunSph(1, 1);
                oceanCoeffs = barystaticLoad * oceanKernelSle(2:length(coeffsId) + 1, 1);
        end

        coeffs = coeffsKernel \ (oceanCoeffs - oceanCoeffsNoUnkown);
        graceSph(coeffsId) = coeffs;
    end

end

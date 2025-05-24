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
%       The default data centre is 'CSR'.
%       When the first argument is a cell array, it is interpreted as
%       {Pcenter, Rlevel, Ldata}.
%       Data type: char
%   Rlevel - Release level of the solution
%       Either 'RL04','RL05', or 'RL06' (or numbers).
%       The default release level is 'RL06'.
%       Currently, only RL06 is guaranteed to work.
%       Data type: char | [numeric]
%   Ldata - Bandwidth of the date product
%       In the case where there are more than one product from a data
%       centre (such as BA 60 or BB 96 standard L2 products) this allows
%       you to choose between them.
%       The default L is 60.
%       Data type: [numeric]
%   Lsle - Bandwidth used to solve the sea level equation
%       The default Lsle is 96.
%       Data type: [numeric]
%   GiaModel - Name of GIA model
%       The default GIA model is 'ice6gd'.
%       Data type: char
%   OceanDomain - Ocean domain
%       The input format should be whatever KERNELCP_NEW accepts, e.g. a
%       GeoDomain object.
%       The default ocean domain is the global ocean (ALLOCEANS) with a
%       0.5° buffer.
%   IncludeC20 - Whether to also recompute C20
%       The default option is false, as C20 is usually replaced with SLR
%       values.
%       Data type: logical | ([numeric])
%   ReplaceWithGAD - Whether to remove GAD instead of GAC
%       The default option is true (see Sun et al., 2016).
%       Data type: logical | ([numeric])
%   Method - Method to solve the sea level equation
%       - 'fingerprint': Solve the sea level equation using the fingerprint
%           method (see SOLVESLE).
%       - 'uniform': Assume a uniform barystatic load over the ocean
%           domain.
%       The default method is 'fingerprint'.
%       Data type: char
%   ForceNew - Logical flag to force reprocess of the data
%       The default option is false.
%       Data type: logical | ([numeric])
%	SaveData - Logical flag to save the data
%		- true: Save the data to disk.
%		- false: Do not save the data to disk.
%		The default option is true.
%		Data types: logical | ([numeric])
%	BeQuiet - Logical flag to print messages
%		- true: Suppress all messages.
%		- false: Print all messages (usually for debugging).
%		The default option is false.
%		Data types: logical | ([numeric])
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
%   2025/05/24, williameclee@arizona.edu (@williameclee)

function varargout = solvedegree1(varargin)
    %% Initialisation
    % Parse inputs
    [pcenter, rlevel, Ldata, Ltruncation, Lsle, giaModel, oceanDomain, ...
         includeC20, replaceWGad, method, timelim, ...
         forceNew, saveData, beQuiet, onlyId] = ...
        parseinputs(varargin);

    [dataPath, ~, outputId] = ...
        getoutputfile(pcenter, rlevel, Ldata, replaceWGad, giaModel, ...
        includeC20, Ltruncation, Lsle, oceanDomain, method);

    if onlyId
        varargout = {outputId};
        return
    end

    if exist(dataPath, 'file') && ~forceNew
        load(dataPath, 'coeffs', 'dates')

        if exist('coeffs', 'var') && exist('dates', 'var')

            if beQuiet <= 1
                fprintf('%s loaded %s\n', upper(mfilename), dataPath)
            end

            varargout = formatoutput(coeffs, dates, timelim, nargout, ...
                {pcenter, rlevel, Ldata}, beQuiet);

            return
        end

    end

    [coeffs, dates] = solvedegree1Core(pcenter, rlevel, Ldata, Ltruncation, Lsle, replaceWGad, ...
        giaModel, includeC20, oceanDomain, method, beQuiet);

    if saveData
        save(dataPath, 'coeffs', 'dates')

        if beQuiet <= 1
            fprintf('%s saved %s\n', upper(mfilename), dataPath)
        end

    end

    varargout = formatoutput(coeffs, dates, timelim, nargout, ...
        {pcenter, rlevel, Ldata}, beQuiet);

end

%% Subfunctions
% Heart of the programme
function [coeffs, dates] = ...
        solvedegree1Core(pcenter, rlevel, Ldata, Ltruncation, Lsle, rwGad, giaModel, includeC20, oceanDomain, method, beQuiet)
    %% Loading data
    wbar = waitbar(0, 'Loading GRACE data', ...
        "Name", upper(mfilename), "CreateCancelBtn", 'setappdata(gcbf,''canceling'',1)');

    [gracePlmt, ~, dates] = grace2plmt_new(pcenter, rlevel, Ldata, ...
        "Unit", 'SD', "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);
    gracePlmt = ensureplmdegree(gracePlmt, Ltruncation);
    gracePlmt = ensureplmdegree(gracePlmt, Lsle);

    % Add back GAC and remove GAD instead (Sun et al., 2016)
    if rwGad
        gacPlmt = aod1b2plmt(pcenter, rlevel, 'GAC', Lsle, ...
            "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
        gacPlmt = ensureplmdegree(gacPlmt, Lsle);
        gadPlmt = aod1b2plmt(pcenter, rlevel, 'GAD', Lsle, ...
            "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
        gadPlmt = ensureplmdegree(gadPlmt, Lsle);
        gracePlmt(:, :, 3:4) = gracePlmt(:, :, 3:4) ...
            + gacPlmt(:, :, 3:4) - gadPlmt(:, :, 3:4);
    end

    % Remove GIA signal for l >= 2 (Sun et al., 2016)
    giaPlmt = gia2plmt(dates, giaModel, "L", Lsle, ...
        "OutputFormat", 'timefirst', "BeQuiet", beQuiet);
    gracePlmt(:, 3:end, 3:4) = gracePlmt(:, 3:end, 3:4) - giaPlmt(:, 3:end, 3:4);

    %% Preparing/preallocating variables
    waitbar(0, wbar, 'Preparing variables');

    if getappdata(wbar, 'canceling')
        delete(wbar);
        warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
        'Processing cancelled');
        fclose(logFid);
        return
    end

    % Whether to also reestimate C20
    coeffs = nan([length(dates), 3 + includeC20]);

    % Precompute kernels
    coeffsId = [2, 3, addmup(Lsle) + 3]';

    if includeC20
        coeffsId = [coeffsId; 4];
    end

    oceanKernelSle = kernelcp_new(Lsle, oceanDomain, ...
        "BeQuiet", beQuiet);
    landKernelSle = eye(size(oceanKernelSle)) - oceanKernelSle;
    coeffsKernel = ...
        oceanKernelSle(2:1 + length(coeffsId), 2:1 + length(coeffsId));
    [~, ~, ~, ~, ~, oceanFunPlm] = ...
        geoboxcap(Lsle, oceanDomain, "BeQuiet", beQuiet);
    kernelOrder = kernelorder(Lsle);

    %% Solving the degree-1 coefficients
    % Decently fast, no need to use parfor
    % Using partfor may run into issues with RAM
    for iDate = 1:length(dates)
        waitbar(iDate / length(dates), wbar, ...
            sprintf('Solving degree-1 coefficients (%d/%d)', iDate, length(dates)));

        if getappdata(wbar, 'canceling')
            delete(wbar);
            warning(sprintf('%s:ProcessCancelledByUser', upper(mfilename)), ...
            'Processing cancelled');
            fclose(logFid);
            return
        end

        coeffs(iDate, :) = ...
            solvedegree1Iter(squeeze(gracePlmt(iDate, :, 3:4)), ...
            oceanDomain, Lsle, coeffsKernel, coeffsId, ...
            oceanKernelSle, landKernelSle, oceanFunPlm, kernelOrder, ...
            method, beQuiet);
    end

    %% Postprocessing and output
    waitbar(1, wbar, 'Postprocessing results');

    % Restore GAC/GAD
    if rwGad
        coeffs(:, 1) = ...
            coeffs(:, 1) - squeeze(gacPlmt(:, 2, 3) - gadPlmt(:, 2, 3));
        coeffs(:, 2) = ...
            coeffs(:, 2) - squeeze(gacPlmt(:, 3, 3) - gadPlmt(:, 3, 3));
        coeffs(:, 3) = ...
            coeffs(:, 3) - squeeze(gacPlmt(:, 3, 4) - gadPlmt(:, 3, 4));

        if includeC20
            coeffs(:, 4) = ...
                coeffs(:, 4) - squeeze(gacPlmt(:, 4, 3) - gadPlmt(:, 4, 3));
        end

    end

    % Restore GIA signal for C20
    if includeC20
        coeffs(:, 4) = coeffs(:, 4) + ...
            squeeze(giaPlmt(:, 4, 3));
    end

    delete(wbar);

end

function [coeffs, graceSph] = ...
        solvedegree1Iter(graceSph, oceanDomain, Lsle, coeffsKernel, coeffsId, ...
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
                    "RotationFeedback", true, "BeQuiet", beQuiet + (beQuiet == 1));
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

% Parse input arguments
function varargout = parseinputs(inputs)
    % Set default parameters
    deOpt.Pcenter = 'CSR';
    deOpt.Rlevel = 'RL06';
    dfOpt.Ldata = 60;
    dfOpt.Ltruncation = 45;
    dfOpt.Lsle = 96;
    dfOpt.GiaModel = 'ice6gd';
    dfOpt.OceanDomain = GeoDomain('alloceans', "Buffer", 0.5); % Sun et al. (2016)
    dfOpt.IncludeC20 = false;
    dfOpt.RwGAD = true;
    dfOpt.Method = 'fingerprint';
    % Construct input parser
    ip = inputParser;
    addOptional(ip, 'Pcenter', deOpt.Pcenter, ...
        @(x) ischar(validatestring(x, {'CSR', 'GFZ', 'JPL'})));
    addOptional(ip, 'Rlevel', deOpt.Rlevel, ...
        @(x) ischar(validatestring(x, {'RL04', 'RL05', 'RL06'})));
    addOptional(ip, 'Ldata', dfOpt.Ldata, ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addOptional(ip, 'Ltruncation', dfOpt.Ltruncation, ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addOptional(ip, 'Lsle', dfOpt.Lsle, ...
        @(x) isscalar(x) && isnumeric(x) && x > 0);
    addOptional(ip, 'GiaModel', dfOpt.GiaModel, ...
        @(x) ischar(x));
    addOptional(ip, 'OceanDomain', dfOpt.OceanDomain, ...
        @(x) isa(x, 'GeoDomain'));
    addOptional(ip, 'IncludeC20', dfOpt.IncludeC20, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'ReplaceWithGAD', dfOpt.RwGAD, ...
        @(x) islogical(x) || isnumeric(x));
    addOptional(ip, 'Method', dfOpt.Method, ...
        @(x) ischar(validatestring(x, {'fingerprint', 'uniform'})));
    addOptional(ip, 'TimeRange', [], ...
        @(x) isempty(x) || isdatetime(x) || isnumeric(x));
    addParameter(ip, 'ForceNew', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'SaveData', true, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'BeQuiet', 0.5, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    addParameter(ip, 'OnlyId', false, ...
        @(x) (isnumeric(x) || islogical(x)) && isscalar(x));
    parse(ip, inputs{:});

    pcenter = ip.Results.Pcenter;
    rlevel = ip.Results.Rlevel;
    Ldata = ip.Results.Ldata;
    Ltruncation = ip.Results.Ltruncation;
    Lsle = ip.Results.Lsle;
    giaModel = ip.Results.GiaModel;
    oceanDomain = ip.Results.OceanDomain;
    includeC20 = logical(ip.Results.IncludeC20);
    replaceWGad = logical(ip.Results.ReplaceWithGAD);
    method = ip.Results.Method;
    timelim = ip.Results.TimeRange;
    forceNew = logical(ip.Results.ForceNew);
    saveData = logical(ip.Results.SaveData);
    beQuiet = double(ip.Results.BeQuiet) * 2;
    onlyId = logical(ip.Results.OnlyId);

    varargout = ...
        {pcenter, rlevel, Ldata, Ltruncation, Lsle, giaModel, oceanDomain, ...
         includeC20, replaceWGad, method, timelim, ...
         forceNew, saveData, beQuiet, onlyId};
end

% Get the input and output file names
function [outputPath, outputFile, deg1Id] = ...
        getoutputfile(pcenter, rlevel, Ldata, replaceWGad, giaModel, ...
        includeC20, Ltruncation, Lsle, oceanDomain, method)
    outputFolder = fullfile(getenv('GRACEDATA'), 'Degree1', 'new');

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
        fprintf('%s created folder %s\n', upper(mfilename), outputFolder);
    end

    gadFlag = '';

    if replaceWGad
        gadFlag = 'GAD';
    end

    includeC20Flag = '';

    if includeC20
        includeC20Flag = '-WC20';
    end

    switch method
        case 'fingerprint'
            methodFlag = '';
        case 'uniform'
            methodFlag = '-uniform';
    end

    productId = sprintf('%s%s%d', pcenter, rlevel, Ldata);
    deg1Id = sprintf('%s-%s%s-Ld%d_Ls%d-%s%s', ...
        gadFlag, giaModel, includeC20Flag, ...
        Ltruncation, Lsle, oceanDomain.Id, methodFlag);

    outputFile = sprintf('%s_%s.mat', productId, deg1Id);

    outputPath = fullfile(outputFolder, outputFile);
end

% Format the output
function output = formatoutput(coeffs, dates, timelim, nOut, product, beQuiet)

    if ~isempty(timelim)
        isTimeRange = dates >= timelim(1) & dates <= timelim(2);
        dates = dates(isTimeRange);
        coeffs = coeffs(isTimeRange, :);
    end

    if nOut <= 2
        output = {coeffs, dates};
        return
    end

    [gracePlmt, ~, dates] = grace2plmt_new(product, ...
        "Unit", 'SD', "TimeRange", timelim, "OutputFormat", 'timefirst', "TimeFormat", 'datetime', ...
        "BeQuiet", beQuiet);

    gracePlmt(:, 2, 3) = coeffs(:, 1);
    gracePlmt(:, 3, 3) = coeffs(:, 2);
    gracePlmt(:, 3, 4) = coeffs(:, 3);

    if size(coeffs, 2) == 4
        gracePlmt(:, 4, 3) = coeffs(:, 4);
    end

    output = {coeffs, dates, gracePlmt};

end

% Make sure the plm has the right degree
function plm = ensureplmdegree(plm, L)

    if size(plm, 2) == addmup(L)
        return
    end

    if size(plm, 2) > addmup(L)
        plm = plm(:, 1:addmup(L), :);
        return
    end

    [order, degree] = addmon(L);
    plm(:, 1:addmup(L), 1) = ...
        repmat(reshape(degree, [1, length(degree)]), [size(plm, 1), 1]);
    plm(:, 1:addmup(L), 2) = ...
        repmat(reshape(order, [1, length(degree)]), [size(plm, 1), 1]);

end

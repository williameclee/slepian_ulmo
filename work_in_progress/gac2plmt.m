%% GAC2PLMT
% [potcoffs, dates] = gac2plmt(Pcenter,Rlevel,units,forcenew)
%
% This program reads in the Level-2 GRACE geoid products from either the CSR or
% GFZ data centers, does some processing, and saves them as a plmt matrix
% in a .mat file.  In particular, the coefficients are reordered to our
% prefered lmcosi format, they are referenced to the WGS84 ellipsoid,
% the C2,0 coefficients are replaced with more accurate measurements from
% satellite laser ranging, and the degree one coefficients are
% substituted with those from Swenson et al. (2008).  You have the option
% of leaving them as geopotential
% or converting them to surface mass density using the method of
% Wahr et al. 1998, based on Love numbers (see PLM2POT).
%
% INPUT:
%
% Pcenter     'CSR' data center at the Center for Space Research
%             'GFZ' data center at the GeoForschungsZentrum Potsdam
% Rlevel      The release level of the solution you want.
%              Either 'RL04' or 'RL05'
% units       'POT' or 'SD' for whether you want geopotential or surface
%               mass density
% forcenew    Wether or not you want to force new generation of a save file
%              (1) or just use the one we already have (0) [default].
%
% OUTPUT:
%
% Returns these variables and saves them in a .mat file:
%    potcoffs       potential coefficients [nmonths x addmup(Ldata) x 6]
%                    these could also be in surface mass density
%    cal_errors     calibrated errors [nmonths x addmup(Ldata) x 4]
%    dates       time stamps in Matlab time
%
% NOTE:
%   The RL05 solutions from CSR do not have any standard deviations
%    given, formal or calibrated.  These values will be reported as 0.
%
%   SLR data available from the GRACE Tellus website:
%   http://grace.jpl.nasa.gov/data/J2/ notably
%   ftp://ftp.csr.utexas.edu/pub/slr/degree_2/C20_RL04_2010_12.txt
%  The header was removed and the file renamed for easy use.
%  Updated files keep getting posted in the same location.
%
%  SLR data available from the GRACE Tellus website:
%   http://grace.jpl.nasa.gov/data/degree1/ notably
%   ftp://podaac.jpl.nasa.gov/allData/tellus/L2/degree_1/
%  The header was removed and the file renamed for easy use.
%  Updated files keep getting posted in the same location.
%
% EXAMPLE: to make a new save file when you have added more months
% [potcoffs,cal_errors,dates]=grace2plmt('CSR','RL05','SD',1);
%
%
%
% Last modified by charig-at-princeton.edu, 02/27/2014
% Last modified by fjsimons-at-alum.mit.edu, 05/17/2011

function varargout = gac2plmt(Pcenter, Rlevel, units, forcenew)

    % Determine parameters and set defaults
    defval('Pcenter', 'CSR')
    defval('Rlevel', 'RL06')
    defval('units', 'SD')
    defval('forcenew', true)

    if ~isempty(getenv('ORIGINALGRACEDATA'))
        inputFolder = fullfile(getenv('ORIGINALGRACEDATA'), ...
            Rlevel, Pcenter);
    elseif ~isempty(getenv('GRACEDATA'))
        inputFolder = fullfile(getenv('GRACEDATA'), ...
            'raw', Rlevel, Pcenter);
    else
        inputFolder = fullfile(getenv('IFILES'), ...
            'GRACE', 'raw', Rlevel, Pcenter);
    end

    if ~isempty(getenv('GRACEDATA'))
        outputFolder = fullfile(getenv('GRACEDATA'));
    else
        outputFolder = fullfile(getenv('IFILES'), 'GRACE');
    end

    % And the name of that save file
    if strcmp(units, 'SD')
        outputPath = fullfile(outputFolder, ...
            sprintf('%s_%s_GAC_SD.mat', Pcenter, Rlevel));
    else
        outputPath = fullfile(outputFolder, ...
            sprintf('%s_%s_GAC.mat', Pcenter, Rlevel));
    end

    % If this file already exists, load it.  Otherwise, or if we force it, make
    % a new one (e.g. you added extra months to the database).
    if exist(outputPath, 'file') == 2 && ~forcenew
        load(outputPath, 'potcoffs', 'dates');
        fprintf('%s loading %s\n', upper(mfilename), outputPath)
        varargout = {potcoffs, dates};
        return
    end

    switch Pcenter
        case 'GFZ'
        case 'CSR'

            switch Rlevel
                case 'RL06'
                    inputFileList = ...
                        ls2cell(fullfile(inputFolder, 'GAC-2*'));
                    Ldata = 180;
            end

        case 'JPL'
    end

    % Initialize
    nMonths = length(inputFileList);
    dates = zeros(1, nMonths);
    potcoffs = nan(nMonths, addmup(Ldata), 4);
    cal_errors = nan(nMonths, addmup(Ldata), 4);

    % Loop over the months
    parfor iMonth = 1:nMonths
        % load geopotential coefficients
        inputPath = fullfile(inputFolder, inputFileList{iMonth});
        [dates(iMonth), potcoffs(iMonth, :, :), cal_errors(iMonth, :, :)] = ...
            gac2plm(inputPath, Pcenter, Ldata, "Unit", units, "DateFormat", 'datenum');

        date = datetime(dates(iMonth), "ConvertFrom", 'datenum', "Format", 'yyyy-MM-dd');
        fprintf('%s processed %s (month %i) \n', upper(mfilename), date, iMonth);
    end

    % Save
    save(outputPath, 'potcoffs', 'cal_errors', 'dates');

    % Collect output
    varargout = {dates, potcoffs, cal_errors};
end

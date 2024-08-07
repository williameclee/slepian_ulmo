%% NATLANTIC
% XY = NATLANTIC(upscale, inclang)
% NATLANTIC(...) % Only makes a plot
%
% Finds the coordinates of some of the worlds' oceans so we can combine
% its localization kernel with those for the missing continents to turn
% into Slepian eigenfunctions for all of the world's oceans.
%
% INPUT:
%
% upscale  0 The standard, default values
%          N Splined values at N times the upscaleolution
% inclang  Inclination angle, determines size of polar caps
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the Atlantic ocean
%
% Last modified by williameclee-at-arizona.edu, Jul 2nd, 2024

function varargout = natlantic(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse the inputs
    [upscale, inclang, buf, moreBuf, forceNew, lonOrigin] = ...
        parsecoastinputs(varargin);
    oceanTitle = 'North Atlantic Ocean';
    oceanParts = 'North Atlantic Ocean';

    %% Check if the data file exists
    [dataFile, ~, dataFileExists] = coastfilename(mfilename, ...
        'Upscale', upscale, 'Inclang', inclang, ...
        'Buffer', buf, 'MoreBuffer', moreBuf);

    if dataFileExists && ~forceNew
        load(dataFile, 'XY')

        varargout = returncoastoutput(nargout, XY, oceanTitle);

        return
    end

    %% Compute the data
    [oceanPoly, oceanLatlim, oceanLonlim] = ...
        findoceanboundary(oceanParts, inclang, lonOrigin);

    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LongitudeOrigin', lonOrigin);

    coastPoly = buffer4oceans(coastPoly, buf, ...
        'MoreBuffer', moreBuf, 'LongitudeOrigin', lonOrigin);

    coastPoly = manualadjustment(coastPoly, buf, lonOrigin);

    oceanPoly = subtract(oceanPoly, coastPoly);

    XY = closecoastline(oceanPoly.Vertices);

    % Convert the data to the right format
    [X, Y] = poly2cw(XY(:, 1), XY(:, 2));
    % [Y, X] = flatearthpoly(Y, X, lonOrigin);
    XY = removeduplicatevertices([X, Y]);

    %% Save and return required data
    save(dataFile, 'XY', '-v7.3')
    fprintf('%s saving %s\n', upper(mfilename), dataFile)

    varargout = returncoastoutput(nargout, XY, oceanTitle);

end

%% Subfunctions
function coastPoly = manualadjustment(coastPoly, buf, lonOrigin)

    if buf >= 0.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [45, 50], 'Lonlim', [-68, -59.5], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [48, 51], 'Lonlim', [-60, -58], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [52, 55], 'Lonlim', [-6, -2], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [49.5, 50.5], 'Lonlim', [-2, 1], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [57, 59], 'Lonlim', [8, 12], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 1
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [18, 19], 'Lonlim', [-76, -75], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 2
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [54, 58], 'Lonlim', [0, 6], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 2.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [60, 62], 'Lonlim', [-2, 0], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 3
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [20, 30], 'Lonlim', [-95, -85], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [10, 20], 'Lonlim', [-85, -68], ...
            'LongitudeOrigin', lonOrigin);
    end
    if buf >= 3.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [13, 15], 'Lonlim', [-80, -65], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [60, 65], 'Lonlim', [-10, 0], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 4
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [60, 61], 'Lonlim', [-56, -52], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [56, 60], 'Lonlim', [-12, -10], ...
            'LongitudeOrigin', lonOrigin);
    end
    if buf >= 4.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [0, 2], 'Lonlim', [-8, 5], ...
            'LongitudeOrigin', lonOrigin);
    end

end

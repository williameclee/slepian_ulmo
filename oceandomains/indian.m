%% INDIAN
% XY = INDIAN(upscale, inclang)
% INDIAN(...) % Only makes a plot
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

function varargout = indian(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse the inputs
    [upscale, inclang, buf, moreBuf, forceNew, lonOrigin] = ...
        parsecoastinputs(varargin);
    oceanTitle = 'Indian Ocean';
    oceanParts = 'Indian Ocean';

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

    % Find the coastline
    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LongitudeOrigin', lonOrigin);

    coastPoly = buffer4oceans(coastPoly, buf, ...
        'MoreBuffer', moreBuf, 'LongitudeOrigin', lonOrigin);

    coastPoly = manualadjustment(coastPoly, buf, lonOrigin);

    % Subtract the coast and earthquake mask from the ocean
    oceanPoly = subtract(oceanPoly, coastPoly);

    XY = closecoastline(oceanPoly.Vertices);

    % Think twice before upscaling!
    % Bezier splines may lead to holes between oceans
    if upscale ~= 0 && upscale ~= 1
        XY = bezier(XY, upscale);
    end

    % Convert the data to the right format
    [X, Y] = poly2cw(XY(:, 1), XY(:, 2));
    % [Y, X] = flatearthpoly(Y, X, lonOrigin);
    XY = removeduplicatevertices([X, Y]);

    %% Save and return required data
    save(dataFile, 'XY', '-v7.3')
    fprintf('%s saving %s\n', upper(mfilename), dataFile)

    varargout = returncoastoutput(nargout, XY, oceanTitle);

end

function coastPoly = manualadjustment(coastPoly, buf, lonOrigin)

    coastPoly = addlandregion(coastPoly, ...
        'Latlim', [-16, -15], 'Lonlim', [127, 129], ...
        'LongitudeOrigin', lonOrigin);

    if buf >= 0.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [13, 28], 'Lonlim', [34, 45], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [24, 30], 'Lonlim', [45, 56], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [10, 13], 'Lonlim', [43, 45.5], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [1, 3.5], 'Lonlim', [100, 103], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 1
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [10, 14], 'Lonlim', [45, 48], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-9.5, 8], 'Lonlim', [115, 116], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [4, 8], 'Lonlim', [98, 100], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 1.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [6, 8], 'Lonlim', [97, 99], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 3.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-24, -20], 'Lonlim', [34, 40], ...
            'LongitudeOrigin', lonOrigin);
    end

end

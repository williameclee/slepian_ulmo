%% SATLANTIC
% XY = SATLANTIC(upscale, inclang)
% SATLANTIC(...) % Only makes a plot
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

function varargout = satlantic(varargin)
    %% Initialisation
    % Suppress warnings
    warning('off', 'MATLAB:polyshape:repairedBySimplify');
    % Parse the inputs
    [upscale, inclang, buf, moreBuf, forceNew, lonOrigin] = ...
        parsecoastinputs(varargin, 'DefaultInclang', 90);
    oceanTitle = 'South Atlantic Ocean';
    oceanParts = 'South Atlantic Ocean';

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

    oceanLatlim = oceanLatlim + [-1, 1];
    oceanLatlim = min(max(oceanLatlim, -90), 90);
    oceanLonlim = oceanLonlim + [-1, 1];

    [~, coastPoly] = gshhscoastline('l', 'Buffer', buf, ...
        'LatLim', oceanLatlim, 'LonLim', oceanLonlim, ...
        'LongitudeOrigin', lonOrigin);

    coastPoly = buffer4oceans(coastPoly, buf, ...
        'MoreBuffer', moreBuf, 'LongitudeOrigin', lonOrigin);

    coastPoly = manualadjustment(coastPoly, buf, lonOrigin);

    oceanPoly = subtract(oceanPoly, coastPoly);
    oceanPoly = subtract(oceanPoly, intersect(oceanPoly, coastPoly));
    oceanPoly = subtract(oceanPoly, intersect(oceanPoly, coastPoly));

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

    if buf >= 0
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-3, 0], 'Lonlim', [-52, -48.5], ...
            'LongitudeOrigin', lonOrigin);
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-2, -1], 'Lonlim', [-49, -48], ...
            'LongitudeOrigin', lonOrigin);
    end
    if buf >= 0.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-42, -41], 'Lonlim', [-65, -64], ...
            'LongitudeOrigin', lonOrigin);
    end

    if buf >= 4.5
        coastPoly = addlandregion(coastPoly, ...
            'Latlim', [-61, -58], 'Lonlim', [-68, -61], ...
            'LongitudeOrigin', lonOrigin);
    end

end

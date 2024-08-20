%% SLEPT2RESID_FITSLEPT
% This function is a helper function for SLEPT2RESID.
% This allows us to fit the Slepian coefficients to a polynomial in parallel

% Last modified by
%   williameclee-at-arizona.edu, 7/16/2024

function varargout = slept2resid_fitslept(slept, fitwhat, errors, ...
        isSpTerm, freq, phase, spFreq, spPhase, dateNml, dateExtraNml, ...
        moredates, G1, G2, G3, GSpec1, GSpec2, GSpec3, nMonth)
    P2ftest = 0;
    P3ftest = 0;
    % If we have a priori error information, create a weighting matrix,
    % and change the G and d matrices to reflect this. Since each
    % coefficient has its own weighting, we have to invert them
    % separately.
    weight = diag(1 ./ errors);
    data = slept;
    G1wt = weight * G1;
    G2wt = weight * G2;
    G3wt = weight * G3;
    dataWt = weight * data;

    extra = zeros([length(dateExtraNml), 1]);

    % This is in case you request a single special term to be looked at
    % in a special way
    if isSpTerm
        GSp1wt = weight * GSpec1;
        GSp2wt = weight * GSpec2;
        GSp3wt = weight * GSpec3;
        nOmega = length(spFreq);
        theta = spPhase;
        omega = spFreq;
    else
        nOmega = length(fitwhat(2:end));

        if nOmega > 0
            omega = freq;
            theta = phase;
        end

    end

    %% First order polynomial
    % Do the fitting by minimizing least squares
    if ~isSpTerm
        % The linear fit with periodics
        coeff1 = (G1wt' * G1wt) \ (G1wt' * dataWt);
    else
        % If there was a special one, substitute
        coeff1 = (GSp1wt' * GSp1wt) \ (GSp1wt' * dataWt);
    end

    data1fit = coeff1(1) + coeff1(2) * dateNml;

    % Add the periodic components (if they exist)
    if nOmega > 0
        startP = length(coeff1) - 2 * nOmega + 1;
        % The amplitude of the periodic components
        amp1 = [coeff1(startP:startP + nOmega - 1), ...
                    coeff1(startP + nOmega:end)];
        amp1 = sqrt(amp1(:, 1) .^ 2 + amp1(:, 2) .^ 2);

        % The phase in time of the periodic components
        phase1 = [coeff1(startP:startP + nOmega - 1), ...
                      coeff1(startP + nOmega:end)];
        phase1 = atan2(phase1(:, 1), phase1(:, 2));

        % Adding the periodic components to the fit
        data1fit = data1fit ...
            + sum(repmat(amp1, 1, nMonth) ...
            .* sin(theta' + repmat(phase1, 1, nMonth)), 1);
    end

    % Compute the residual time series for this coefficient
    resid1 = data - data1fit';

    % Here's the definition of the signal at lm vs time
    sleptSig = data1fit;
    % Here's the definition of the residual at lm vs time
    sleptRes = resid1;

    % Do extra dates if you have them
    if moredates
        extraFn1 = coeff1(1) + coeff1(2) * (dateExtraNml);

        if nOmega > 0
            th_extras = repmat(omega, length(dateExtraNml), 1) ...
                * 2 * pi .* repmat((dateExtraNml)', 1, nOmega);
            % Evaluate at the missing dates
            extraFn1 = extraFn1 ...
                + sum(repmat(amp1, 1, 1) .* sin(th_extras' + repmat(phase1, 1, 1)), 1);
        end

        extra = extraFn1;
    end

    % Get the residual sum of squares for later F tests
    rss1 = sum(resid1 .^ 2);

    if fitwhat(1) == 1
        ftests = [0, P2ftest, P3ftest];
        varargout = {sleptSig, sleptRes, ftests, extra};
        return
    end

    %% Second order polynomial
    % Now repeat that all again with second order polynomial, if we want

    if ~isSpTerm
        coeff2 = (G2wt' * G2wt) ...
            \ (G2wt' * dataWt);
    else
        coeff2 = (GSp2wt' * GSp2wt) ...
            \ (GSp2wt' * dataWt);
    end

    data2fit = coeff2(1) + coeff2(2) * (dateNml) ...
        + coeff2(3) * (dateNml) .^ 2;

    if nOmega > 0
        startP = length(coeff2) - 2 * nOmega + 1;
        amp2 = [coeff2(startP:startP + nOmega - 1), ...
                    coeff2(startP + nOmega:end)];
        amp2 = sqrt(amp2(:, 1) .^ 2 + amp2(:, 2) .^ 2);

        phase2 = [coeff2(startP:startP + nOmega - 1), ...
                      coeff2(startP + nOmega:end)];
        phase2 = atan2(phase2(:, 1), phase2(:, 2));

        data2fit = data2fit + ...
            sum(repmat(amp2, 1, nMonth) ...
            .* sin(theta' + repmat(phase2, 1, nMonth)), 1);
    end

    % Compute the residual time series for this coefficient
    resid2 = data - data2fit';

    % Do extra dates if you have them
    if moredates
        extraFn2 = coeff2(1) + coeff2(2) * (dateExtraNml) ...
            + coeff2(3) * (dateExtraNml) .^ 2;

        if nOmega > 0
            th_extras = repmat(omega, length(dateExtraNml), 1) ...
                * 2 * pi .* repmat((dateExtraNml)', 1, nOmega);
            extraFn2 = extraFn2 ...
                + sum(repmat(amp1, 1, 1) ...
                .* sin(th_extras' + repmat(phase1, 1, 1)), 1);
        end

    end

    % Get the residual sum of squares
    rss2 = sum(resid2 .^ 2);
    % Calculate an F-score for this new fit
    fratioP2 = (rss1 - rss2) / 1 / (rss2 / (length(slept) - length(coeff2)));
    fscore = finv(.95, 1, length(slept) - length(coeff2));

    if fratioP2 > fscore
        P2ftest = 1;
        % We pass, so update the signal and residuals with this new fit
        sleptRes = resid2;
        sleptSig = data2fit;

        if moredates
            extra = extraFn2;
        end

    else
        P2ftest = 0;
    end

    if fitwhat(1) == 2
        ftests = [0, P2ftest, P3ftest];
        varargout = {sleptSig, sleptRes, ftests, extra};
        return
    end

    %% Third order polynomial
    % Now repeat that all again with third order polynomial, if we want
    if ~isSpTerm
        coeff3 = (G3wt' * G3wt) \ (G3wt' * dataWt);
    else
        coeff3 = (GSp3wt' * GSp3wt) \ (GSp3wt' * dataWt);
    end

    data3fit = coeff3(1) + coeff3(2) * (dateNml) ...
        + coeff3(3) * (dateNml) .^ 2 ...
        + coeff3(4) * (dateNml) .^ 3;

    if nOmega > 0
        startP = length(coeff3) - 2 * nOmega + 1;
        amp3 = [coeff3(startP:startP + nOmega - 1), ...
                    coeff3(startP + nOmega:end)];
        amp3 = sqrt(amp3(:, 1) .^ 2 + amp3(:, 2) .^ 2);

        phase3 = [coeff3(startP:startP + nOmega - 1), ...
                      coeff3(startP + nOmega:end)];
        phase3 = atan2(phase3(:, 1), phase3(:, 2));

        data3fit = data3fit + ...
            sum(repmat(amp3, 1, nMonth) ...
            .* sin(theta' + repmat(phase3, 1, nMonth)), 1);
    end

    resid3 = data - data3fit';

    % Do extra dates if you have them
    if moredates
        extraFn3 = coeff3(1) + coeff3(2) * (dateExtraNml) ...
            + coeff3(3) * (dateExtraNml) .^ 2 ...
            + coeff3(4) * (dateExtraNml) .^ 3;

        if nOmega > 0
            th_extras = repmat(omega, length(dateExtraNml), 1) ...
                * 2 * pi .* repmat((dateExtraNml)', 1, nOmega);
            extraFn3 = extraFn3 ...
                + sum(repmat(amp1, 1, 1) ...
                .* sin(th_extras' + repmat(phase1, 1, 1)), 1);
        end

    end

    % Get the residual sum of squares
    rss3 = sum(resid3 .^ 2);
    % Calculate an F-score for this new fit
    fratioP3 = (rss1 - rss3) / 2 / (rss3 / (length(slept) - length(coeff3)));
    fscore = finv(.95, 1, length(slept) - length(coeff3));

    if fratioP3 > fscore
        P3ftest = 1;
        % We pass, so update the signal and residuals with this new fit
        sleptRes = resid3;
        sleptSig = data3fit';

        if moredates
            extra = extraFn3';
        end

    else
        P3ftest = 0;
    end

    % Make the matrix ftests
    ftests = [0, P2ftest, P3ftest];
    varargout = {sleptSig, sleptRes, ftests, extra};
end

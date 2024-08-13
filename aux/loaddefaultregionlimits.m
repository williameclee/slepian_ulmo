function [latLim, lonLim] = loaddefaultregionlimits(region)

    switch region
        case 'africa'
            latLim = [-50, 50];
            lonLim = [-30, 60];
        case 'namerica'
            latLim = [10, 90];
            lonLim = [-170, -50];
        case 'samerica'
            latLim = [-65, 20];
            lonLim = [-90, -30];
        case 'antarctica'
            latLim = [-90, -60];
            lonLim = [-180, 180];
        case 'australia'
            latLim = [-50, -5];
            lonLim = [105, 165];
        case 'eurasia'
            latLim = [5, 90];
            lonLim = [-15, 190];
        case 'greenland'
            latLim = [56, 88];
            lonLim = [-91, 10];
        case 'oceans'
            latLim = [-90, 90];
            lonLim = [0, 360];
        case 'atlantic'
            latLim = [-90, 75];
            lonLim = [-100, 25];
        case 'satlantic'
            latLim = [-80, 5];
            lonLim = [-75, 25];
        case 'natlantic'
            latLim = [-5, 70];
            lonLim = [-100, 15];
        case 'pacific'
            latLim = [-90, 80];
            lonLim = [110, 300];
        case 'spacific'
            latLim = [-90, 5];
            lonLim = [110, 300];
        case 'npacific'
            latLim = [-5, 80];
            lonLim = [110, 290];
        case 'indian'
            latLim = [-80, 35];
            lonLim = [15, 150];
    end

end

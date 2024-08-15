product = {'CSR', 'RL06', 60};

domainNames = {'indian', 'npacific'};
incs = [60, 90];
bufs = [0, 1];
L = 18;

for iDom = 1:length(domainNames)
    domainName = domainNames{iDom};

    for iInc = 1:length(incs)
        inc = incs(iInc);

        for iBuf = 1:length(bufs)
            buf = bufs(iBuf);
            domain = GeoDomain(domainName, "Buffer", buf, "Inclination", inc);
            S0 = grace2slept(product, domain.Lonlat, 0, L, ...
                [], [], [], [], 'SD', 1);
            S = grace2slept_new(product, domain, L, "Unit", 'SD', "ForceNew", true, "SaveData", false);

            assert(isequal(S0, S))
        end

    end

end

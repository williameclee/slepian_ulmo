domainNames = {'arctic', 'npacific'};
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

            K0 = kernelcp(L, domain.Lonlat);

            K = kernelcp_new(L, domain, ...
                "ForceNew", true, "SaveData", false);
            assert(isequal(K0, K))
        end

    end

end

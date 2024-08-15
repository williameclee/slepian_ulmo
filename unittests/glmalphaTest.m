domainNames = {'npacific', 'indian'};
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

            G0 = glmalpha(domain.Lonlat, L);

            G = glmalpha_new(domain, L, ...
                "ForceNew", true, "SaveData", false);

            try
                assert(isequal(G0, G))
            catch
				disp(domainName)
				disp(inc)
				disp(buf)
                assert(isequal(G0, G))
            end

        end

    end

end

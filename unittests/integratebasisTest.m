domain = GeoDomain('pacific');
L = 18;

[G, V] = glmalpha_new(domain, L);

I = integratebasis_new(G, domain, "ForceNew", true);
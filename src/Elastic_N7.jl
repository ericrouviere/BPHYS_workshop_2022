module Elastic_N7

export makeRandStrain, computeEnergy, cooperativity, prepareAllosterySitesStrains, 
    computeFitness


function makeRandStrain(site, r)
    # generate a random strain,
    # rotations and translations and normalize.
    strain = randn(length(site))
    removeTranslationRotation!(strain, r[site])
    normalize!(strain)
    return strain
end

function computeEnergy(net, strain, site)
    E, F, Δr = computeResponse(net.r, net.H, strain, site)
    return E
end

cooperativity(E10, E01, E11) = E10 + E01 - E11

function prepareAllosterySitesStrains(W, L)
    # build sites and strains vectors for computing allostery.
    actSiteNodes = collect(Int(floor(W/2)-1) : Int(floor(W/2))+2)
    allSiteNodes = W*L .- actSiteNodes .+ 1
    actSite = sort([actSiteNodes * 2; actSiteNodes * 2 .- 1])
    allSite = sort([allSiteNodes * 2; allSiteNodes * 2 .- 1])
    actStrain = makeRandStrain(actSite, r)
    allStrain = makeRandStrain(allSite, r)
    sites =  fill([actSite; allSite], 3)
    strains =  [ [actStrain; zeros(8)], [zeros(8); allStrain], [actStrain; allStrain]]
    return sites, strains
end

function computeFitness(net::Network, strains, sites)
    # compute the allosteric fitness for a network.
    E10 = computeEnergy(net, strains[1], sites[2])
    E01 = computeEnergy(net, strains[2], sites[2])
    E11 = computeEnergy(net, strains[3], sites[3])
    fitness = cooperativity(E10, E01, E11)
    return fitness
end

end

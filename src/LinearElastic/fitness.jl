

stiffness(E) = E
softness(E) = -E
softness2(E) = 1/E
logStiffness(E) = log(10, E)
logSoftness(E) = -log(10, E)

specificity(E1, E2) = log(10,E1) - log(10, E2) # select for log gap in stiffnesses



function computeFitness(Q::Network, pokes::Pokes, assay::String)
   
    energies = computeEnergiesFast(Q, pokes)
    if assay == "Stiffness"
        ϕ = stiffness(energies...)
    elseif assay == "Softness"
        ϕ = softness(energies...)
    elseif assay == "Softness2"
        ϕ = softness2(energies...)
    elseif assay == "LogStiffness"
        ϕ = logStiffness(energies...)
    elseif assay == "LogSoftness"
        ϕ = logSoftness(energies...)
    elseif assay == "Specificity"
        ϕ = specificity(energies...)
    else
        error(assay*" : Not supported assay")
    end
    return ϕ
end







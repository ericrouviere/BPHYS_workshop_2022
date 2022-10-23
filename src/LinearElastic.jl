module LinearElastic

using LinearAlgebra, LinearSolve, PyPlot, Statistics, Random, JLD, Distributed


push!(LOAD_PATH,"./")

include("LinearElastic/types.jl")
include("LinearElastic/buildnetwork.jl")
include("LinearElastic/physics.jl")
include("LinearElastic/evolve.jl")
include("LinearElastic/fitness.jl")
include("LinearElastic/mutate.jl")
include("LinearElastic/tools.jl")
include("LinearElastic/analysis.jl")
include("LinearElastic/plotfuncs.jl")
include("LinearElastic/ensemble.jl")

export

    #from types.jl
    Network,
    Pokes,
    NetParams,
    EvoParams,

    # from buildnetwork.jl
    buildStructure,
    buildSpringTable,
    seq2Springs,
    randSeq,
    buildNetwork,
    buildPopulation,
    buildStructAdjMatTable,
    makePokesListFlux,

    # from physics.jl
    buildStructureMatrix,
    computeHessian!,
    computeHessian,
    updateHessian!,
    computeResponse!,
    computeResponse,
    computeResponse_LM,
    getTranslationRotation,
    removeTranslationRotation!,
    computeEnergies,
    computeEnergiesFast,
    createNormalStrains,
    fillBMatrix!,

    # from evolve.jl
    evolveMC,
    evolvePop,
    evolvePop!,
    updateEnvironment,
    samplingDetails,

    # from fitness.jl
    computeFitness,
    stiffness,
    softness,

    # from mutate.jl
    mutateNode!,
    updateSpringMatrix!,
    mutateAtRate!,

    # from tools.jl
    r2xy,
    xy2r,
    networkCopy!,
    pop2Seqs,
    convertSettings2NetParams,
    convertSettings2EvoParams,
    compressPokesList,
    expandPokesList,
    saveEvoData,

    # from analysis.jl
    computeDMS,
    countMutSensitivity,
    computeEnergyDMS,

    # from plotfuncs.jl
    plotNetwork!,
    plotNetwork,
    plotDisplacment!,
    plotBondStiffness!,
    
    # from ensemble.jl
    evolvePopEnsemble,
    computeEnergyDMSEns,
    countMutSensitivityEns

end

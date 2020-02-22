### execution path
set ExecutionPath {
    ParticlePropagator
    Merger
    Efficiency
    DecayFilter
    TrackSmearing
    TimeSmearing
    TreeWriter
}

### propagate inside the solenoid
module ParticlePropagator ParticlePropagator {
    set InputArray Delphes/stableParticles
    set OutputArray stableParticles
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons
    set MuonOutputArray muons
    set Radius 100.e-2
    set HalfLength 200.e-2
    set Bz 0.2
}

### merge all tracks
module Merger Merger {
    add InputArray ParticlePropagator/chargedHadrons
    add InputArray ParticlePropagator/electrons
    add InputArray ParticlePropagator/muons
    set OutputArray tracks
}

### efficiency and acceptance
module Efficiency Efficiency {
    set InputArray Merger/tracks
    set OutputArray tracks
    set EfficiencyFormula { 1.0 }
}

### decays
module DecayFilter DecayFilter {
    set InputArray Efficiency/tracks
    set OutputArray tracks
}

### tracking resolution
module TrackSmearing TrackSmearing {
    set InputArray DecayFilter/tracks
    set OutputArray tracks
    set PResolutionFormula { 0.01 }
    set CtgThetaResolutionFormula { 0.0 }
    set PhiResolutionFormula { 0.0 }
    set D0ResolutionFormula { 0.0 }
    set DZResolutionFormula { 0.0 }
}

### time resolution
module TimeSmearing TimeSmearing {
    set InputArray TrackSmearing/tracks
    set OutputArray tracks
    set TimeResolution 0.020e-9
}

### tree definition
module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass
    add Branch Delphes/allParticles particles GenParticle
    add Branch TimeSmearing/tracks tracks Track
}

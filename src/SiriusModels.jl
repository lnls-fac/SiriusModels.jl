"""Models module."""
module SiriusModels
    include("SI/StorageRing.jl")
    export StorageRing

    using Track.Auxiliary
    using PrecompileTools

    @setup_workload begin
        @compile_workload begin
            m = StorageRing.create_accelerator()
            m.radiation_state = Auxiliary.full
            m.cavity_state = Auxiliary.on
            m.vchamber_state = Auxiliary.on
        end
    end

end # module SiriusModels

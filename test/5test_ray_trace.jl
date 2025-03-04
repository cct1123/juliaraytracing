using Test
include("../MathModule.jl")

include("../GeometryModule.jl")
include("../SourceModule.jl")
include("../MaterialModule.jl")

include("../ObjectModule.jl")
include("../DetectorModule.jl")
include("../IllustrationModule.jl")

include("../RayTrace.jl")
using .RayTrace

using Main.GeometryModule: Volume, Surface

@testset "RayTrace Tests" begin
    # Test ray-surface intersection
    @testset "Surface Intersection" begin
        nothing
    end

    # Test ray-volume interaction
    @testset "Volume Interaction" begin
        nothing
    end

    # Test ray-object intersection
    @testset "Object Intersection" begin
        nothing
    end

    # Test dielectric material scattering
    @testset "Dielectric Scattering" begin
        nothing
    end
end

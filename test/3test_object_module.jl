
using Test
include("../MathModule.jl")
include("../GeometryModule.jl")
include("../MaterialModule.jl")
include("../ObjectModule.jl")
using .ObjectModule

# something modules for test
using .GeometryModule: Surface, Volume, make_box, make_plane
using .MaterialModule: Material, Dielectric

@testset "ObjectModule Tests" begin
    # Create dummy instances of Volume, Surface, and Material
    volume = make_box(10.0, 10.0, 10.0)  
    surface = make_plane(10.0, 10.0)
    material = Dielectric(1.0)

    # Test Object creation with Volume
    obj1 = Object(volume, material)
    @test obj1.geometry isa Volume
    @test obj1.material isa Material

    # Test Object creation with Surface
    obj2 = Object(surface, material)
    @test obj2.geometry isa Surface
    @test obj2.material isa Material
end

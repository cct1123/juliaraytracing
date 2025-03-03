using Test
include("../MaterialModule.jl")
using .MaterialModule

@testset "MaterialModule Tests" begin
    # Test Dielectric material creation
    d1 = Dielectric(1.5, 0.01)
    @test d1 isa Dielectric
    @test d1.n == 1.5
    @test d1.α == 0.01
    
    # Test default Dielectric materials
    @test air.n == 1.0
    @test diamond.n == 2.4
    @test diamond_ib_heavy.α == 50.0E-4
    @test diamond_ib_mild.α == 10.0E-4
    @test diamond_iia.α == 0.1E-4
    @test pdms.n == 1.4
    
    # Test Metal material creation
    m1 = Metal(5.8E7)
    @test m1 isa Metal
    @test m1.σ == 5.8E7
    
    # Test SheetMetal material creation
    sm1 = SheetMetal(5.8E7, 1.0E-6, 0.01)
    @test sm1 isa SheetMetal
    @test sm1.σ == 5.8E7
    @test sm1.δ == 1.0E-6
    @test sm1.thickness == 0.01
end

using Test
using LinearAlgebra


include("../MathModule.jl")  
using .MathModule  

include("../SourceModule.jl")
using .SourceModule

# Test isotropic distribution
@testset "Isotropic Distribution" begin
    dist = isotropic_distribution()  # Get isotropic distribution function
    dir1 = dist()  # Sample direction
    dir2 = dist()  # Sample another direction
    @test typeof(dir1) == Vector{Float64}
    @test length(dir1) == 3
    @test dir1 != dir2  # The directions should be different
end

# Test collimated distribution
@testset "Collimated Distribution" begin
    test_direction = [1.0, 0.0, 0.0]  # Define a fixed direction
    dist = collimated_distribution(test_direction)
    dir1 = dist()  # Sample direction (should always be the same)
    @test dir1 == normalize(test_direction)
end

# Test cone distribution
@testset "Cone Distribution" begin
    main_direction = [0.0, 0.0, 1.0]  # Direction pointing in z-axis
    spread_angle = 0.1  # Small spread angle
    dist = cone_distribution(main_direction, spread_angle)
    dir1 = dist()  # Sample direction
    @test typeof(dir1) == Vector{Float64}
    @test length(dir1) == 3
    @test abs(acos(dot(dir1, main_direction)) - spread_angle) < 0.1  # Direction should be within spread angle
end

# Test dipole distribution
@testset "Dipole Distribution" begin
    dipole_moment = [0.0, 0.0, 1.0]  # Dipole pointing in z-direction
    dist = dipole_distribution(dipole_moment)
    dir1 = dist()  # Sample direction
    @test typeof(dir1) == Vector{Float64}
    @test length(dir1) == 3
    @test abs(dot(dir1, dipole_moment)) < cos(π/6)  # Should be perpendicular to dipole moment
end

# Test PointSource structure
@testset "PointSource Structure" begin
    frame = Frame()  # Assuming Frame has a default constructor
    wavelength = 500.0  # Wavelength in micrometers
    intensity = 1.0  # Intensity of light source
    distribution = isotropic_distribution()  # Isotropic distribution function
    point_source = PointSource(frame, wavelength, intensity, distribution)
    
    @test typeof(point_source) == PointSource
    @test point_source.frame == frame
    @test point_source.wavelength == wavelength
    @test point_source.intensity == intensity
    @test typeof(point_source.distribution()) == Vector{Float64}
end

# Test EnsembleSource structure
@testset "EnsembleSource Structure" begin
    frame = Frame()  # Assuming Frame has a default constructor
    wavelength = 500.0  # Wavelength in micrometers
    intensity = 1.0  # Intensity of light source
    distribution = isotropic_distribution()  # Isotropic distribution function
    point_source = PointSource(frame, wavelength, intensity, distribution)
    
    ensemble = EnsembleSource([point_source, point_source])
    @test typeof(ensemble) == EnsembleSource
    @test length(ensemble.objects) == 2
    @test ensemble.objects[1] == point_source
end
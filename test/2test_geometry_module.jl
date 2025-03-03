# Test script for GeometryModule

using Test

include("../MathModule.jl")   
include("../GeometryModule.jl")
using .GeometryModule

using .MathModule: Frame

# Test data setup
const center = [0.0, 0.0, 0.0]
const radius = 1.0
const height = 2.0
const lx = 2.0
const ly = 2.0
const lz = 2.0

# Test surface shapes
@testset "Surface Shapes" begin
    # Test sphere shape
    shape_sp = shape_sphere(center, radius)
    @test shape_sp([0.0, 0.0, 1.0]) ≈ 0.0   # Point on surface

    # Test hemisphere shape
    shape_hem = shape_hemisphere(center, radius)
    @test shape_hem([0.0, 0.0, 1.0]) ≈ 0.0   # Point on surface

    # Test paraboloid shape
    shape_parab = shape_paraboloid(center, 1.0)
    @test shape_parab([1.0, 0.0, 1.0]) ≈ 0.75   # Point on surface
end

# Test volume creation functions
@testset "Volume Creation" begin
    # Test making a box
    box = make_box(lx, ly, lz)
    @test length(box.surfaces) == 6  # Ensure 6 surfaces for the box

    # Test making a cone
    cone = make_cone(1.0, 2.0, height)
    @test length(cone.surfaces) == 2  # Ensure 2 surfaces for the cone

    # Test making a sphere
    sphere = make_sphere(radius)
    @test length(sphere.surfaces) == 2  # Ensure 2 surfaces for the sphere

    # Test making a cylinder
    cylinder = make_cylinder(1.0, height)
    @test length(cylinder.surfaces) == 3  # Ensure 3 surfaces for the cylinder (top, bottom and shell)
end

# Test geometric transformations
@testset "Geometric Transformations" begin
    # Define a box and apply rotation
    box = make_box(lx, ly, lz)
    rotated_box = BoundVolume(Frame(center, rm_xflip), box.surfaces, box.bounds)
    
    @test rotated_box.frame.orientation == rm_xflip  # Check if the rotation is applied
end

println("All tests passed!")

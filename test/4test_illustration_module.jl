using Test
include("../MathModule.jl")
include("../GeometryModule.jl")
include("../SourceModule.jl")
include("../DetectorModule.jl")
include("../IllustrationModule.jl")
using .IllustrationModule

using .MathModule: Frame, rotation_matrix
using .SourceModule: Ray
using .GeometryModule: ImplicitSurface, Surface, shape_paraboloid, shape_mexican_hat, border_circle, border_bounds, Volume, make_cylinder, surface_normal
using .DetectorModule: TrajectoryNode, Trajectory

# Test Ray creation and ray drawing
@testset "Ray Drawing Test" begin
    # Create a few sample rays
    ray1 = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
    ray2 = Ray([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    
    # Draw rays and ensure traces are returned
    rays = [ray1, ray2]
    traces = draw_rays(rays)
    @test length(traces) == 3
end

# Test Surface drawing with default resolution
@testset "Surface Drawing Test" begin
    # Mock Surface object for testing
    focal_length = 2.0  # Focal length of the paraboloid
    bound_radius = 4.5
    identity_matrix = rotation_matrix([0.0, 0.0, 1.0], 0.0) 
    labframe = Frame([0.0,10.0,0.0], identity_matrix)
    x_range = (-5.0, 5.0)  # X range
    y_range = (-5.0, 5.0)  # Y range
    bounds = [x_range, y_range]
    surface = ImplicitSurface(
        labframe, 
        shape_paraboloid([0.0, 0.0, 0.0], focal_length), 
        bounds,
        border_circle(bound_radius))

    # Test surface drawing
    surface_traces = draw_surface(surface, 50, 50; color="random", display=false)
    @test length(surface_traces) > 0
end

# Test Normals drawing for a surface
@testset "Surface Normal Drawing Test" begin
    flip_matrix = rotation_matrix([1.0, 0.0, 0.0], π*1.0) 
    objframe3 = Frame([0.0,-5.0,0.0], flip_matrix)
    x_range = (-5.0, 5.0)  # X range
    y_range = (-5.0, 5.0)  # Y range
    bounds = [x_range, y_range]
    surface = ImplicitSurface(
        objframe3, 
        shape_mexican_hat([0.0, 0.0, 10.0], 4.0, 1.0, 1.0),
        bounds,
        border_bounds()
        )

    # Test normal drawing
    normal_traces = draw_snormals(surface, 50; ray_length=0.5, arrow_scale=0.3, color="random")
    @test length(normal_traces) > 0
end

# Test Volume Drawing
@testset "Volume Drawing Test" begin
    volume = make_cylinder(10.0, 50.0)
    volume_traces = draw_volume(volume, 5, 5, 5; display=false)
    @test length(volume_traces) > 0
end

# Test Trajectory Drawing
@testset "draw_trajectory Test" begin
    # Define a mock trajectory tree
    root_ray = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0], amplitude=1.0)
    child_ray1 = Ray([1.0, 1.0, 0.0], [1.0, -1.0, 0.0], amplitude=0.7)
    child_ray2 = Ray([1.0, 1.0, 0.0], [-1.0, 1.0, 0.0], amplitude=0.5)
    child_ray11 = Ray([3.0, -2.0, 0.0], [0.0, 1.0, 1.0], amplitude=0.5)

    
    root_node = TrajectoryNode(root_ray, [])
    child_node2 = TrajectoryNode(child_ray2, [])
    child_node1 = TrajectoryNode(child_ray1, [])
    child_node11 = TrajectoryNode(child_ray11, [])
    
    push!(root_node.children, child_node1, child_node2)
    push!(child_node1.children, child_node11)
    traj = Trajectory(root_node)

    # Call the function
    traces = draw_trajectory(traj; color="rgba(255,0,0, 1)", display=false)

    # Verify output
    @test length(traces) == 4 
end
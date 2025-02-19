using LinearAlgebra
using PlotlyJS
using ForwardDiff, Roots
using Random  # To generate random colors

include("math.jl")

const Vec3 = Vector{Float64};
const Vec2 = Vector{Float64};
const Point = Vec3;

# Define the Ray structure ===============================================================================
mutable struct Ray
    origin::Point # 3d vector origin position
    direction::Vec3 # 3d vector direction
    # polarization:: Vec2 # 2d vector polariation looking into the propagation direction, s-wave , p-wave components
end
# ==============================================================================================================

# Define the structure for light sources ===============================================================================
function isotropic_distribution()
    θ = 2π * rand()     # Azimuthal angle (0 to 2π)
    ϕ = acos(2 * rand() - 1)  # Polar angle (cosine-weighted for uniform sphere)
    return [sin(ϕ) * cos(θ), sin(ϕ) * sin(θ), cos(ϕ)]
end

function collimated_distribution(direction::Vec3)
    return () -> normalize(direction)  # Always emits in the same direction
end

function cone_distribution(main_direction::Vec3, spread_angle::Float64)
    function sample_direction()
        θ = spread_angle * sqrt(rand())  # Cone angle (scaled probability)
        ϕ = 2π * rand()  # Uniform azimuth

        # Local coordinate system
        z = normalize(main_direction)
        x = normalize(cross(z, [0, 1, 0]))  # Ensure perpendicular x
        y = cross(z, x)  # Ensure perpendicular y

        # Convert to Cartesian
        dir = cos(θ) * z + sin(θ) * (cos(ϕ) * x + sin(ϕ) * y)
        return normalize(dir)
    end
    return sample_direction
end

custom_dist = () -> normalize(randn(3))  # Gaussian-distributed emission

abstract type Source end

struct PointSource <: Source
    # define the point source structure
    center::Point # 3d vector center position
    intensity::Float64 # intensity of the light source
    distribution::Function  # Function that returns a direction
end

struct EnsembleSource <: Source
    objects::Vector{Source}
end

function emit_ray(source::Source)::Ray
    direction = source.distribution()
    return Ray(source.center, direction)
end

function emit_rays(source::Source, num_rays::Int)::Vector{Ray}
    rays = Ray[]
    for _ in 1:num_rays
        push!(rays, emit_ray(source))
    end
    return rays
end
# ==============================================================================================================

# Define the structure for surfaces and objects ===============================================================================
# Abstract Object for the rays to hit 

struct  Frame 
    origin:: Point
    orientation:: Matrix{Float64} # a 3x3 rotation matrix that transform the lab frame to the object frame
    
end

struct Surface
    # define the surface structure
    frame:: Frame # the frame of the surface
    shape:: Function # bound functions specifying the 3d shape
    bounds:: Vector{Tuple{Float64, Float64}}
    border:: Function # bound functions specifying the boundary line
end

# example shape function: Define a spherical surface as the shape function
function shape_sphere(center::Vec3, radius::Float64)
    return (p) -> (p - center)⋅(p - center) - radius^2  # Sphere equation: ||p - center|| = radius
end

function shape_hemisphere(center::Vec3, radius::Float64)
    function result(p)
        incap = radius^2-(p[1] - center[1])^2-(p[2] - center[2])^2
        return incap < 0 ? NaN : (p[3]- center[3]) - sqrt(incap) 
    end
    return (p) -> result(p)
end

# Define a paraboloid shape function
# function shape_paraboloid(center::Vec3, f::Float64)
#     return (p::Vec3) -> (p[3]- center[3]) - ((p[1]- center[1])^2 + (p[2]- center[2])^2) / (4 * f)  # Paraboloid equation: z = (x^2 + y^2) / (4f)
# end

function shape_paraboloid(center::Vec3, f::Float64)
    return (p) -> (p[3]- center[3]) - ((p[1]- center[1])^2 + (p[2]- center[2])^2) / (4 * f)  # Paraboloid equation: z = (x^2 + y^2) / (4f)
end

function shape_mexican_hat(center::Vec3, A::Float64, B::Float64, C::Float64)
    return (p) -> (p[3] - center[3]) - A * (1 - ((p[1] - center[1])^2 + (p[2] - center[2])^2) / B^2) * exp(-((p[1] - center[1])^2 + (p[2] - center[2])^2) / (2 * C^2))
end


function bound_circle(center::Vec2, radius::Float64)
    return (p) -> (p - center)⋅(p - center) - radius^2  # Circle equation: ||p - center|| = radius
end

function bound_none(center::Vec2, radius::Float64)
    return (p) -> -1  # no bound line defined
end

function surface_normal(p, surface::Surface)

    # Calculate the gradient of the surface shape function
    grad_S = ForwardDiff.gradient(surface.shape, p)
    
    # Normalize the gradient to get the normal
    return grad_S ./ norm(grad_S)
end

function draw_surface(surface:: Surface, nx::Int=100, ny::Int=100, color="random")
    # plot the surface
    xarray = range(surface.bounds[1][1], stop=surface.bounds[1][2], length=nx)
    yarray = range(surface.bounds[2][1], stop=surface.bounds[2][2], length=ny)
    xmesh = zeros(nx, ny)
    ymesh = zeros(nx, ny)
    zmesh = zeros(nx, ny)
    for ii in 1:1:nx, jj in 1:1:ny
    
        bounded = surface.border([xarray[ii], yarray[jj]]) < 0
        xval = xarray[ii]
        yval = yarray[jj]
        zval = bounded ? -surface.shape([xarray[ii], yarray[jj], 0.0]) : NaN

        xval, yval, zval = surface.frame.orientation * [xval, yval, zval]
        xmesh[ii, jj] = xval + surface.frame.origin[1]
        ymesh[ii, jj] = yval + surface.frame.origin[2]
        zmesh[ii, jj] = zval + surface.frame.origin[3]
    end
    # Generate a random color (RGB format)
    if color=="random"
        color = string("rgb(", round(255 * rand()), ",", round(255 * rand()), ",", round(255 *  rand()), ")")
    end;

    surface_data = PlotlyJS.surface(
        x=xmesh,
        y=ymesh,
        z=zmesh,
        surfacecolor=zmesh*0.0,
        colorscale=[(0, color), (1, color)],  # Fix the color scale to the same random color
        shading="smooth",  # Smoothing the shading on the surface
        showscale=false,
        opacity=0.8
        )
    return [surface_data]
end

function draw_snormals(surface:: Surface, num:: Int, color="random")
    xarray = range(surface.bounds[1][1], stop=surface.bounds[1][2], length=nx)
    yarray = range(surface.bounds[2][1], stop=surface.bounds[2][2], length=ny)

    rays_snorm = Ray[]
    for ii in 1:1:num
        bounded = false
        while ~bounded
            xx = rand(surface.bounds[1][1]:surface.bounds[1][2])
            yy = rand(surface.bounds[2][1]:surface.bounds[2][2])
            bounded = surface.border([xx, yy]) < 0
        end
        zz = -surface.shape([xx,yy,0.0])
        snorm = surface_normal([xx, yy, zz], surface)
        snorm = surface.frame.orientation * snorm
        xx, yy, zz = surface.frame.orientation * [xx, yy, zz] + surface.frame.origin
        push!(rays_snorm, Ray([xx, yy, zz], snorm))
    end
    return draw_rays(rays_snorm)
end

# Function to generate and plot rays with arrows
function draw_rays(rays::Vector{Ray}; ray_length=0.5, arrow_scale=0.3)
    traces = GenericTrace[]

    # Line traces for the rays
    for ray in rays
        start = ray.origin
        stop = start + normalize(ray.direction) * ray_length
        trace_rays = scatter3d(
            x=[start[1], stop[1]], y=[start[2], stop[2]], z=[start[3], stop[3]],
            mode="lines",
            line=attr(width=3),
        )
        push!(traces, trace_rays)
    end

    # Arrow traces using cones
    origins = hcat([ray.origin + normalize(ray.direction) * ray_length for ray in rays]...)  # Arrow positions
    directions = hcat([normalize(ray.direction) for ray in rays]...)  # Arrow directions

    trace_arrows = cone(
        x=origins[1, :], y=origins[2, :], z=origins[3, :],  # Positions
        u=directions[1, :], v=directions[2, :], w=directions[3, :],  # Directions
        sizemode="absolute",
        sizeref=arrow_scale,  # Controls arrow size
        anchor="tail",
        showscale=false,
        opacity=0.8
    )
    push!(traces, trace_arrows)

    return traces
end
# ==============================================================================================================


if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    # example usage of the drawing functions
    focal_length = 2.0  # Focal length of the paraboloid
    hem_radius = 4.0  # Radius of the hemispherical boundary
    bound_radius = 4.5
    x_range = (-5.0, 5.0)  # X range
    y_range = (-5.0, 5.0)  # Y range
    nx = 200 # number of points along x
    ny = 200 # number of points along y
    num_sn = 200 # number of surface normals
    identity_matrix = rotation_matrix([0.0, 0.0, 1.0], 0.0) 
    rotx_matrix = rotation_matrix([1.0, 0.0, 0.0], π/4.0) 
    labframe = Frame([0.0,10.0,0.0], identity_matrix)
    objframe = Frame([10.0,0.0,0.0], rotx_matrix)
    bounds = [(x_range[1], x_range[2]), (y_range[1], y_range[2])]

    parasurface = Surface(
        labframe, 
        # objframe,
        # shape_paraboloid([0.0, 0.0, 0.0], focal_length), 
        shape_mexican_hat([0.0, 0.0, 0.0], 4.0, 1.0, 1.0),
        bounds,
        bound_circle([0.0, 0.0], bound_radius))



    hemisurface = Surface(
        objframe, 
        # shape_paraboloid([0.0, 0.0, 0.0], focal_length), 
        shape_hemisphere([0.0, 0.0, 0.0], hem_radius), 
        bounds,
        bound_circle([0.0, 0.0], bound_radius))

    traces1_sf = draw_surface(parasurface, nx, ny)
    traces1_sn = draw_snormals(parasurface, num_sn)
    traces2_sf = draw_surface(hemisurface, nx, ny)
    traces2_sn = draw_snormals(hemisurface, num_sn)

    layout = Layout(
        title="Surface Plot",
        width=800,   # Set width in pixels
        height=600,  # Set height in pixels
        scene=attr(
            xaxis=attr(visible=false, showgrid=false, zeroline=false),
            yaxis=attr(visible=false, showgrid=false, zeroline=false),
            zaxis=attr(visible=false, showgrid=false, zeroline=false),
            bgcolor="rgba(0,0,0,0)",  # Transparent background
            aspectmode="data"
        ),
        paper_bgcolor="rgba(0,0,0,0)",  # Transparent outer background
        showlegend=false
    )
    # Display the plot
    fig = plot(vcat(traces1_sf,traces1_sn, traces2_sf, traces2_sn) , layout)
    display(fig)
end
using LinearAlgebra
using PlotlyJS
using ForwardDiff, Roots
using Random 

include("math.jl")

const Vec3 = Vector{Float64};
const Vec2 = Vector{Float64};
const Point = Vec3;
const invalid_vector =[NaN, NaN, NaN]

# Define the Ray structure ===============================================================================
mutable struct Ray
    origin::Point # 3d vector origin position
    direction::Vec3 # 3d vector direction
    # polarization:: Vec2 # 2d vector polariation looking into the propagation direction, s-wave , p-wave components
end
# ==============================================================================================================

# Define the structure for light sources ===============================================================================
function isotropic_distribution()
    function sample_direction()
        θ = 2π * rand()     # Azimuthal angle (0 to 2π)
        ϕ = acos(2 * rand() - 1)  # Polar angle (cosine-weighted for uniform sphere)
        return [sin(ϕ) * cos(θ), sin(ϕ) * sin(θ), cos(ϕ)]
    end
    return () -> sample_direction()
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

# Define the structure for surfaces  ===============================================================================
struct  Frame 
    origin:: Point
    orientation:: Matrix{Float64} # a 3x3 rotation matrix that transform the lab frame to the object frame
    
end

abstract type Surface end

struct ImplicitSurface <: Surface 
    # define the surface structure
    frame:: Frame # the frame of the surface
    shape:: Function # bound functions specifying the 3d shape, f(x, y, z)=0
    bounds:: Vector{Tuple{Float64, Float64}}
    border:: Function # bound functions specifying the boundary line
end

function shape_plane(center::Vec3)
    return (p) -> p[3] - center[3]
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


function shape_coneflat(center::Vec3, height::Float64, radius::Float64, angle::Float64)
    function result(p)
        x, y, z = p .- center  # Shift coordinates relative to the center
        radial_dist = sqrt(x^2 + y^2)
        r_top = radius
        r_bottom = radius + height * tan(angle)

        if radial_dist <= radius
            return z-height  # Outside the top
        elseif radial_dist <= r_bottom
            return  z-height+(radial_dist-radius)/tan(angle) # Conical region
        else
            return z  # Flat bottom outside cone
        end
    end
    return (p) -> result(p)
end

function border_circle(radius::Float64)
    return (p) -> p⋅p - radius^2  # Circle equation: ||p - center|| = radius
end

function border_bounds()
    return (p) -> -1  #  bound line for surface is defined by bounds
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

function draw_snormals(surface:: Surface, num:: Int,  ray_length::Float64=0.5, arrow_scale::Float64=0.3, color="random")
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
    return draw_rays(rays_snorm, ray_length, arrow_scale, color)
end

# Function to generate and plot rays with arrows
function draw_rays(rays::Vector{Ray}, ray_length::Float64=0.5, arrow_scale::Float64=0.3, color::String="random")
    traces = GenericTrace[]

    color = color == "random" ? "rgb($(rand(0:255)), $(rand(0:255)), $(rand(0:255)))" : color
    # Line traces for the rays
    for ray in rays
        start = ray.origin
        stop = start + normalize(ray.direction) * ray_length
        trace_rays = scatter3d(
            x=[start[1], stop[1]], y=[start[2], stop[2]], z=[start[3], stop[3]],
            mode="lines",
            line=attr(color=color,width=3)
        )
        push!(traces, trace_rays)
    end

    # Arrow traces using cones
    origins = hcat([ray.origin + normalize(ray.direction) * ray_length for ray in rays]...)  # Arrow positions
    directions = hcat([normalize(ray.direction) for ray in rays]...)  # Arrow directions

    trace_arrows = cone(
        x=origins[1, :], y=origins[2, :], z=origins[3, :],  # Positions
        u=directions[1, :], v=directions[2, :], w=directions[3, :],  # Directions
        sizemode="raw",
        sizeref=arrow_scale,  # Controls arrow size
        anchor="tail",
        showscale=false,
        colorscale=[(0, color), (1, color)],
        opacity=0.8
    )
    push!(traces, trace_arrows)

    return traces
end
# ==============================================================================================================


#  Define the structure for material and objects ==========================================================================
abstract type Material end

struct Dielectric <: Material
    n:: Float64
end

struct Metal <: Material
    fuzz:: Float64
    albedo:: Vec3
end

abstract type Object end # object that can be hit by rays
struct BoundObject <: Object
    frame:: Frame # the frame of the object
    surfaces:: Vector{Surface}
    material:: Material
    bounds:: Vector{Tuple{Float64, Float64}}
end

struct Sheet <: Object
    frame:: Frame # the frame of the object
    surfaces:: Vector{Surface}
    material:: Material
end

struct Blocker <: Object
    frame:: Frame # the frame of the object
    surfaces:: Vector{Surface}
end

function draw_object(
    obj::Object, 
    num_sx::Int, 
    num_sy::Int, 
    num_sn::Int,
    ;
    drawnormals::Bool=true, 
    arrow_scale::Float64=1.0,
    color_sf::String="random", 
    color_sn::String="random",
    display::Bool=true)
    objsize = sum([upper - lower for (lower, upper) in obj.bounds])/3.0
    arrsize = objsize/20.0*arrow_scale
    color_sf = color_sf == "random" ? "rgb($(rand(0:255)), $(rand(0:255)), $(rand(0:255)))" : color_sf
    color_sn = color_sn == "random" ? "rgb($(rand(0:255)), $(rand(0:255)), $(rand(0:255)))" : color_sn
    plottraces = GenericTrace[]
    for ss in obj.surfaces
        # rotate and shift the surface by the object's coordinates first
        frame_sf_lab = Frame(obj.frame.origin+ss.frame.origin, obj.frame.orientation*ss.frame.orientation)
        ss_tolab = ImplicitSurface(
            frame_sf_lab, 
            ss.shape, 
            ss.bounds, 
            ss.border
        ) # only work ImplicitSurface for now TODO: generalize it

        append!(plottraces, vcat(
            draw_surface(ss_tolab, num_sx, num_sy, color_sf), 
            draw_snormals(ss_tolab, num_sn, 0.5*arrsize, 0.3*arrsize, color_sn)
            )) 
    end

    if display
        layout = PlotlyJS.Layout(
        title="Box Surface Plot",
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
        fig = PlotlyJS.plot(plottraces , layout);
        PlotlyJS.display(fig)   
        return NaN
    else
        return plottraces
    end
end
# ================================================================================================================

# how the rays interact with the objects and surfaces=============================================================================



# ================================================================================================================

# =============================================================================================================
# ==============================================================================================================


if abspath(@__FILE__) == abspath(PROGRAM_FILE)
    focal_length = 2.0  # Focal length of the paraboloid
    hem_radius = 4.0  # Radius of the hemispherical boundary
    bound_radius = 4.5
    x_range = (-5.0, 5.0)  # X range
    y_range = (-5.0, 5.0)  # Y range
    nx = 50 # number of points along x
    ny = 50 # number of points along y
    num_sn = 100 # number of surface normals
    identity_matrix = rotation_matrix([0.0, 0.0, 1.0], 0.0) 
    rotx_matrix = rotation_matrix([1.0, 0.0, 0.0], π/4.0) 
    roty45_matrix = rotation_matrix([0.0, 1.0, 0.0], -π/4.0) 
    flip_matrix = rotation_matrix([1.0, 0.0, 0.0], π*1.0) 
    
    labframe = Frame([0.0,10.0,0.0], identity_matrix)
    objframe = Frame([10.0,0.0,0.0], rotx_matrix)
    objframe2 = Frame([0.0,-5.0,0.0], roty45_matrix)
    objframe3 = Frame([0.0,-5.0,0.0], flip_matrix)
    bounds = [x_range, y_range]
    
    parasurface = ImplicitSurface(
        labframe, 
        shape_paraboloid([0.0, 0.0, 0.0], focal_length), 
        bounds,
        border_circle(bound_radius))
    
        
    parasurface = ImplicitSurface(
        labframe, 
        shape_paraboloid([0.0, 0.0, 0.0], focal_length), 
        bounds,
        border_circle(bound_radius))
        
    hemisurface = ImplicitSurface(
        objframe, 
        # shape_paraboloid([0.0, 0.0, 0.0], focal_length), 
        shape_hemisphere([0.0, 0.0, 0.0], hem_radius), 
        bounds,
        border_circle(bound_radius))
    
    coneflatsurface = ImplicitSurface(
        objframe2,
        shape_coneflat([0.0, 0.0, -5.0], 4.0, 2.0, π/180*30),
        bounds,
        border_bounds()
        # border_circle([0.0, 0.0], bound_radius)
        )
    
    flatsurfarce1 = ImplicitSurface(
        objframe2, 
        # objframe,
        shape_plane([0.0, 0.0, 5.0]),
        bounds,
        border_bounds()
        # border_circle([0.0, 0.0], bound_radius)
        )
    
    mhatsurface = ImplicitSurface(
        objframe3, 
        shape_mexican_hat([0.0, 0.0, 10.0], 4.0, 1.0, 1.0),
        bounds,
        border_bounds())
    
    
    plottraces = vcat(
        draw_surface(parasurface, nx, ny),
        draw_snormals(parasurface, num_sn),
        draw_surface(hemisurface, nx, ny),
        draw_snormals(hemisurface, num_sn),
        draw_surface(coneflatsurface, nx, ny),
        draw_snormals(coneflatsurface, num_sn),
        draw_surface(flatsurfarce1, nx, ny),
        draw_snormals(flatsurfarce1, num_sn),
        draw_surface(mhatsurface, nx, ny),
        draw_snormals(mhatsurface, num_sn)
    )
    
    
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
    fig = plot(plottraces , layout);
    display(fig)
end
module GeometryModule # module begin----------------------
include("./MathModule.jl")  # Include only if necessary
using Main.MathModule

# Define the structure for surfaces  ===============================================================================

abstract type Surface end
abstract type Volume end

struct ImplicitSurface <: Surface 
    # define the surface structure
    frame:: Frame # the frame of the surface
    shape:: Function # bound functions specifying the 3d shape, f(x, y, z)=0
    bounds:: Vector{Tuple{Float64, Float64}}
    border:: Function # bound functions specifying the boundary line
end

struct BoundVolume <: Volume 
    # define the surface structure
    frame:: Frame # the frame of the surface
    surfaces:: Vector{ImplicitSurface} 
    bounds:: Vector{Tuple{Float64, Float64}}
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
        incap = radius^2 - (p[1] - center[1])^2 - (p[2] - center[2])^2
        return (p[3] - center[3]) - sqrt(max(0, incap))
    end
    return (p) -> result(p)
end


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

function border_outcircle(radius::Float64)
    return (p) -> -p⋅p + radius^2  # Circle equation: ||p - center|| = radius
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






end # module end----------------------

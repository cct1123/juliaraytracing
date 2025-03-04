module GeometryModule # module begin----------------------
using LinearAlgebra
using ForwardDiff
if abspath(PROGRAM_FILE) == @__FILE__
    println("ONLY for DEVELPMENT: Script $(@__FILE__) is running directly!")
    include("./MathModule.jl") 
    using .MathModule
    
    # Do something here
else
    using Main.MathModule
end

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

# some surface shapes ---------------------
function surface_normal(p, surface::Surface)
    # Calculate the gradient of the surface shape function
    grad_S = ForwardDiff.gradient(surface.shape, p)

    # Normalize the gradient to get the normal
    return grad_S ./ norm(grad_S)
end

function shape_plane(center::Vec3)
    return (p) -> p[3] - center[3]
end

# example shape function: Define a spherical surface as the shape function
function shape_sphere(center::Vec3, radius::Float64)
    return (p) -> dot((p - center), (p - center)) - radius^2  # Sphere equation: ||p - center|| = radius
end

function shape_hemisphere(center::Vec3, radius::Float64)
    function result(p)
        incap = radius^2 - (p[1] - center[1])^2 - (p[2] - center[2])^2
        return (p[3] - center[3]) - sqrt(max(0, incap))
    end
    return (p) -> result(p)
end

function shape_cylinder(center::Vec3, radius::Float64, height::Float64)
    return (p) -> begin
        # Cylinder equation: (x^2 + y^2) = r^2, with z within [0, height]
        x, y, z = p - center
        if z >= 0 && z <= height
            return x^2 + y^2 - radius^2 
        elseif z > height
            return z - height
        elseif z < 0
            return -z
        end
    end
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

# function border_polygon(base_vertices::Vector{Vector{Float64}})
#     n = length(base_vertices)
#     return (p) -> begin
#         # Check if the point is inside or on the polygon using the winding number approach
#         wn = 0  # Winding number
#         for i in 1:n
#             v1 = base_vertices[i]
#             v2 = base_vertices[mod(i, n) + 1]
            
#             # Check if the point is on the edge
#             if is_point_on_edge(v1, v2, p)
#                 return true  # Point is on the border
#             end
            
#             # Check the winding number
#             if (v1[2] <= p[2] && v2[2] > p[2]) || (v1[2] > p[2] && v2[2] <= p[2])
#                 if cross_product(v2 - v1, p - v1) > 0
#                     wn += 1
#                 else
#                     wn -= 1
#                 end
#             end
#         end
        
#         return wn != 0  # Point is inside the polygon if the winding number is non-zero
#     end
# end

# # Helper function to check if the point is on the edge
# function is_point_on_edge(v1, v2, p)
#     cross_prod = cross_product(v2 - v1, p - v1)
#     return isapprox(cross_prod, 0.0) && dot(v2 - v1, p - v1) >= 0 && dot(v1 - v2, p - v2) >= 0
# end

# Helper function to compute the 2D cross product
function cross_product(v1, v2)
    return v1[1] * v2[2] - v1[2] * v2[1]
end

function border_bounds()
    return (p) -> -1  #  bound line for surface is defined by bounds
end

# some common Volume----------------------------------
const cen000 = [0.0, 0.0, 0.0]
const borderbds = border_bounds()

function make_plane(
    lx::Float64, ly::Float64, 
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
    )  
    return ImplicitSurface(
        Frame(center, orientation),
        shape_plane(cen000),
        [(-lx/2.0, lx/2.0), (-ly/2.0, ly/2.0)],
        borderbds   
    )
end

function make_circleplane(
    radius::Float64, 
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
    )  
    # make a circle
    return ImplicitSurface(
        Frame(center, orientation),
        shape_plane(cen000),
        [(-radius, radius), (-radius, radius)],
        border_circle(radius)   
    )
end

function make_box(
    lx::Float64, ly::Float64, lz::Float64, 
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
    )  
    # the origin is defined at the center of the box volume 
    sf_top = ImplicitSurface(
        Frame([0.0, 0.0, lz/2.0], rm_eye),
        shape_plane(cen000),
        [(-lx/2.0, lx/2.0), (-ly/2.0, ly/2.0)],
        borderbds
    )

    sf_bottom = ImplicitSurface(
        Frame([0.0, 0.0, -lz/2.0], rm_xflip),
        shape_plane(cen000),
        [(-lx/2.0, lx/2.0), (-ly/2.0, ly/2.0)],
        borderbds,
        )

    sf_side_r = ImplicitSurface(
        Frame([lx/2.0, 0.0, 0.0], rm_yp90),
        shape_plane(cen000),
        [(-lz/2.0, lz/2.0), (-ly/2.0, ly/2.0)],
        borderbds
        )

    sf_side_l = ImplicitSurface(
        Frame([-lx/2.0, 0.0, 0.0], rm_yn90),
        shape_plane(cen000),
        [(-lz/2.0, lz/2.0), (-ly/2.0, ly/2.0)],
        borderbds
        )

    sf_side_f = ImplicitSurface(
        Frame([0.0, ly/2.0, 0.0], rm_xn90),
        shape_plane(cen000),
        [ (-lx/2.0, lx/2.0), (-lz/2.0, lz/2.0)],
        borderbds
        )

    sf_side_b = ImplicitSurface(
        Frame([0.0, -ly/2.0, 0.0], rm_xp90),
        shape_plane(cen000),
        [ (-lx/2.0, lx/2.0), (-lz/2.0, lz/2.0)],
        borderbds
        )

    box_surfaces = [sf_top, sf_bottom, sf_side_r, sf_side_l, sf_side_f, sf_side_b]

    bounds_xyz = [(-lx*1.05/2.0, lx*1.05/2.0), (-ly*1.05/2.0, ly*1.05/2.0), (-lz*1.05/2.0, lz*1.05/2.0)]
    return BoundVolume(
                Frame(center, orientation),
                box_surfaces,
                bounds_xyz
            )
end

function make_cone(
    r_top::Float64, r_base::Float64, height::Float64, 
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
    )
    # the origin is defined at the center at the base plane
    θ_cone = atan((r_base-r_top)/height)
    bounds_xyz = [(-r_base*1.05, r_base*1.05), (-r_base*1.05, r_base*1.05), (-height*0.05, height*1.05)]
    
    sf_coneflat = ImplicitSurface(
        Frame(cen000, rm_eye),
        shape_coneflat(cen000, height, r_top, θ_cone),
        bounds_xyz[1:2],
        border_circle(r_base)
    )

    sf_base = ImplicitSurface(
        Frame(cen000, rm_xflip),
        shape_plane(cen000),
        bounds_xyz[1:2],
        border_circle(r_base)
    )

    surfaces = [sf_coneflat, sf_base]
    return BoundVolume(
                Frame(center, orientation),
                surfaces,
                bounds_xyz
            )
end

function make_sphere(
    radius::Float64, 
    center::Vec3=cen000
    )
    r_bound = radius*1.05
    shapehemisphere = shape_hemisphere(cen000, radius)
    bordercircle = border_circle(radius)
    bounds_xyz = [
        (-r_bound, r_bound), 
        (-r_bound, r_bound), 
        (-r_bound, r_bound)]
    
    sf_hemi_t = ImplicitSurface(
        Frame(cen000, rm_eye), 
        shapehemisphere, 
        bounds_xyz[1:2],
        bordercircle
    )
    sf_hemi_b = ImplicitSurface(
        Frame(cen000, rm_xflip), 
        shapehemisphere, 
        bounds_xyz[1:2],
        bordercircle
    )
    
    surfaces = [sf_hemi_t, sf_hemi_b]
    return BoundVolume(
        Frame(center, rm_eye), 
        surfaces, 
        bounds_xyz
    )
end

function make_celestial(
    radius::Float64, 
    center::Vec3=cen000
    )
    r_bound = radius*1.05
    bordercircle = border_circle(radius)
    bounds_xyz = [
        (-r_bound, r_bound), 
        (-r_bound, r_bound), 
        (-r_bound, r_bound)]
    
    sf_sphere = ImplicitSurface(
        Frame(cen000, rm_eye), 
        shape_sphere(cen000, radius), 
        bounds_xyz[1:2],
        bordercircle
    )

    return BoundVolume(
        Frame(center, rm_eye), 
        [sf_sphere], 
        bounds_xyz
    )
end

# below are untested functions from ChatGPT-----
# function make_prism(
#     base_vertices::Vector{Vec3}, height::Float64,
#     center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
#     )
#     # Define bounding box
#     min_x = minimum(v[1] for v in base_vertices)
#     max_x = maximum(v[1] for v in base_vertices)
#     min_y = minimum(v[2] for v in base_vertices)
#     max_y = maximum(v[2] for v in base_vertices)
#     bounds_xyz = [(min_x, max_x), (min_y, max_y), (0.0, height)]

#     # Base and top surfaces
#     sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_plane(cen000), bounds_xyz[1:2], border_polygon(base_vertices))
#     sf_top = ImplicitSurface(Frame([0.0, 0.0, height], rm_eye), shape_plane(cen000), bounds_xyz[1:2], border_polygon(base_vertices))
    
#     # Side faces
#     side_surfaces = [
#         ImplicitSurface(Frame([(v1 + v2) / 2..., 0.0], rm_yp90), shape_plane(cen000), [(0.0, height), (-norm(v2 - v1) / 2, norm(v2 - v1) / 2)], borderbds)
#         for (v1, v2) in zip(base_vertices, circshift(base_vertices, -1))
#     ]
    
#     return BoundVolume(Frame(center, orientation), [sf_base, sf_top, side_surfaces], bounds_xyz)
# end

function make_cylinder(
    radius::Float64, height::Float64,
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
    )
    bounds_xyz = [(-radius, radius), (-radius, radius), (0.0, height)]

    sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_plane(cen000), bounds_xyz[1:2], border_circle(radius))
    sf_top = ImplicitSurface(Frame([0.0, 0.0, height], rm_eye), shape_plane(cen000), bounds_xyz[1:2], border_circle(radius))
    sf_side = ImplicitSurface(Frame(cen000, rm_eye), shape_cylinder(cen000, radius, height), bounds_xyz[1:2], borderbds)

    return BoundVolume(Frame(center, orientation), [sf_base, sf_top, sf_side], bounds_xyz)
end

# function make_ring(
#     r_inner::Float64, r_outer::Float64, height::Float64,
#     center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
#     )
#     bounds_xyz = [(-r_outer, r_outer), (-r_outer, r_outer), (0.0, height)]

#     sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_annulus(r_inner, r_outer), bounds_xyz[1:2], border_annulus(r_inner, r_outer))
#     sf_top = ImplicitSurface(Frame([0.0, 0.0, height], rm_eye), shape_annulus(r_inner, r_outer), bounds_xyz[1:2], border_annulus(r_inner, r_outer))
#     sf_outer = ImplicitSurface(Frame(cen000, rm_eye), shape_cylinder(cen000, r_outer, height), bounds_xyz[1:2], borderbds)
#     sf_inner = ImplicitSurface(Frame(cen000, rm_eye), shape_cylinder(cen000, r_inner, height), bounds_xyz[1:2], borderbds)

#     return BoundVolume(Frame(center, orientation), [sf_base, sf_top, sf_outer, sf_inner], bounds_xyz)
# end

# function make_hemisphere(
#     radius::Float64,
#     center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
#     )
#     bounds_xyz = [(-radius, radius), (-radius, radius), (0.0, radius)]

#     sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_plane(cen000), bounds_xyz[1:2], border_circle(radius))
#     sf_cap = ImplicitSurface(Frame(cen000, rm_eye), shape_sphere(cen000, radius), bounds_xyz[1:2], borderbds)
    
#     return BoundVolume(Frame(center, orientation), [sf_base, sf_cap], bounds_xyz)
# end

# function make_spherical_lens(
#     r1::Float64, r2::Float64, thickness::Float64,
#     center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
#     )
#     bounds_xyz = [(-r1, r1), (-r1, r1), (-thickness / 2, thickness / 2)]

#     sf_front = ImplicitSurface(Frame([0.0, 0.0, thickness / 2], rm_eye), shape_sphere(cen000, r1), bounds_xyz[1:2], border_circle(r1))
#     sf_back = ImplicitSurface(Frame([0.0, 0.0, -thickness / 2], rm_xflip), shape_sphere(cen000, r2), bounds_xyz[1:2], border_circle(r2))
    
#     return BoundVolume(Frame(center, orientation), [sf_front, sf_back], bounds_xyz)
# end
# above are untested functions from ChatGPT-----

export ImplicitSurface, BoundVolume, Surface, Volume, 
       shape_plane, shape_sphere, shape_hemisphere, shape_paraboloid, shape_mexican_hat, shape_coneflat, surface_normal,
       border_circle, border_outcircle, border_bounds,
       make_plane, make_box, make_cone, make_sphere,
       make_cylinder#, make_ring, make_hemisphere, make_spherical_lens

end # module end----------------------

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

# some surface shapes ---------------------
function surface_normal(p, surface::Surface)
    return surface_normal(p, surface.shape)
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

# some common Volume----------------------------------
rm_eye = rotation_matrix([0.0, 0.0, 1.0], 0.0) 
rm_xflip = rotation_matrix([1.0, 0.0, 0.0], 1.0*π)
rm_xp90 = rotation_matrix([1.0, 0.0, 0.0], 0.5*π)
rm_xn90 = rotation_matrix([1.0, 0.0, 0.0], -0.5*π)

rm_yflip = rotation_matrix([0.0, 1.0, 0.0], 1.0*π)
rm_yp90 = rotation_matrix([0.0, 1.0, 0.0], 0.5*π)
rm_yn90 = rotation_matrix([0.0, 1.0, 0.0], -0.5*π)

cen000 = [0.0, 0.0, 0.0]
borderbds = border_bounds()

function make_box(
    lx::Float64, ly::Float64, lz::Float64, 
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
    )  
    # the origin is defined at the center of the box volume 
    sf_top = ImplicitSurface(
        Frame([0.0, 0.0, lz/2.0], rm_eye),
        shape_plane(cen000),
        bounds_xy,
        borderbds
    )

    sf_bottom = ImplicitSurface(
        Frame([0.0, 0.0, -lz/2.0], rm_xflip),
        shape_plane(cen000),
        bounds_xy,
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
    return BoundVlolume(
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
    bounds_xyz = [(-r_base*1.05, r_base*1.05), (-r_base*1.05, r_base*1.05), [(-height*0.05, height*1.05)]]
    
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
    return BoundVlolume(
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
    return BoundObject(
        Frame(center, rm_eye), 
        surfaces, 
        bounds_xyz
    )
end

# below are untested functions from ChatGPT-----
function make_prism(
    base_vertices::Vector{Vec3}, height::Float64,
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
)
    # Define bounding box
    min_x = minimum(v[1] for v in base_vertices)
    max_x = maximum(v[1] for v in base_vertices)
    min_y = minimum(v[2] for v in base_vertices)
    max_y = maximum(v[2] for v in base_vertices)
    bounds_xyz = [(min_x, max_x), (min_y, max_y), (0.0, height)]

    # Base and top surfaces
    sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_polygon(base_vertices), bounds_xyz[1:2], border_polygon(base_vertices))
    sf_top = ImplicitSurface(Frame([0.0, 0.0, height], rm_eye), shape_polygon(base_vertices), bounds_xyz[1:2], border_polygon(base_vertices))
    
    # Side faces
    side_surfaces = [
        ImplicitSurface(Frame([(v1 + v2) / 2..., 0.0], rm_yp90), shape_plane(cen000), [(0.0, height), (-norm(v2 - v1) / 2, norm(v2 - v1) / 2)], borderbds)
        for (v1, v2) in zip(base_vertices, circshift(base_vertices, -1))
    ]
    
    return BoundVolume(Frame(center, orientation), [sf_base, sf_top; side_surfaces], bounds_xyz)
end

function make_cylinder(
    radius::Float64, height::Float64,
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
)
    bounds_xyz = [(-radius, radius), (-radius, radius), (0.0, height)]

    sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_plane(cen000), bounds_xyz[1:2], border_circle(radius))
    sf_top = ImplicitSurface(Frame([0.0, 0.0, height], rm_eye), shape_plane(cen000), bounds_xyz[1:2], border_circle(radius))
    sf_side = ImplicitSurface(Frame(cen000, rm_eye), shape_cylinder(radius, height), bounds_xyz[1:2], borderbds)

    return BoundVolume(Frame(center, orientation), [sf_base, sf_top, sf_side], bounds_xyz)
end

function make_ring(
    r_inner::Float64, r_outer::Float64, height::Float64,
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
)
    bounds_xyz = [(-r_outer, r_outer), (-r_outer, r_outer), (0.0, height)]

    sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_annulus(r_inner, r_outer), bounds_xyz[1:2], border_annulus(r_inner, r_outer))
    sf_top = ImplicitSurface(Frame([0.0, 0.0, height], rm_eye), shape_annulus(r_inner, r_outer), bounds_xyz[1:2], border_annulus(r_inner, r_outer))
    sf_outer = ImplicitSurface(Frame(cen000, rm_eye), shape_cylinder(r_outer, height), bounds_xyz[1:2], borderbds)
    sf_inner = ImplicitSurface(Frame(cen000, rm_eye), shape_cylinder(r_inner, height), bounds_xyz[1:2], borderbds)

    return BoundVolume(Frame(center, orientation), [sf_base, sf_top, sf_outer, sf_inner], bounds_xyz)
end

function make_hemisphere(
    radius::Float64,
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
)
    bounds_xyz = [(-radius, radius), (-radius, radius), (0.0, radius)]

    sf_base = ImplicitSurface(Frame(cen000, rm_xflip), shape_plane(cen000), bounds_xyz[1:2], border_circle(radius))
    sf_cap = ImplicitSurface(Frame(cen000, rm_eye), shape_sphere(cen000, radius), bounds_xyz[1:2], borderbds)
    
    return BoundVolume(Frame(center, orientation), [sf_base, sf_cap], bounds_xyz)
end

function make_spherical_lens(
    r1::Float64, r2::Float64, thickness::Float64,
    center::Vec3=cen000, orientation::Matrix{Float64}=rm_eye
)
    bounds_xyz = [(-r1, r1), (-r1, r1), (-thickness / 2, thickness / 2)]

    sf_front = ImplicitSurface(Frame([0.0, 0.0, thickness / 2], rm_eye), shape_sphere(cen000, r1), bounds_xyz[1:2], border_circle(r1))
    sf_back = ImplicitSurface(Frame([0.0, 0.0, -thickness / 2], rm_xflip), shape_sphere(cen000, r2), bounds_xyz[1:2], border_circle(r2))
    
    return BoundVolume(Frame(center, orientation), [sf_front, sf_back], bounds_xyz)
end
# above are untested functions from ChatGPT-----



end # module end----------------------

# rtweekend.jl
include("rtweekend.jl")
using .RTWeekend
# Helper function to check if a vector is near zero
function near_zero(v::Vec3)
    s = 1e-8
    return abs(v.x) < s && abs(v.y) < s && abs(v.z) < s
end
# "hittable.jl"
# include("hittable.jl")
abstract type Material end
mutable struct HitRecord
    p::Vec3                  # Intersection point
    normal::Vec3             # Surface normal
    t::Float64               # Ray parameter at intersection
    front_face::Bool         # Whether the ray hits the front face
    material::Material       # Material at the intersection point
end

function set_face_normal!(rec::HitRecord, ray::Ray, outward_normal::Vec3)
    rec.front_face = dot(ray.direction, outward_normal) < 0
    rec.normal = rec.front_face ? outward_normal : -outward_normal
end


# Hittable: Abstract type for objects that can be hit by a ray
abstract type Hittable end

function hit(obj::Hittable, ray::Ray, ray_t::Interval)::Union{HitRecord, Nothing}
    error("hit() not implemented for this Hittable type")
end

# sphere hittable object
struct Sphere <: Hittable
    center::Vec3
    radius::Float64
    material::Material
end

function hit(sphere::Sphere, ray::Ray, ray_t::Interval)::Union{HitRecord, Nothing}
    oc = sphere.center - ray.origin
    a = dot(ray.direction, ray.direction)
    h = dot(ray.direction, oc)
    c = dot(oc, oc) - sphere.radius^2
    discriminant = h^2 - a * c

    if discriminant < 0
        return nothing  # No intersection
    end

    sqrtd = sqrt(discriminant)

    # Find the nearest root within the acceptable range
    root = (h - sqrtd) / a
    if !surrounds(ray_t, root)
        root = (h + sqrtd) / a
        if !surrounds(ray_t, root)
            return nothing  # No valid intersection
        end
    end

    # Compute intersection point and normal
    p = at(ray, root)
    outward_normal = (p - sphere.center) / sphere.radius
    rec = HitRecord(p, outward_normal, root, true, sphere.material)
    set_face_normal!(rec, ray, outward_normal)

    return rec
end


# HittableList: A list of hittable objects
struct HittableList <: Hittable
    objects::Vector{Hittable}
end

HittableList() = HittableList(Vector{Hittable}())

function add!(list::HittableList, obj::Hittable)
    push!(list.objects, obj)
end

function clear!(list::HittableList)
    empty!(list.objects)
end

function hit(list::HittableList, ray::Ray, ray_t::Interval)::Union{HitRecord, Nothing}
    closest_so_far = ray_t.max
    hit_anything = nothing

    for obj in list.objects
        temp_rec = hit(obj, ray, Interval(ray_t.min, closest_so_far))
        if temp_rec !== nothing
            closest_so_far = temp_rec.t
            hit_anything = temp_rec
        end
    end

    return hit_anything
end


# some scattering and reflection

# Scatter function for materials
function scatter(material::Material, ray_in::Ray, rec::HitRecord)
    # Default behavior: no scattering
    return false, Vec3(0.0, 0.0, 0.0), Ray(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0))
end

struct Lambertian <: Material
    albedo::Vec3  # Surface color (reflectance)
end

function scatter(material::Lambertian, ray_in::Ray, rec::HitRecord)
    scatter_direction = rec.normal + random_unit_vector()

    # Catch degenerate scatter direction
    if near_zero(scatter_direction)
        scatter_direction = rec.normal
    end

    ray_scattered = Ray(rec.p, scatter_direction)
    attenuation = material.albedo
    return true, attenuation, ray_scattered
end

struct Metal <: Material
    albedo::Vec3  # Surface color (reflectance)
end

function reflect(v::Vec3, n::Vec3)
    return v - 2 * dot(v, n) * n
end

function scatter(material::Metal, ray_in::Ray, rec::HitRecord)
    ray_reflected = reflect(normalize(ray_in.direction), rec.normal)
    ray_scattered = Ray(rec.p, ray_reflected)
    attenuation = material.albedo
    return (dot(ray_scattered.direction, rec.normal) > 0), attenuation, ray_scattered
end



# camera.jl
include("camera.jl")
function ray_color(camera::Camera, ray::Ray, world::HittableList, depth::Int)
    if depth <= 0
        return Vec3(0.0, 0.0, 0.0) # No light contribution at max depth
    end

    rec = hit(world, ray, Interval(0.001, RTWeekend.infinity))
    if rec !== nothing
        scattered, attenuation, ray_scattered = scatter(rec.material, ray, rec)
        if scattered
            return attenuation * ray_color(camera, ray_scattered, world, depth - 1)
        else
            return Vec3(0.0, 0.0, 0.0)
        end
    end

    # Background gradient
    unit_direction = normalize(ray.direction)
    t = 0.5 * (unit_direction.y + 1.0)
    return (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0)
end


#main.jl
function main()
    # World setup
    world = HittableList()

    material_ground = Lambertian(Vec3(0.8, 0.8, 0.0))
    material_center = Lambertian(Vec3(0.1, 0.2, 0.5))
    material_left = Metal(Vec3(0.8, 0.8, 0.8))
    material_right = Metal(Vec3(0.8, 0.6, 0.2))

    add!(world, Sphere(Vec3(0.0, -100.5, -1.0), 100.0, material_ground))
    add!(world, Sphere(Vec3(0.0, 0.0, -1.2), 0.5, material_center))
    add!(world, Sphere(Vec3(-1.0, 0.0, -1.0), 0.5, material_left))
    add!(world, Sphere(Vec3(1.0, 0.0, -1.0), 0.5, material_right))

    # Camera setup
    cam = Camera(aspect_ratio=16.0 / 9.0, image_width=400, samples_per_pixel=100)

    # Render the scene
    render(cam, world, "images/ch10_image.ppm")
    ppm_to_jpg("images/ch10_image.ppm", "images/ch10_image.jpg")
end

main()
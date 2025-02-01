# main.jl
include("rtweekend.jl")
using .RTWeekend


# HitRecord: Stores information about a ray-object intersection
mutable struct HitRecord
    p::Vec3
    normal::Vec3
    t::Float64
    front_face::Bool
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

# Sphere: A hittable sphere
struct Sphere <: Hittable
    center::Vec3
    radius::Float64
end

function hit(sphere::Sphere, ray::Ray, ray_t::Interval)::Union{HitRecord, Nothing}
    oc = sphere.center - ray.origin
    a = dot(ray.direction, ray.direction)
    h = dot(ray.direction, oc)
    c = dot(oc, oc) - sphere.radius^2
    discriminant = h^2 - a*c

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
    rec = HitRecord(p, outward_normal, root, true)
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

# Ray Color Function
function ray_color(ray::Ray, world::HittableList)
    rec = hit(world, ray, Interval(0.001, RTWeekend.infinity))
    if rec !== nothing
        return 0.5 * Vec3(rec.normal.x + 1.0, rec.normal.y + 1.0, rec.normal.z + 1.0)
    end

    # Background gradient
    unit_direction = normalize(ray.direction)
    t = 0.5 * (unit_direction.y + 1.0)
    return (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0)
end

# Main Function
function main()
    # Image setup
    aspect_ratio = 16.0 / 9.0
    image_width = 400
    image_height = max(1, trunc(Int, image_width / aspect_ratio))
    
    # Camera setup
    viewport_height = 2.0
    viewport_width = viewport_height * (image_width / image_height)
    focal_length = 1.0
    
    origin = Vec3(0.0, 0.0, 0.0)
    horizontal = Vec3(viewport_width, 0.0, 0.0)
    vertical = Vec3(0.0, -viewport_height, 0.0)
    lower_left_corner = origin - horizontal/2.0 - vertical/2.0 - Vec3(0.0, 0.0, focal_length)
    
    # World setup
    world = HittableList()
    add!(world, Sphere(Vec3(0.0, 0.0, -1.0), 0.5))
    add!(world, Sphere(Vec3(0.0, -100.5, -1.0), 100.0))
    
    # Render to PPM
    ppmfile_name = "images/ch6.8_image.ppm"
    open(ppmfile_name, "w") do file
        println(file, "P3")
        println(file, image_width, " ", image_height)
        println(file, "255")
        
        for j in 0:image_height-1
            for i in 0:image_width-1
                u = i / (image_width - 1)
                v = j / (image_height - 1)
                direction = lower_left_corner + horizontal*u + vertical*v - origin
                ray = Ray(origin, direction)
                col = ray_color(ray, world)
                
                ir = trunc(Int, 255.999 * col.x)
                ig = trunc(Int, 255.999 * col.y)
                ib = trunc(Int, 255.999 * col.z)
                
                println(file, ir, " ", ig, " ", ib)
            end
        end
    end
    println("Rendering complete. Image saved as $ppmfile_name")
    ppm_to_jpg(ppmfile_name, "images/ch6.8_image.jpg")
end

# Run the ray tracer
main()
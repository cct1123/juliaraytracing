include("ch4.1_ray.jl")
include("ch3_vec3.jl")
include("ppm_to_jpg.jl")


mutable struct HitRecord
    p::Vec3
    normal::Vec3
    t::Float64
    front_face::Bool
end

abstract type Hittable end

struct Sphere <: Hittable
    center::Vec3
    radius::Float64
end

function set_face_normal(rec::HitRecord, ray::Ray, outward_normal::Vec3)
    # Determine if the ray is hitting the front or back face
    rec.front_face = dot(ray.direction, outward_normal) < 0
    # Set the normal to face against the ray
    rec.normal = rec.front_face ? outward_normal : -outward_normal
    return rec
end

function hit(sphere::Sphere, ray::Ray, t_min::Float64, t_max::Float64)::Union{HitRecord, Nothing}
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
    if root <= t_min || root >= t_max
        root = (h + sqrtd) / a
        if root <= t_min || root >= t_max
            return nothing  # No valid intersection
        end
    end

    # Compute intersection point and outward normal
    p = at(ray, root)
    outward_normal = (p - sphere.center) / sphere.radius

    # Create a HitRecord and set the face normal
    rec = HitRecord(p, outward_normal, root, true)  # Default front_face to true
    rec = set_face_normal(rec, ray, outward_normal)

    return rec
end

function ray_color(ray::Ray, world::Vector{Hittable})
    closest_so_far = Inf
    hit_record = nothing

    # Check for intersections with all objects in the world
    for obj in world
        temp_rec = hit(obj, ray, 0.001, closest_so_far)
        if temp_rec !== nothing
            closest_so_far = temp_rec.t
            hit_record = temp_rec
        end
    end

    # If we hit something, color it using the normal
    if hit_record !== nothing
        return 0.5 * Vec3(hit_record.normal.x + 1.0, hit_record.normal.y + 1.0, hit_record.normal.z + 1.0)
    end

    # Otherwise, return the background gradient
    unit_direction = normalize(ray.direction)
    t = 0.5 * (unit_direction.y + 1.0)
    return (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0)
end


function main()
    # Image setup
    aspect_ratio = 16.0 / 9.0
    image_width = 400
    image_height = trunc(Int, image_width / aspect_ratio)
    
    # Camera setup
    viewport_height = 2.0
    viewport_width = aspect_ratio * viewport_height
    focal_length = 1.0
    
    origin = Vec3(0.0, 0.0, 0.0)
    horizontal = Vec3(viewport_width, 0.0, 0.0)
    vertical = Vec3(0.0, viewport_height, 0.0)
    lower_left_corner = origin - horizontal/2.0 - vertical/2.0 - Vec3(0.0, 0.0, focal_length)
    
    # World setup
    world = Hittable[
        Sphere(Vec3(0.0, 0.0, -1.0), 0.5),
        Sphere(Vec3(0.0, -100.5, -1.0), 100.0)
    ]
    
    # Render to PPM
    ppmfile_name = "images/ch6.4_image.ppm"
    open(ppmfile_name, "w") do file
        println(file, "P3")
        println(file, image_width, " ", image_height)
        println(file, "255")
        
        for j in image_height-1:-1:0
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
    ppm_to_jpg(ppmfile_name, "images/ch6.4_image.jpg")
end

main()
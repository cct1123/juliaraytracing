include("ppm_to_jpg.jl")
include("ch4.1_ray.jl")
include("ch3_vec3.jl")
include("ch6.3_hittable_abstraction.jl")
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

function hit(list::HittableList, ray::Ray, t_min::Float64, t_max::Float64)::Union{HitRecord, Nothing}
    closest_so_far = t_max
    hit_anything = nothing

    for obj in list.objects
        temp_rec = hit(obj, ray, t_min, closest_so_far)
        if temp_rec !== nothing
            closest_so_far = temp_rec.t
            hit_anything = temp_rec
        end
    end

    return hit_anything
end

function ray_color(ray::Ray, world::HittableList)
    rec = hit(world, ray, 0.001, Inf)
    if rec !== nothing
        return 0.5 * Vec3(rec.normal.x + 1.0, rec.normal.y + 1.0, rec.normal.z + 1.0)
    end

    # Background gradient
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
    world = HittableList()
    add!(world, Sphere(Vec3(0.0, 0.0, -1.0), 0.5))
    add!(world, Sphere(Vec3(0.0, -100.5, -1.0), 100.0))
    
    # Render to PPM
    ppmfile_name = "images/ch6.5_image.ppm"
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
    ppm_to_jpg(ppmfile_name, "images/ch6.5_image.jpg")
end

main()
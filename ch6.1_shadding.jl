include("ch4.1_ray.jl")
include("ch3_vec3.jl")

function hit_sphere(center::Vec3, radius::Float64, ray::Ray)
    oc = center - ray.origin
    a = dot(ray.direction, ray.direction)
    b = -2.0 * dot(ray.direction, oc)
    c = dot(oc, oc) - radius^2
    discriminant = b^2 - 4*a*c

    if discriminant < 0
        return -1.0  # No intersection
    else
        return (-b - sqrt(discriminant)) / (2.0*a)  # Return the smallest t
    end
end

function ray_color(ray::Ray)
    t = hit_sphere(Vec3(0.0, 0.0, -1.0), 0.5, ray)
    if t > 0.0
        # Compute the intersection point
        intersection_point = at(ray, t)
        # Compute the normal vector (unit length)
        normal = normalize(intersection_point - Vec3(0.0, 0.0, -1.0))
        # Map normal components from [-1, 1] to [0, 1] for color
        return 0.5 * Vec3(normal.x + 1.0, normal.y + 1.0, normal.z + 1.0)
    end

    # Background gradient (unchanged)
    unit_direction = normalize(ray.direction)
    t = 0.5 * (unit_direction.y + 1.0)
    return (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0)
end

include("ch4.1_ray.jl")
include("ch3_vec3.jl")

# Ensure images directory exists
if !isdir("images")
    mkdir("images")
end

function hit_sphere(center::Vec3, radius::Float64, ray::Ray)
    oc = center - ray.origin
    a = dot(ray.direction, ray.direction)
    b = -2.0 * dot(ray.direction, oc)
    c = dot(oc, oc) - radius^2
    discriminant = b^2 - 4*a*c

    if discriminant < 0
        return -1.0  # No intersection
    else
        return (-b - sqrt(discriminant)) / (2.0*a)  # Return the smallest t
    end
end

function ray_color(ray::Ray)
    t = hit_sphere(Vec3(0.0, 0.0, -1.0), 0.5, ray)
    if t > 0.0
        # Compute the intersection point
        intersection_point = at(ray, t)
        # Compute the normal vector (unit length)
        normal = normalize(intersection_point - Vec3(0.0, 0.0, -1.0))
        # Map normal components from [-1, 1] to [0, 1] for color
        return 0.5 * Vec3(normal.x + 1.0, normal.y + 1.0, normal.z + 1.0)
    end

    # Background gradient (unchanged)
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
    
    # Render to PPM
    open("images/ch6_image.ppm", "w") do file
        println(file, "P3")
        println(file, image_width, " ", image_height)
        println(file, "255")
        
        for j in image_height-1:-1:0
            for i in 0:image_width-1
                u = i / (image_width - 1)
                v = j / (image_height - 1)
                direction = lower_left_corner + horizontal*u + vertical*v - origin
                ray = Ray(origin, direction)
                col = ray_color(ray)
                
                ir = trunc(Int, 255.999 * col.x)
                ig = trunc(Int, 255.999 * col.y)
                ib = trunc(Int, 255.999 * col.z)
                
                println(file, ir, " ", ig, " ", ib)
            end
        end
    end
end

main()
include("rtweekend.jl")
using .RTWeekend
include("hittable.jl")
# rtweekend.jl
function random_double()
    rand()
end

function random_double(min::Float64, max::Float64)
    min + (max - min) * random_double()
end

function random_vec3(min::Float64, max::Float64)
    return Vec3(random_double(min, max), random_double(min, max), random_double(min, max))
end

function length_squared(v::Vec3)
    return norm(v)^2
end

function unit_vector(v::Vec3)
    return normalize(v)
end

function random_unit_vector()
    while true
        p = random_vec3(-1.0, 1.0)
        len_sq = length_squared(p)
        if 1e-160 < len_sq <= 1.0
            return unit_vector(p)
        end
    end
end

function random_on_hemisphere(normal::Vec3)
    on_unit_sphere = random_unit_vector()
    if dot(on_unit_sphere, normal) > 0.0
        return on_unit_sphere
    else
        return -on_unit_sphere
    end
end

function my_clamp(x, min_val, max_val)
    if x < min_val
        return min_val
    elseif x > max_val
        return max_val
    else
        return x
    end
end

# camera.jl
struct Camera
    aspect_ratio::Float64      # Ratio of image width over height
    image_width::Int           # Rendered image width in pixel count
    image_height::Int          # Rendered image height (calculated)
    samples_per_pixel::Int     # Number of random samples per pixel
    center::Vec3               # Camera center
    pixel00_loc::Vec3          # Location of pixel (0, 0)
    pixel_delta_u::Vec3        # Offset to pixel to the right
    pixel_delta_v::Vec3        # Offset to pixel below
    pixel_samples_scale::Float64 # Color scale factor for a sum of pixel samples
end

function Camera(; aspect_ratio=16.0 / 9.0, image_width=400, samples_per_pixel=100)
    # Calculate image height and ensure it's at least 1
    image_height = max(1, trunc(Int, image_width / aspect_ratio))

    # Camera center
    center = Vec3(0.0, 0.0, 0.0)

    # Determine viewport dimensions
    focal_length = 1.0
    viewport_height = 2.0
    viewport_width = viewport_height * (image_width / image_height)

    # Calculate vectors across the horizontal and down the vertical viewport edges
    viewport_u = Vec3(viewport_width, 0.0, 0.0)
    viewport_v = Vec3(0.0, -viewport_height, 0.0)

    # Calculate the horizontal and vertical delta vectors
    pixel_delta_u = viewport_u / image_width
    pixel_delta_v = viewport_v / image_height

    # Calculate the location of the upper-left pixel
    viewport_upper_left = center - Vec3(0.0, 0.0, focal_length) - viewport_u/2 - viewport_v/2
    pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v)

    # Color scale factor for a sum of pixel samples
    pixel_samples_scale = 1.0 / samples_per_pixel

    # Create and return the Camera object
    Camera(aspect_ratio, image_width, image_height, samples_per_pixel, center, pixel00_loc, pixel_delta_u, pixel_delta_v, pixel_samples_scale)
end

function sample_square()
    # Returns a random point in the square surrounding a pixel
    Vec3(random_double() - 0.5, random_double() - 0.5, 0.0)
end

function get_ray(camera::Camera, i::Int, j::Int)
    # Get a ray for the pixel at location (i, j) with random sampling
    pixel_sample = camera.pixel00_loc + (i + random_double() - 0.5) * camera.pixel_delta_u + (j + random_double() - 0.5) * camera.pixel_delta_v
    ray_direction = pixel_sample - camera.center
    Ray(camera.center, ray_direction)
end

function ray_color(camera::Camera, ray::Ray, world::HittableList, depth::Int)
    if depth <= 0
        return Vec3(0.0, 0.0, 0.0) # No light contribution at max depth
    end

    rec = hit(world, ray, Interval(0.001, RTWeekend.infinity))
    if rec !== nothing
        direction = rec.normal + random_unit_vector() # Lambertian reflection
        return 0.5 * ray_color(camera, Ray(rec.p, direction), world, depth - 1)
    end

    # Background gradient
    unit_direction = normalize(ray.direction)
    t = 0.5 * (unit_direction.y + 1.0)
    return (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0)
end

function linear_to_gamma(linear_component::Float64)
    if linear_component > 0
        return sqrt(linear_component)
    else
        return 0.0
    end
end

function write_color(file, pixel_color::Vec3, samples_per_pixel::Int)
    scale = 1.0 / samples_per_pixel
    r = pixel_color.x * scale
    g = pixel_color.y * scale
    b = pixel_color.z * scale

    # Apply gamma correction
    r = linear_to_gamma(r)
    g = linear_to_gamma(g)
    b = linear_to_gamma(b)

    # Clamp the color to [0, 0.999] and convert to [0, 255]
    ir = trunc(Int, 256 * my_clamp(r, 0.0, 0.999))
    ig = trunc(Int, 256 * my_clamp(g, 0.0, 0.999))
    ib = trunc(Int, 256 * my_clamp(b, 0.0, 0.999))

    println(file, ir, " ", ig, " ", ib)
end

function render(camera::Camera, world::HittableList, filename::String)
    open(filename, "w") do file
        println(file, "P3")
        println(file, camera.image_width, " ", camera.image_height)
        println(file, "255")

        for j in 0:camera.image_height-1
            for i in 0:camera.image_width-1
                pixel_color = Vec3(0.0, 0.0, 0.0)
                for s in 1:camera.samples_per_pixel
                    ray = get_ray(camera, i, j)
                    pixel_color += ray_color(camera, ray, world, 50) # Max depth = 50
                end
                write_color(file, pixel_color, camera.samples_per_pixel)
            end
        end
    end
    println("Rendering complete. Image saved as $filename")
end

# main.jl
function main()
    # World setup
    world = HittableList()
    add!(world, Sphere(Vec3(0.0, 0.0, -1.0), 0.5))
    add!(world, Sphere(Vec3(0.0, -100.5, -1.0), 100.0))

    # Camera setup
    cam = Camera(aspect_ratio=16.0 / 9.0, image_width=400, samples_per_pixel=100)

    # Render the scene
    render(cam, world, "images/ch9.5_image.ppm")
    ppm_to_jpg("images/ch9.5_image.ppm", "images/ch9.5_image.jpg")
end

main()
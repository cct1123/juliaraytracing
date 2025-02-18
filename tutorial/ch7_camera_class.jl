# main.jl
include("rtweekend.jl")
using .RTWeekend
include("hittable.jl")


# camera.jl
struct Camera
    aspect_ratio::Float64  # Ratio of image width over height
    image_width::Int       # Rendered image width in pixel count
    image_height::Int      # Rendered image height (calculated)
    center::Vec3           # Camera center
    pixel00_loc::Vec3      # Location of pixel (0, 0)
    pixel_delta_u::Vec3    # Offset to pixel to the right
    pixel_delta_v::Vec3    # Offset to pixel below
end

function Camera(; aspect_ratio=16.0 / 9.0, image_width=400)
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

    # Create and return the Camera object
    Camera(aspect_ratio, image_width, image_height, center, pixel00_loc, pixel_delta_u, pixel_delta_v)
end

function get_ray(camera::Camera, i::Int, j::Int)
    # Get a ray for the pixel at location (i, j)
    pixel_center = camera.pixel00_loc + (i * camera.pixel_delta_u) + (j * camera.pixel_delta_v)
    ray_direction = pixel_center - camera.center
    Ray(camera.center, ray_direction)
end

function ray_color(camera::Camera, ray::Ray, world::HittableList)
    rec = hit(world, ray, Interval(0.001, RTWeekend.infinity))
    if rec !== nothing
        return 0.5 * Vec3(rec.normal.x + 1.0, rec.normal.y + 1.0, rec.normal.z + 1.0)
    end

    # Background gradient
    unit_direction = normalize(ray.direction)
    t = 0.5 * (unit_direction.y + 1.0)
    return (1.0 - t) * Vec3(1.0, 1.0, 1.0) + t * Vec3(0.5, 0.7, 1.0)
end

function render(camera::Camera, world::HittableList, filename::String)
    open(filename, "w") do file
        println(file, "P3")
        println(file, camera.image_width, " ", camera.image_height)
        println(file, "255")

        for j in 0:camera.image_height-1
            for i in 0:camera.image_width-1
                ray = get_ray(camera, i, j)
                col = ray_color(camera, ray, world)
                
                ir = trunc(Int, 255.999 * col.x)
                ig = trunc(Int, 255.999 * col.y)
                ib = trunc(Int, 255.999 * col.z)
                
                println(file, ir, " ", ig, " ", ib)
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
    cam = Camera(aspect_ratio=16.0 / 9.0, image_width=400)

    # Render the scene
    render(cam, world, "images/ch7_image.ppm")
    ppm_to_jpg("images/ch7_image.ppm", "images/ch7_image.jpg")
end

main()
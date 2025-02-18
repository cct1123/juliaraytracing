struct Camera
    aspect_ratio::Float64      # Ratio of image width over height
    image_width::Int           # Rendered image width in pixel count
    image_height::Int          # Rendered image height (calculated)
    samples_per_pixel::Int     # Number of random samples per pixel
    max_depth::Int             # Maximum number of ray bounces
    vfov::Float64              # Vertical field of view (degrees)
    lookfrom::Vec3             # Camera position
    lookat::Vec3               # Point the camera is looking at
    vup::Vec3                  # Camera-relative "up" direction
    center::Vec3               # Camera center
    pixel00_loc::Vec3          # Location of pixel (0, 0)
    pixel_delta_u::Vec3        # Offset to pixel to the right
    pixel_delta_v::Vec3        # Offset to pixel below
    u::Vec3                    # Camera basis vector (right)
    v::Vec3                    # Camera basis vector (up)
    w::Vec3                    # Camera basis vector (opposite view direction)
    pixel_samples_scale::Float64 # Color scale factor for a sum of pixel samples
    defocus_angle::Float64     # Variation angle of rays through each pixel
    focus_dist::Float64        # Distance from camera to focus plane
    defocus_disk_u::Vec3       # Defocus disk horizontal radius
    defocus_disk_v::Vec3       # Defocus disk vertical radius
end

function Camera(;
    aspect_ratio=16.0 / 9.0,
    image_width=400,
    samples_per_pixel=100,
    max_depth=50,
    vfov=90.0,
    lookfrom=Vec3(0.0, 0.0, 0.0),
    lookat=Vec3(0.0, 0.0, -1.0),
    vup=Vec3(0.0, 1.0, 0.0),
    defocus_angle=0.0,  # Default: no defocus blur
    focus_dist=10.0     # Default: focus distance
)
    # Calculate image height and ensure it's at least 1
    image_height = max(1, trunc(Int, image_width / aspect_ratio))

    # Camera center
    center = lookfrom

    # Determine viewport dimensions
    theta = deg2rad(vfov)  # Convert FOV from degrees to radians
    h = tan(theta / 2)     # Half-height of the viewport
    viewport_height = 2 * h * focus_dist
    viewport_width = viewport_height * (image_width / image_height)

    # Calculate the u, v, w unit basis vectors for the camera coordinate frame
    w = normalize(lookfrom - lookat)  # Opposite view direction
    u = normalize(cross(vup, w))      # Right direction
    v = cross(w, u)                   # Up direction

    # Calculate the vectors across the horizontal and down the vertical viewport edges
    viewport_u = viewport_width * u    # Vector across viewport horizontal edge
    viewport_v = viewport_height * -v  # Vector down viewport vertical edge

    # Calculate the horizontal and vertical delta vectors from pixel to pixel
    pixel_delta_u = viewport_u / image_width
    pixel_delta_v = viewport_v / image_height

    # Calculate the location of the upper-left pixel
    viewport_upper_left = center - (focus_dist * w) - (viewport_u / 2) - (viewport_v / 2)
    pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v)

    # Calculate the defocus disk basis vectors
    defocus_radius = focus_dist * tan(deg2rad(defocus_angle / 2))
    defocus_disk_u = u * defocus_radius
    defocus_disk_v = v * defocus_radius

    # Color scale factor for a sum of pixel samples
    pixel_samples_scale = 1.0 / samples_per_pixel

    # Create and return the Camera object
    Camera(
        aspect_ratio,
        image_width,
        image_height,
        samples_per_pixel,
        max_depth,
        vfov,
        lookfrom,
        lookat,
        vup,
        center,
        pixel00_loc,
        pixel_delta_u,
        pixel_delta_v,
        u,
        v,
        w,
        pixel_samples_scale,
        defocus_angle,
        focus_dist,
        defocus_disk_u,
        defocus_disk_v
    )
end

function sample_square()
    # Returns a random point in the square surrounding a pixel
    Vec3(random_double(0.0,1.0) - 0.5, random_double(0.0,1.0) - 0.5, 0.0)
end


# function get_ray(camera::Camera, i::Int, j::Int)
#     # Get a ray for the pixel at location (i, j) with random sampling
#     pixel_sample = camera.pixel00_loc + (i + random_double(0.0,1.0) - 0.5) * camera.pixel_delta_u + (j + random_double(0.0,1.0) - 0.5) * camera.pixel_delta_v
#     ray_direction = pixel_sample - camera.center
#     Ray(camera.center, ray_direction)
# end

function get_ray(camera::Camera, i::Int, j::Int)
    # Get a random point on the defocus disk
    pixel_sample = camera.pixel00_loc + (i + random_double(0.0, 1.0) - 0.5) * camera.pixel_delta_u + (j + random_double(0.0, 1.0) - 0.5) * camera.pixel_delta_v
    ray_origin = camera.defocus_angle <= 0 ? camera.center : camera.center + (camera.defocus_disk_u * random_in_unit_disk().x) + (camera.defocus_disk_v * random_in_unit_disk().y)
    ray_direction = pixel_sample - ray_origin
    Ray(ray_origin, ray_direction)
end

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

function progress_bar(i, n)
    percent = floor(Int, (i / n) * 100)  # Calculate percentage
    bar_length = 100  # Length of the bar (in characters)
    progress = floor(Int, (i / n) * bar_length)  # Calculate progress
    # Create the bar string
    bar = "|" * repeat("█", progress) * repeat(" ", bar_length - progress) * "|"
    # Print the progress bar
    print("\r", "Rendering ", bar, " ", percent, "%")
    flush(stdout)  # Ensure the output is immediately printed
    if i == n
        println("") # move to the next line
    end
end


function render(camera::Camera, world::HittableList, filename::String)
    open(filename, "w") do file
        println(file, "P3")
        println(file, camera.image_width, " ", camera.image_height)
        println(file, "255")

        for j in 0:camera.image_height-1
            progress_bar(j, camera.image_height-1)
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
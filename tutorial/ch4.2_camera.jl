include("ch4.1_ray.jl")
include("ch3.1_color.jl")

# Function to calculate the color based on ray direction
function ray_color(r::Ray)
    a = 0.5 * (r.direction.y + 1.0)  # Calculate a value based on the y-component of the direction
    white = Color(1.0, 1.0, 1.0)  # White color
    sky_color = Color(0.5, 0.7, 1.0)  # Sky color

    # Blend between white and sky color based on 'a'
    return (1.0 - a) * white + a * sky_color
end

# Main rendering function

function render_scene(filename, image_width::Int)
    # Image setup
    aspect_ratio = 16.0 / 9.0
    image_height = round(Int, image_width / aspect_ratio)
    image_height = max(image_height, 1)

    # Camera setup
    focal_length = 1.0
    viewport_height = 2.0
    viewport_width = viewport_height * (Float64(image_width) / image_height)
    camera_center = Point3(0.0, 0.0, 0.0)

    # Viewport vectors
    viewport_u = Vec3(viewport_width, 0.0, 0.0)
    viewport_v = Vec3(0.0, -viewport_height, 0.0)

    # Pixel deltas
    pixel_delta_u = normalize(viewport_u) / Float64(image_width)
    pixel_delta_v = normalize(viewport_v) / Float64(image_height)

    # Upper left pixel
    viewport_upper_left = camera_center - Vec3(0.0, 0.0, focal_length) - viewport_u / 2.0 - viewport_v / 2.0
    pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v)

    open(filename, "w") do file
        # Render header
        println(file, "P3")
        println(file, "$image_width $image_height")
        println(file, "255")

        # Render image pixels
        for j in 1:image_height
            print("\rScanlines remaining: $(image_height - j) ")
            for i in 1:image_width
                pixel_center = pixel00_loc + Float64(i - 1) * pixel_delta_u + Float64(j - 1) * pixel_delta_v
                ray_direction = pixel_center - camera_center
                r = Ray(camera_center, ray_direction)

                pixel_color = ray_color(r)
                write_color(file, pixel_color)
            end
        end
    end
    println("\rDone.")
end

# Run the rendering function
filename = "images/render_scene.ppm"
render_scene(filename, 400)
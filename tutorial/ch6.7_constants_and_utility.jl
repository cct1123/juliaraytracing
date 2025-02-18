# main.jl
include("rtweekend.jl")
using .RTWeekend

# # Include other necessary files
# include("hittable.jl")
# include("sphere.jl")
# include("hittable_list.jl")
include("ch6.5_list_of_hittable.jl")

# Main rendering function
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
    ppmfile_name = "images/ch6.7_image.ppm"
    open(ppmfile_name, "w") do file
        println(file, "P3")
        println(file, image_width, " ", image_height)
        println(file, "255")
        
        for j in 0:image_height-1  # Iterate from top to bottom
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
    ppm_to_jpg(ppmfile_name, "images/ch6.7_image.jpg")
end

main()
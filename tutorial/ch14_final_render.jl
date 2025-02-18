# rtweekend.jl
include("rtweekend.jl")
using .RTWeekend

include("hittable.jl")

include("camera.jl")

# main.jl
function main()
    # World setup
    world = HittableList()

    # Ground plane (a large Lambertian sphere)
    ground_material = Lambertian(Vec3(0.5, 0.5, 0.5))
    add!(world, Sphere(Vec3(0.0, -1000.0, 0.0), 1000.0, ground_material))

    # Generate random small spheres
    for a in -11:10
        for b in -11:10
            choose_mat = random_double()
            center = Vec3(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double())

            # Ensure spheres don't overlap with the large spheres
            if norm(center - Vec3(4.0, 0.2, 0.0)) > 0.9
                if choose_mat < 0.8
                    # Diffuse material
                    albedo = Vec3(random_double(), random_double(), random_double()) * Vec3(random_double(), random_double(), random_double())
                    sphere_material = Lambertian(albedo)
                elseif choose_mat < 0.95
                    # Metal material
                    albedo = Vec3(random_double(0.5, 1.0), random_double(0.5, 1.0), random_double(0.5, 1.0))
                    fuzz = random_double(0.0, 0.5)
                    sphere_material = Metal(albedo, fuzz)
                else
                    # Glass material
                    sphere_material = Dielectric(1.5)
                end
                add!(world, Sphere(center, 0.2, sphere_material))
            end
        end
    end

    # Add three large spheres
    material1 = Dielectric(1.5)  # Glass sphere
    add!(world, Sphere(Vec3(0.0, 1.0, 0.0), 1.0, material1))

    material2 = Lambertian(Vec3(0.4, 0.2, 0.1))  # Diffuse sphere
    add!(world, Sphere(Vec3(-4.0, 1.0, 0.0), 1.0, material2))

    material3 = Metal(Vec3(0.7, 0.6, 0.5), 0.0)  # Metal sphere
    add!(world, Sphere(Vec3(4.0, 1.0, 0.0), 1.0, material3))

    # Camera setup
    cam = Camera(
        aspect_ratio=16.0 / 9.0,
        image_width=1200,
        samples_per_pixel=50,  # High sample count for a clean image
        max_depth=50,
        vfov=20.0,  # Narrow field of view for a zoomed-in effect
        lookfrom=Vec3(13.0, 2.0, 3.0),  # Camera position
        lookat=Vec3(0.0, 0.0, 0.0),     # Point the camera is looking at
        vup=Vec3(0.0, 1.0, 0.0),        # Camera "up" direction
        defocus_angle=0.6,              # Defocus blur angle
        focus_dist=10.0                 # Focus distance
    )

    # Render the scene
    render(cam, world, "images/final_render.ppm")
    ppm_to_jpg("images/final_render.ppm", "images/final_render.jpg")
end

main()
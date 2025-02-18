# rtweekend.jl
include("rtweekend.jl")
using .RTWeekend
# Helper function to check if a vector is near zero
function near_zero(v::Vec3)
    s = 1e-8
    return abs(v.x) < s && abs(v.y) < s && abs(v.z) < s
end
# "hittable.jl"
include("hittable.jl")



# camera.jl
include("camera.jl")



# Refraction function
function refract(uv::Vec3, n::Vec3, etai_over_etat::Float64)
    cos_theta = min(dot(-uv, n), 1.0)
    r_out_perp = etai_over_etat * (uv + cos_theta * n)
    r_out_parallel = -sqrt(abs(1.0 - length_squared(r_out_perp))) * n
    return r_out_perp + r_out_parallel
end

# Schlick approximation for reflectance
function reflectance(cosine::Float64, refraction_index::Float64)
    r0 = (1 - refraction_index) / (1 + refraction_index)
    r0 = r0^2
    return r0 + (1 - r0) * (1 - cosine)^5
end


# Dielectric material
struct Dielectric <: Material
    refraction_index::Float64  # Refractive index of the material
end

function scatter(material::Dielectric, ray_in::Ray, rec::HitRecord)
    attenuation = Vec3(1.0, 1.0, 1.0)  # Glass absorbs nothing
    refraction_ratio = rec.front_face ? (1.0 / material.refraction_index) : material.refraction_index

    unit_direction = normalize(ray_in.direction)
    cos_theta = min(dot(-unit_direction, rec.normal), 1.0)
    sin_theta = sqrt(1.0 - cos_theta^2)

    cannot_refract = refraction_ratio * sin_theta > 1.0
    direction = if cannot_refract || reflectance(cos_theta, refraction_ratio) > rand()
        reflect(unit_direction, rec.normal)  # Reflect
    else
        refract(unit_direction, rec.normal, refraction_ratio)  # Refract
    end

    scattered = Ray(rec.p, direction)
    return true, attenuation, scattered
end


# Reflection function
function reflect(v::Vec3, n::Vec3)
    return v - 2 * dot(v, n) * n
end

function main()
    # World setup
    world = HittableList()

    material_ground = Lambertian(Vec3(0.8, 0.8, 0.0))
    material_center = Lambertian(Vec3(0.1, 0.2, 0.5))
    material_left = Dielectric(1.5)  # Glass sphere
    material_bubble = Dielectric(1.0 / 1.5)  # Air bubble inside glass
    material_right = Metal(Vec3(0.8, 0.6, 0.2), 0.0)

    add!(world, Sphere(Vec3(0.0, -100.5, -1.0), 100.0, material_ground))
    add!(world, Sphere(Vec3(0.0, 0.0, -1.2), 0.5, material_center))
    add!(world, Sphere(Vec3(-1.0, 0.0, -1.0), 0.5, material_left))
    add!(world, Sphere(Vec3(-1.0, 0.0, -1.0), 0.4, material_bubble))  # Inner air bubble
    add!(world, Sphere(Vec3(1.0, 0.0, -1.0), 0.5, material_right))

    # Camera setup
    cam = Camera(aspect_ratio=16.0 / 9.0, image_width=400, samples_per_pixel=100)

    # Render the scene
    render(cam, world, "images/ch11_image.ppm")
    ppm_to_jpg("images/ch11_image.ppm", "images/ch11_image.jpg")
end

main()
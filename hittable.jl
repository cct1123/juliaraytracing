abstract type Material end
mutable struct HitRecord
    p::Vec3                  # Intersection point
    normal::Vec3             # Surface normal
    t::Float64               # Ray parameter at intersection
    front_face::Bool         # Whether the ray hits the front face
    material::Material       # Material at the intersection point
end

function set_face_normal!(rec::HitRecord, ray::Ray, outward_normal::Vec3)
    rec.front_face = dot(ray.direction, outward_normal) < 0
    rec.normal = rec.front_face ? outward_normal : -outward_normal
end


# Hittable: Abstract type for objects that can be hit by a ray
abstract type Hittable end

function hit(obj::Hittable, ray::Ray, ray_t::Interval)::Union{HitRecord, Nothing}
    error("hit() not implemented for this Hittable type")
end

# sphere hittable object
struct Sphere <: Hittable
    center::Vec3
    radius::Float64
    material::Material
end

function hit(sphere::Sphere, ray::Ray, ray_t::Interval)::Union{HitRecord, Nothing}
    oc = sphere.center - ray.origin
    a = dot(ray.direction, ray.direction)
    h = dot(ray.direction, oc)
    c = dot(oc, oc) - sphere.radius^2
    discriminant = h^2 - a * c

    if discriminant < 0
        return nothing  # No intersection
    end

    sqrtd = sqrt(discriminant)

    # Find the nearest root within the acceptable range
    root = (h - sqrtd) / a
    if !surrounds(ray_t, root)
        root = (h + sqrtd) / a
        if !surrounds(ray_t, root)
            return nothing  # No valid intersection
        end
    end

    # Compute intersection point and normal
    p = at(ray, root)
    outward_normal = (p - sphere.center) / sphere.radius
    rec = HitRecord(p, outward_normal, root, true, sphere.material)
    set_face_normal!(rec, ray, outward_normal)

    return rec
end


# HittableList: A list of hittable objects
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

function hit(list::HittableList, ray::Ray, ray_t::Interval)::Union{HitRecord, Nothing}
    closest_so_far = ray_t.max
    hit_anything = nothing

    for obj in list.objects
        temp_rec = hit(obj, ray, Interval(ray_t.min, closest_so_far))
        if temp_rec !== nothing
            closest_so_far = temp_rec.t
            hit_anything = temp_rec
        end
    end

    return hit_anything
end



# some scattering and reflection

function reflect(v::Vec3, n::Vec3)
    return v - 2 * dot(v, n) * n
end


# Scatter function for materials
function scatter(material::Material, ray_in::Ray, rec::HitRecord)
    # Default behavior: no scattering
    return false, Vec3(0.0, 0.0, 0.0), Ray(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, 0.0))
end

struct Lambertian <: Material
    albedo::Vec3  # Surface color (reflectance)
end

function scatter(material::Lambertian, ray_in::Ray, rec::HitRecord)
    scatter_direction = rec.normal + random_unit_vector()

    # Catch degenerate scatter direction
    if near_zero(scatter_direction)
        scatter_direction = rec.normal
    end

    ray_scattered = Ray(rec.p, scatter_direction)
    attenuation = material.albedo
    return true, attenuation, ray_scattered
end

# Metal material with fuzziness
struct Metal <: Material
    albedo::Vec3  # Surface color (reflectance)
    fuzz::Float64 # Fuzziness factor (0 = no fuzz, 1 = max fuzz)
end

function scatter(material::Metal, ray_in::Ray, rec::HitRecord)
    reflected = reflect(normalize(ray_in.direction), rec.normal)
    scattered_direction = reflected + material.fuzz * random_unit_vector()
    scattered = Ray(rec.p, scattered_direction)
    attenuation = material.albedo
    return (dot(scattered.direction, rec.normal) > 0), attenuation, scattered
end


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

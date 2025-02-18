function hit_sphere(center::Vec3, radius::Float64, ray::Ray)
    oc = center - ray.origin
    a = dot(ray.direction, ray.direction)  # Same as before
    h = dot(ray.direction, oc)            # Replace b with h = dot(direction, oc)
    c = dot(oc, oc) - radius^2            # Same as before
    discriminant = h^2 - a*c              # Simplified discriminant

    if discriminant < 0
        return -1.0  # No intersection
    else
        return (h - sqrt(discriminant)) / a  # Simplified t calculation
    end
end
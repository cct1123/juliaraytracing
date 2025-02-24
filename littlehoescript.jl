# Physical Ray Tracing Script Example
# ====================================
# This script demonstrates a simplified physical ray tracer.
# Each ray carries its position (origin), direction, amplitude, phase,
# and (optionally) wavelength. At a dielectric boundary, we use Fresnel
# coefficients to decide probabilistically whether the ray is reflected
# or transmitted. The ray is then updated in place.
#
# The goal is to capture interference effects at a detector plane by
# ensuring that each ray's amplitude and phase are adjusted according
# to physical principles.

using LinearAlgebra, Random, Roots

# Define a basic 3D vector type alias.
const Vec3 = Vector{Float64}

#----------------------------
# Structures and Transformations
#----------------------------

# Define a Frame for coordinate transformations (if needed)
struct Frame
    origin::Vec3
    orientation::Matrix{Float64}  # 3x3 rotation matrix
end

# A mutable Ray carries its state, so we can update it in place.
mutable struct Ray
    origin::Vec3         # Position in space
    direction::Vec3      # Normalized direction vector
    amplitude::Float64   # Field amplitude
    phase::Float64       # Phase (radians)
    wavelength::Float64  # Wavelength (meters)
end

# A mutable HitRecord to store information at an interaction point.
mutable struct HitRecord
    t::Float64      # Time (or distance) along the ray where hit occurs
    p::Vec3         # Intersection point
    normal::Vec3    # Surface normal at the hit
    aligned::Bool   # Whether the ray direction aligns with the surface normal
    amplitude::Float64  # Amplitude at hit point (if needed)
    phase::Float64      # Phase at hit point (if needed)
end

# Transform ray to a new frame (e.g., to simplify intersection tests)
function transform_to(ray::Ray, frame::Frame)
    rotT = transpose(frame.orientation)
    new_origin = rotT * (ray.origin - frame.origin)
    new_direction = rotT * ray.direction
    return Ray(new_origin, new_direction, ray.amplitude, ray.phase, ray.wavelength)
end

# Transform ray back to the original frame
function transform_from(ray::Ray, frame::Frame)
    new_origin = frame.orientation * ray.origin + frame.origin
    new_direction = frame.orientation * ray.direction
    return Ray(new_origin, new_direction, ray.amplitude, ray.phase, ray.wavelength)
end

#----------------------------
# Reflection, Refraction, and Fresnel Calculations
#----------------------------

# Reflect the incident vector about a surface normal.
function reflect(normal::Vec3, incident::Vec3)
    cosI = -dot(normal, incident)
    return incident + 2 * cosI * normal
end

# Refract the incident vector. Returns an invalid vector (NaN) if total internal reflection occurs.
const invalid_vector = [NaN, NaN, NaN]
function refract(normal::Vec3, incident::Vec3, n1::Float64, n2::Float64)
    n_ratio = n1 / n2
    cosI = -dot(normal, incident)
    sinT2 = n_ratio^2 * (1 - cosI^2)
    if sinT2 > 1.0
        return invalid_vector  # Total Internal Reflection (TIR)
    end
    cosT = sqrt(1 - sinT2)
    return n_ratio * incident + (n_ratio * cosI - cosT) * normal
end

# Fresnel reflectance using Schlick's approximation.
function fresnel_reflectance(normal::Vec3, incident::Vec3, n1::Float64, n2::Float64)
    r0 = ((n1 - n2) / (n1 + n2))^2
    cosX = -dot(normal, incident)
    # Adjust cosX for rays going from a denser medium.
    if n1 > n2
        n_ratio = n1 / n2
        sinT2 = n_ratio^2 * (1 - cosX^2)
        if sinT2 > 1.0
            return 1.0  # Total internal reflection
        end
        cosX = sqrt(1 - sinT2)
    end
    return r0 + (1 - r0) * (1 - cosX)^5
end

#----------------------------
# Example: Interaction at a Dielectric Boundary
#----------------------------

# This function simulates a ray's interaction at a dielectric boundary.
# It uses a random number to choose between reflection and transmission.
# The ray is updated in place based on the Fresnel coefficients.
function interact_at_boundary(ray::Ray, surface_normal::Vec3, n1::Float64, n2::Float64)
    # Compute Fresnel reflection coefficient.
    R = fresnel_reflectance(surface_normal, ray.direction, n1, n2)
    if rand() < R
        # Reflection: update ray direction, amplitude, and phase.
        new_direction = reflect(surface_normal, ray.direction)
        new_amplitude = ray.amplitude * R
        new_phase = ray.phase + π   # Reflection typically induces a π phase shift.
        println("Reflection chosen (R = $R)")
    else
        # Transmission: update ray direction using refraction.
        new_direction = refract(surface_normal, ray.direction, n1, n2)
        if new_direction == invalid_vector
            # Fall back to reflection if total internal reflection occurs.
            new_direction = reflect(surface_normal, ray.direction)
            new_amplitude = ray.amplitude * R
            new_phase = ray.phase + π
            println("Total internal reflection; using reflection")
        else
            new_amplitude = ray.amplitude * (1 - R)
            new_phase = ray.phase  # Transmission may have little or no phase shift.
            println("Transmission chosen (T = $(1 - R))")
        end
    end
    # For demonstration, we update the ray's origin to a new point along the ray.
    # In practice, this should be the actual intersection point.
    ray.origin += ray.direction  # Move ray forward (a placeholder update)
    ray.direction = new_direction
    ray.amplitude = new_amplitude
    ray.phase = new_phase
    return ray
end

#----------------------------
# Main Simulation Loop Example
#----------------------------

# Simulate a ray's path through multiple interactions.
function simulate_ray_path(ray::Ray, n1::Float64, n2::Float64, max_bounces::Int)
    for bounce in 1:max_bounces
        # In a complete simulation, you would first test for surface intersection.
        # Here, we assume each bounce happens at a fixed surface with a known normal.
        surface_normal = [0.0, 0.0, 1.0]  # Example: horizontal surface.
        interact_at_boundary(ray, surface_normal, n1, n2)
        println("Bounce $bounce: origin=$(ray.origin), direction=$(ray.direction), amplitude=$(ray.amplitude), phase=$(ray.phase)")
        # Terminate if the ray's amplitude falls below a threshold.
        if ray.amplitude < 1e-3
            break
        end
    end
end

#----------------------------
# Example Initialization and Run
#----------------------------

# Initialize a ray with:
# - Origin: starting at (0, 0, -5)
# - Direction: pointing towards the positive z-axis (normalized)
# - Amplitude: 1.0 (normalized energy)
# - Phase: 0.0 (initial phase)
# - Wavelength: 500 nm (500e-9 meters)
initial_ray = Ray([0.0, 0.0, -5.0], normalize([0.0, 0.0, 1.0]), 1.0, 0.0, 500e-9)

# Simulate the ray's path through up to 10 interactions at a dielectric boundary.
# n1 = 1.0 (e.g., air), n2 = 1.5 (e.g., glass)
simulate_ray_path(initial_ray, 1.0, 1.5, 10)

#----------------------------
# Detector Considerations
#----------------------------
# In a full simulation, you would trace millions of rays.
# At the detector plane, each ray's amplitude and phase contribute to the
# overall field as a complex sum:
#
#     E_total = Σ (amplitude * exp(im * phase))
#
# This interference calculation will reveal the intensity pattern at the detector.
# Using amplitude adjustment allows you to statistically capture the Fresnel effects
# without needing to create an exponential number of ray instances.

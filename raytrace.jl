module RayTrace # module begin ----------------------------
using LinearAlgebra
if abspath(PROGRAM_FILE) == @__FILE__
    println("ONLY for DEVELPMENT: Script $(@__FILE__) is running directly!")
    include("./MathModule.jl")  # Include only if necessary
    using .MathModule
    include("./GeometryModule.jl")
    using .GeometryModule
    include("./MaterialModule.jl")
    using .MaterialModule
    include("./SourceModule.jl")
    using .SourceModule
    include("./ObjectModule.jl")
    using .ObjectModule
    include("./DetectorModule.jl")
    using .DetectorModule
else
    using Main.MathModule
    using Main.GeometryModule
    using Main.MaterialModule
    using Main.SourceModule
    using Main.ObjectModule
    using Main.DetectorModule
end

mutable struct Setup
    sources:: Vector{Source}
    objects:: Vector{Object}
    detectors:: Vector{Detector}
end


mutable struct HitRecord
    t::Float64               # time of intersection
    p::Vec3                  # Intersection point
    normal::Vec3             # Surface normal
    aligned::Bool         # whether the ray direction align with surface normal
    material:: Material
    function HitRecord(t, p, normal, aligned, material=AIR)
        return new(t, p, normal, aligned, material)
    end
end


function hit(surface::Surface, ray::Ray; t_min::Float64=1e-6, t_max::Float64=1e6)::Union{HitRecord, Nothing}
    """
    assume the ray direction is normalized, 
    Arguments:
        surface: Surface with frame, shape function, bounds and border

    return: 
        t_hit: time required to travel from ray origin to hit point
        p_hit: hit point on the surface in the surface frame
    """
    # ray_original = deepcopy(ray) 

    # initiate the hit HitRecord
    hitrecord = nothing
    # transform ray into the surface frame
    transform_ray_to!(ray, surface.frame)
    # println("trying to hit a surface")
    # TODO: before searching for root, estimate the hit time and hit point first using the frame and bounds of the surface
    t_hit = find_intersection(surface.shape, ray.origin, ray.direction, t_min, t_max)
    if t_hit !== nothing
        # println("so I hitted something at t=", t_hit)
        p_hit = ray.origin + t_hit * ray.direction # in the surface frame
        # check whether the hit point is within the surface bounds
        xbounded = (surface.bounds[1][1] <= p_hit[1] <= surface.bounds[1][2])
        ybounded = (surface.bounds[2][1] <= p_hit[2] <= surface.bounds[2][2])
        bounded = xbounded && ybounded
        
        if bounded
            # further determin if the hit point is inside the border
            bounded = surface.border(p_hit[1:2]) < 0
            if bounded
                # find the surface normal at the hit point
                snormal = surface_normal(p_hit, surface)
                # println("Hiting surface: ", surface)
                # println("hit time at ", t_hit)
                # println("hit point at ", p_hit)
                # println("surface normal ", snormal)
                # determine if the ray is aligned with the surface normal, i.e. going from inside to outside
                aligned = dot(ray.direction, snormal) > 0

                hitrecord = HitRecord(t_hit, p_hit, snormal, aligned)
                # bring the p_hit and snormal back to the original ray frame
                hitrecord.p = transform_point_from(p_hit, surface.frame)
                hitrecord.normal = transform_direction_from(snormal, surface.frame)
            end
        end
    end
    # transform ray back to the original frame
    transform_ray_from!(ray, surface.frame)

    # println("orignial ray:", ray_original)
    # println("ray after hitting:", ray)
    return hitrecord
end

function hit(volume::Volume, ray::Ray; t_min=1e-6, t_max=1e6)::Union{HitRecord, Nothing}
    hitrecord = nothing
    hitrecord_surface = nothing
    t_hit = t_max
    # transform ray to the volume frame
    transform_ray_to!(ray, volume.frame)
    for surface in volume.surfaces
        hitrecord_surface = hit(surface, ray;t_min=t_min, t_max=t_max)
        if hitrecord_surface !== nothing
            if hitrecord_surface.t < t_hit
                t_hit = hitrecord_surface.t
                hitrecord = hitrecord_surface
            end
        end
    end
    if hitrecord !== nothing
        # transform hit point and surface normal back to the original frame
        hitrecord.p = transform_point_from(hitrecord.p, volume.frame)
        hitrecord.normal = transform_direction_from(hitrecord.normal, volume.frame)
    end
    # transform ray back to the original frame
    transform_ray_from!(ray, volume.frame)
    return hitrecord
end

function hit(object::Object, ray::Ray; t_min=1e-6, t_max=1e6)::Union{HitRecord, Nothing}
    hitrecord = hit(object.geometry, ray; t_min=t_min, t_max=t_max)
    if hitrecord !== nothing
        hitrecord.material = object.material
    end
    return hitrecord
end

function scatter(ray_in::Ray, p_hit::Vector{Float64}, snormal::Vector{Float64}, material::Dielectric)::Vector{Ray}
    # the normal is pointing from material to the outside
    aligned = dot(ray_in.direction, snormal) > 0
    # println("aligned?: ", aligned)
    n1 = aligned ? material.n : 1.0
    n2 = aligned ? 1.0 : material.n
    atten = aligned ? exp(-material.α*norm(p_hit-ray_in.origin)) : 1.0
    ray_in.amplitude = ray_in.amplitude*atten
    bigt, refracted, bigr, reflected = fresnel(-snormal, ray_in.direction, n1, n2)

    return [Ray(p_hit, refracted;amplitude=ray_in.amplitude*sqrt(bigt)), 
            Ray(p_hit, reflected;amplitude=ray_in.amplitude*sqrt(bigr))]
end

function scatter(ray_in::Ray, p_hit::Vector{Float64}, snormal::Vector{Float64}, material1::Dielectric, material2::Dielectric)::Vector{Ray}
    # the ray is hitting at a boundary between two objects
    # the normal is pointing from material1 to material2
    aligned = dot(ray_in.direction, snormal) > 0
    
    n1 = aligned ? material1.n : material2.n
    n2 = aligned ? material2.n : material1.n

    # println("Hi I hit the boundary of two objects")
    # println("Is it going from n1 to n2?", aligned)
    # println("n1 = $(n1), n2 = $(n2)")
    atten = aligned ? exp(-material1.α*norm(p_hit-ray_in.origin)) : exp(-material2.α*norm(p_hit-ray_in.origin))
    ray_in.amplitude = ray_in.amplitude*atten
    bigt, refracted, bigr, reflected = fresnel(-snormal, ray_in.direction, n1, n2)

    return [Ray(p_hit, refracted;amplitude=ray_in.amplitude*sqrt(bigt)), 
            Ray(p_hit, reflected;amplitude=ray_in.amplitude*sqrt(bigr))]
end



DELTA_T = 1E-9 # to count for the floating point error, which should be determined by the wavelength 

function propagate_ray!(objects::Vector{Object}, ray::Ray, hitrecord_shortest::Vector{Union{Nothing, HitRecord}}, tjnode::TrajectoryNode, depth::Int64=0; 
                        t_min::Float64=1e-6, t_max::Float64=1e6, amp_terminate=1e-6, depth_max::Int64=100)
    # WARNING: t_min should not be 0 otherwise the ray always hits the same surface!
    # to determine which object the ray hits first --------------------------------
    t_hit_shortest = t_max
    
    # intialize hitrecord_shortest
    for jj in 1:1:length(hitrecord_shortest)
        hitrecord_shortest[jj] = nothing
    end
    # println(hitrecord_shortest)
    # ray_original = deepcopy(ray) 
    for object in objects
        hitrecord = hit(object, ray; t_min=t_min, t_max=t_max)
        if hitrecord !== nothing
            if hitrecord.t < t_hit_shortest 
                t_hit_shortest = hitrecord.t
                hitrecord_shortest[1] = hitrecord
            elseif abs(hitrecord.t-t_hit_shortest) < DELTA_T # TODO the determined float might be slightly different by root finding
                # the ray might hit the same place but with a different object
                hitrecord_shortest[2] = hitrecord
            end
        end
    end
    # to determine which object the ray hits first --------------------------------
    if hitrecord_shortest[1] === nothing
        # the ray doesn't hit any shit
        return depth
    end
    # println("the ray hit at t=", hitrecord_shortest[1].t)
    rays_out = Ray[]
    if hitrecord_shortest[2] === nothing
        # only hit one object for sure
        # rays_out = scatter(ray, hitrecord_shortest[1])
        rays_out = scatter(ray, hitrecord_shortest[1].p, hitrecord_shortest[1].normal, hitrecord_shortest[1].material)
        # println("rays out = ", rays_out)
    else
        # might hit a boundary between two objects
        if abs(hitrecord_shortest[1].t-hitrecord_shortest[2].t) > DELTA_T
            # only hit one object at the closest, the overlapping boudary is not the closest
            # rays_out = scatter(ray, hitrecord_shortest[1])
            rays_out = scatter(ray, hitrecord_shortest[1].p, hitrecord_shortest[1].normal, hitrecord_shortest[1].material)
        else
            # hit the boundary of two objects at the closest
            # rays_out = scatter(ray, hitrecord_shortest[1], hitrecord_shortest[2])
            rays_out = scatter(ray, hitrecord_shortest[1].p, hitrecord_shortest[1].normal, hitrecord_shortest[1].material, hitrecord_shortest[2].material)
        end
    end
    
    for ray_new in rays_out
        # println("ray_new = ", ray_new)
        node_new = TrajectoryNode(ray_new, [])
        push!(tjnode.children, node_new)
        if (ray_new.amplitude > amp_terminate) && (depth < depth_max)
            propagate_ray!(objects, ray_new, hitrecord_shortest, node_new, depth + 1; amp_terminate=amp_terminate)
            # otherwise ray is terminated
        elseif depth >= depth_max
            # kill the ray if the depth is over depth_max
            ray_new.amplitude = 0.0
        end
    end
    return depth + 1
end



function detect!(detector::TrajectoryRecorder, trajectory::Trajectory)
    detector.count += 1
    if detector.count % detector.interval == 0
        push!(detector.trajectories, trajectory)
    end
end

function detect!(detector::Camera, trajectory::Trajectory)
    leafnodes = find_leaf_nodes(trajectory)
    for leafnode in leafnodes
        ray = leafnode.ray
        hitrecord = hit(detector.geometry, ray)
        if hitrecord !== nothing
            x, y, z =hitrecord.p
            @assert z ≈ 0.0
            idx_pixel_x = x // detector.pixel_size[1] 
            idx_pixel_y = y // detector.pixel_size[2]
            reminder_x = x % detector.pixel_size[1]
            reminder_y = y % detector.pixel_size[2]       
            cosθ = dot(hitrecord.normal, ray.direction)
            camera.pixels[idx_pixel_x+1, idx_pixel_y+1] += abs(ray.amplitude*cosθ)^2
            ray_new = Ray(hitrecord.p, ray.direction; amplitude = 0.0)
            node_new = TrajectoryNode(ray_new, [])
            push!(leafnode.children, node_new)
        end
    end
end

function detect!(detector::Objective, trajectory::Trajectory)
    leafnodes = find_leaf_nodes(trajectory)
    for leafnode in leafnodes
        ray = leafnode.ray
        hitrecord = hit(detector.surface, ray)
        if hitrecord !== nothing
            if abs(dot(hitrecord.normal, ray.direction)) > detector.NA
                detector.collection += abs(ray.amplitude)^2
                ray_new = Ray(hitrecord.p, ray.direction; amplitude = 0.0)
                node_new = TrajectoryNode(ray_new, [])
                push!(leafnode.children, node_new)
            end
            
        end
    end
end

export Setup, HitRecord, hit, scatter, propagate_ray!, DELTA_T, detect!
end # module end ------------------------------------------------
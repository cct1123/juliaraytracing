using GLMakie, LinearAlgebra

# Define the Ray struct
struct Ray
    origin::Vector{Float64}  # 3D point (origin)
    direction::Vector{Float64}  # 3D direction (unit vector)
end

# Function to generate and plot rays
function draw_rays(rays::Vector{Ray}; ray_length=1.0)
    fig = Figure(resolution=(800, 600))
    ax = Axis3(fig[1, 1], title="3D Rays Visualization", aspect=:data)

    for ray in rays
        start = ray.origin
        stop = start + normalize(ray.direction) * ray_length
        lines!(ax, [start[1], stop[1]], [start[2], stop[2]], [start[3], stop[3]], linewidth=2)
        scatter!(ax, [start[1]], [start[2]], [start[3]], color=:red, markersize=10)  # Mark origin
    end

    display(fig)
end

# Generate some rays from a point source
function generate_rays(center::Vector{Float64}, num_rays::Int)
    rays = Vector{Ray}()
    for _ in 1:num_rays
        dir = normalize(randn(3))  # Random unit direction
        push!(rays, Ray(center, dir))
    end
    return rays
end

# Example: A point source at (0,0,0) emitting rays in random directions
rays = generate_rays([0.0, 0.0, 0.0], 1000)
draw_rays(rays)
sleep(10)

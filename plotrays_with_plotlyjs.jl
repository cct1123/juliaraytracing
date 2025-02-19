using PlotlyJS, LinearAlgebra

# Define the Ray struct
struct Ray
    origin::Vector{Float64}  # 3D point (origin)
    direction::Vector{Float64}  # 3D direction (unit vector)
end

# Function to generate and plot rays
function draw_rays(rays::Vector{Ray}; ray_length=1.0)
    # Create lists for the origin points and ray ends
    origins_x, origins_y, origins_z = Float64[], Float64[], Float64[]
    ends_x, ends_y, ends_z = Float64[], Float64[], Float64[]
    traces = GenericTrace[]
    for ray in rays
        start = ray.origin
        stop = start + normalize(ray.direction) * ray_length
        trace_rays = scatter3d(
        x=[start[1], stop[1]], y=[start[2], stop[2]], z=[start[3], stop[3]],
        mode="lines",
        name="Ray Directions"
        )
        push!(traces, trace_rays)
    end

    layout = Layout(
        title="3D Rays Visualization",
        scene=attr(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            aspectmode="data"
        ),
        showlegend=false
    )
    
    # Create the plot with both origin and ray lines
    html = plot(traces, layout)
    savefig(html, "rays_plot.html")  # Save the plot as an HTML file

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

# Example: A point source at (0, 0, 0) emitting rays in random directions
rays = generate_rays([0.0, 0.0, 0.0], 1000)
draw_rays(rays)
sleep(10)

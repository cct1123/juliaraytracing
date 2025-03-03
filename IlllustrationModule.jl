module IlllustrationModule # module begin ============
 

using PlotlyJS

LAYOUT_DEFAULT = PlotlyJS.Layout(
    title="Plot",
    width=800,   # Set width in pixels
    height=600,  # Set height in pixels
    scene=attr(
        xaxis=attr(visible=false, showgrid=false, zeroline=false),
        yaxis=attr(visible=false, showgrid=false, zeroline=false),
        zaxis=attr(visible=false, showgrid=false, zeroline=false),
        bgcolor="rgba(0,0,0,0)",  # Transparent background
        aspectmode="data"
    ),
    paper_bgcolor="rgba(0,0,0,0)",  # Transparent outer background
    showlegend=false
    )

# Function to generate and plot rays with arrows
function draw_rays(rays::Vector{Ray}, ray_length::Float64=0.5, arrow_scale::Float64=0.3, color::String="random")
    traces = GenericTrace[]

    color = color == "random" ? "rgb($(rand(0:255)), $(rand(0:255)), $(rand(0:255)))" : color
    # Line traces for the rays
    for ray in rays
        start = ray.origin
        stop = start + normalize(ray.direction) * ray_length
        trace_rays = scatter3d(
            x=[start[1], stop[1]], y=[start[2], stop[2]], z=[start[3], stop[3]],
            mode="lines",
            line=attr(color=color,width=3)
        )
        push!(traces, trace_rays)
    end

    # Arrow traces using cones
    origins = hcat([ray.origin + normalize(ray.direction) * ray_length for ray in rays]...)  # Arrow positions
    directions = hcat([normalize(ray.direction) for ray in rays]...)  # Arrow directions
    if length(origins) > 0
        trace_arrows = cone(
            x=origins[1, :], y=origins[2, :], z=origins[3, :],  # Positions
            u=directions[1, :], v=directions[2, :], w=directions[3, :],  # Directions
            sizemode="raw",
            sizeref=arrow_scale,  # Controls arrow size
            anchor="tail",
            showscale=false,
            colorscale=[(0, color), (1, color)],
            opacity=0.8
        )
        push!(traces, trace_arrows)
    end

    return traces
end

function draw_surface(surface:: Surface, nx::Int=100, ny::Int=100, color="random", display::Bool=false)
    # plot the surface
    xarray = range(surface.bounds[1][1], stop=surface.bounds[1][2], length=nx)
    yarray = range(surface.bounds[2][1], stop=surface.bounds[2][2], length=ny)
    xmesh = zeros(nx, ny)
    ymesh = zeros(nx, ny)
    zmesh = zeros(nx, ny)
    for ii in 1:1:nx, jj in 1:1:ny
    
        bounded = surface.border([xarray[ii], yarray[jj]]) < 0
        xval = xarray[ii]
        yval = yarray[jj]
        zval = bounded ? -surface.shape([xarray[ii], yarray[jj], 0.0]) : NaN

        xval, yval, zval = surface.frame.orientation * [xval, yval, zval]
        xmesh[ii, jj] = xval + surface.frame.origin[1]
        ymesh[ii, jj] = yval + surface.frame.origin[2]
        zmesh[ii, jj] = zval + surface.frame.origin[3]
    end
    # Generate a random color (RGB format)
    if color=="random"
        color = string("rgb(", round(255 * rand()), ",", round(255 * rand()), ",", round(255 *  rand()), ")")
    end;

    surface_data = PlotlyJS.surface(
        x=xmesh,
        y=ymesh,
        z=zmesh,
        surfacecolor=zmesh*0.0,
        colorscale=[(0, color), (1, color)],  # Fix the color scale to the same random color
        shading="smooth",  # Smoothing the shading on the surface
        showscale=false,
        opacity=0.8
        )

    if display
        layout = LAYOUT_DEFAULT
        # Display the plot
        fig = plot([surface_data] , layout);
        display(fig)   
        return NaN  
    end
    return [surface_data]
end


function draw_surface_adaptive(surface::Surface; base_density=50, refinement_factor=3, color="random", display::Bool=false)
    x_range = surface.bounds[1]
    y_range = surface.bounds[2]

    # Initial coarse sampling
    x_coarse = range(x_range[1], stop=x_range[2], length=base_density)
    y_coarse = range(y_range[1], stop=y_range[2], length=base_density)

    x_points = Float64[]
    y_points = Float64[]

    for i in 1:length(x_coarse)-1
        for j in 1:length(y_coarse)-1
            x1, x2 = x_coarse[i], x_coarse[i+1]
            y1, y2 = y_coarse[j], y_coarse[j+1]

            z1 = surface.shape([x1, y1, 0.0])
            z2 = surface.shape([x2, y2, 0.0])

            dz = abs(z2 - z1)  # Height difference (morphology change)

            # More points where dz is large
            num_subdivisions = max(2, round(Int, dz * refinement_factor))
            x_sub = range(x1, stop=x2, length=num_subdivisions)
            y_sub = range(y1, stop=y2, length=num_subdivisions)

            append!(x_points, x_sub)
            append!(y_points, y_sub)
        end
    end

    xmesh = zeros(length(x_points))
    ymesh = zeros(length(y_points))
    zmesh = zeros(length(x_points))

    for i in 1:length(x_points)
        xval, yval = x_points[i], y_points[i]
        bounded = surface.border([xval, yval]) < 0
        zval = bounded ? -surface.shape([xval, yval, 0.0]) : NaN

        xval, yval, zval = surface.frame.orientation * [xval, yval, zval]
        xmesh[i] = xval + surface.frame.origin[1]
        ymesh[i] = yval + surface.frame.origin[2]
        zmesh[i] = zval + surface.frame.origin[3]
    end

    # Generate random color
    if color == "random"
        color = string("rgb(", round(255 * rand()), ",", round(255 * rand()), ",", round(255 * rand()), ")")
    end

    surface_data = PlotlyJS.scatter3d(
        x=xmesh, y=ymesh, z=zmesh,
        mode="markers",
        marker=Dict(:size => 2, :color => color, :opacity => 0.8)
    )

    if display
        layout = LAYOUT_DEFAULT
        fig = plot([surface_data], layout)
        display(fig)
        return NaN
    end
    return [surface_data]
end


function draw_snormals(surface:: Surface, num:: Int,  ray_length::Float64=0.5, arrow_scale::Float64=0.3, color="random")
    xarray = range(surface.bounds[1][1], stop=surface.bounds[1][2], length=nx)
    yarray = range(surface.bounds[2][1], stop=surface.bounds[2][2], length=ny)

    rays_snorm = Ray[]
    for ii in 1:1:num
        bounded = false
        while ~bounded
            xx = rand(surface.bounds[1][1]:surface.bounds[1][2])
            yy = rand(surface.bounds[2][1]:surface.bounds[2][2])
            bounded = surface.border([xx, yy]) < 0
        end
        zz = -surface.shape([xx,yy,0.0])
        snorm = surface_normal([xx, yy, zz], surface)
        snorm = surface.frame.orientation * snorm
        xx, yy, zz = surface.frame.orientation * [xx, yy, zz] + surface.frame.origin
        push!(rays_snorm, Ray([xx, yy, zz], snorm))
    end
    return draw_rays(rays_snorm, ray_length, arrow_scale, color)
end

function draw_object(
    obj::Object, 
    num_sx::Int, 
    num_sy::Int, 
    num_sn::Int,
    ;
    drawnormals::Bool=true, 
    arrow_scale::Float64=1.0,
    color_sf::String="random", 
    color_sn::String="random",
    display::Bool=true)
    objsize = sum([upper - lower for (lower, upper) in obj.bounds])/3.0
    arrsize = objsize/20.0*arrow_scale
    color_sf = color_sf == "random" ? "rgb($(rand(0:255)), $(rand(0:255)), $(rand(0:255)))" : color_sf
    color_sn = color_sn == "random" ? "rgb($(rand(0:255)), $(rand(0:255)), $(rand(0:255)))" : color_sn
    plottraces = GenericTrace[]
    for ss in obj.surfaces
        # rotate and shift the surface by the object's coordinates first
        frame_sf_lab = Frame(obj.frame.origin+ss.frame.origin, obj.frame.orientation*ss.frame.orientation)
        ss_tolab = ImplicitSurface(
            frame_sf_lab, 
            ss.shape, 
            ss.bounds, 
            ss.border
        ) # only work ImplicitSurface for now TODO: generalize it

        append!(plottraces, vcat(
            draw_surface(ss_tolab, num_sx, num_sy, color_sf), 
            draw_snormals(ss_tolab, num_sn, 0.5*arrsize, 0.3*arrsize, color_sn)
            )) 
    end

    if display
        layout = LAYOUT_DEFAULT
        # Display the plot
        fig = PlotlyJS.plot(plottraces , layout);
        PlotlyJS.display(fig)   
        return NaN
    else
        return plottraces
    end
end


# Function to extract ray segments for plotting
function extract_ray_segments(node::TrajectoryNode, segments::Vector{Tuple{Point, Point, Float64}})
    ray = node.ray
    if !isempty(node.children)
        # Use the origin of the first child as the endpoint
        endpoint = node.children[1].ray.origin
        push!(segments, (ray.origin, endpoint, ray.amplitude))
        for child in node.children
            extract_ray_segments(child, segments)
        end
    end
end

# Function to normalize amplitudes to [0,1] for opacity scaling
function normalize_amplitudes(segments)
    if length(segments) == 0
        return segments
    else
        max_amp = maximum(abs(seg[3]) for seg in segments)
        min_amp = minimum(abs(seg[3]) for seg in segments)
        return [(seg[1], seg[2], (seg[3] - min_amp) / (max_amp - min_amp + 1e-6)) for seg in segments]
    end
end

function set_opacity_in_rgba(rgba_string::String, opacity::Float64)
    # Remove "rgba()" and split the string by commas
    rgb_values = match(r"rgba\((\d+),\s*(\d+),\s*(\d+),\s*(\d*\.?\d+)\)", rgba_string).captures

    # Extract RGB values and convert to integers
    r = parse(Int, rgb_values[1])
    g = parse(Int, rgb_values[2])
    b = parse(Int, rgb_values[3])

    # Construct the new RGBA string with user-defined opacity
    new_rgba = "rgba($r, $g, $b, $opacity)"
    
    return new_rgba
end

RAYOPACITY_MIN = 0.05 # Minimum visibility
RAYOPACITY_MAX = 0.8 # Fully visible
# Function to convert amplitude to RGBA color (blue with variable opacity)
function amplitude_to_rgba(amplitude::Float64; color::String="rgba(255, 255, 0, 1)")
    alpha = RAYOPACITY_MIN + (RAYOPACITY_MAX - RAYOPACITY_MIN) * amplitude
    return set_opacity_in_rgba(color, alpha)
end

# Function to plot the trajectory using Plotly
function draw_trajectory(trajectory::Trajectory; color="rgba(255, 255, 0, 1)", sizemode="", ray_width=3.0, arrow_scale=1.0, draw_cone::Bool=true, display::Bool=true)    
    # sizemode:: String = "scale" or "fixed
    segments = Tuple{Point, Point, Float64}[]
    extract_ray_segments(trajectory.root, segments)

    # Normalize amplitude for opacity scaling
    segments = normalize_amplitudes(segments)

    traces = GenericTrace[]
    linewidth = ray_width * 1.0
    sizeref = arrow_scale * 50.0
    # Loop through each segment (ray) to plot rays and cones (arrows)
    for (i, (start_point, end_point, norm_amplitude)) in enumerate(segments)
        rgba_color = amplitude_to_rgba(norm_amplitude; color=color)  # Convert amplitude to RGBA color

        # direction vector from the ray origin to the endpoint
        direction = [end_point[1] - start_point[1], end_point[2] - start_point[2], end_point[3] - start_point[3]]
        magnitude = norm(direction)
        direction /= magnitude
        if sizemode == "scale"
            linewidth = ray_width * magnitude / 500.0 
            sizeref = arrow_scale * magnitude / 10.0
        end

        # Ray trace
        push!(traces, scatter3d(
            x=[start_point[1], end_point[1]], 
            y=[start_point[2], end_point[2]], 
            z=[start_point[3], end_point[3]], 
            mode="lines",
            line=attr(width=linewidth, color=rgba_color),  # Apply RGBA color
            showlegend=false,
        ))
        if draw_cone
            # Add cone at the endpoint of the ray (arrow visualization)
            # Cone trace, head at the endpoint
            cone_trace = cone(
                x=[end_point[1]], y=[end_point[2]], z=[end_point[3]],  # Head of the cone at the endpoint
                u=[direction[1]], v=[direction[2]], w=[direction[3]],  # Direction vector (from the origin to the endpoint)
                sizemode="raw",
                sizeref=sizeref,  # Controls arrow size
                anchor="tip",  # Position the cone head at the endpoint
                showscale=false,
                colorscale=[(0, rgba_color), (1, rgba_color)],  # Set the same color as the ray
            )
            
            push!(traces, cone_trace)
        end
    end


    if display
        layout = LAYOUT_DEFAULT
        layout.paper_bgcolor = "rgba(255,255,255,1)"
        # Display the plot
        fig = PlotlyJS.plot(traces , layout);
        PlotlyJS.display(fig)   
        return NaN
    else
        return traces
    end
end
   
end # module end ============
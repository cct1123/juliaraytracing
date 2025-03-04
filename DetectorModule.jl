module DetectorModule # module begin ----------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    println("ONLY for DEVELPMENT: Script $(@__FILE__) is running directly!")
    include("./MathModule.jl")  # Include only if necessary
    using .MathModule
    include("./GeometryModule.jl")
    using .GeometryModule
    include("./SourceModule.jl")
    using .SourceModule
else
    # include("./MathModule.jl")  # Include only if necessary
    using Main.MathModule
    using Main.GeometryModule: Volume, Surface, ImplicitSurface, make_celestial, make_circleplane
    using Main.SourceModule: Ray
end

mutable struct TrajectoryNode
    ray::Ray
    children::Vector{TrajectoryNode}
    # function TrajectoryNode(ray::Ray, children::Vector{TrajectoryNode} = [])
    #     new(ray, children)
    # end
end

mutable struct Trajectory
    root::TrajectoryNode
end


abstract type Detector end

mutable struct TrajectoryRecorder <: Detector
    interval:: Int64 # record trajectory every interval
    count::Int64
    trajectories:: Vector{Trajectory}
    function TrajectoryRecorder(interval::Int64, count::Int64=0, trajectories::Vector{Trajectory}=Trajectory[])
        new(interval, count, trajectories)
    end
end

reset_trajectoryrecorder!(tr::TrajectoryRecorder) = tr.count = 0

# mutable struct Counter <: Detector
#     count:: Int64
#     geometry:: Union{Volume, Surface} 
# end

mutable struct Camera <: Detector
    # a 2D areal detector with a border defined by a implicit line function
    frame:: Frame # the coordinate frame of the camera, not the image frame!!
    pixel_number:: Tuple{Int64, Int64}
    pixel_size:: Tuple{Float64, Float64}
    pixels:: Matrix{Float64} # record the intensity of each pixel
    surface:: Surface
    function Camera(
            frame::Frame, 
            pixel_number::Tuple{Int64, Int64}, 
            pixel_size::Tuple{Float64, Float64}, 
            # pixels::Union{Matrix{Float64}, Nothing}=nothing, 
            surface::Union{Surface, Nothing}=nothing
        )
        if surface === nothing 
            surface = make_plane(pixel_number[1]*pixel_size[1], pixel_number[2]*pixel_size[2], frame.origin, frame.orientation)
        end
        new(frame, pixel_number, pixel_size, zeros(pixel_number[1], pixel_number[2]), surface)
    end
end

function make_sCMOS(
    center::Vec3=[0.0, 0.0, 0.0],
    orientation::Matrix{Float64}=rm_eye
)
    # ref https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=11421
    return Camera(Frame(center, orientation), (1920, 1080), (5.04, 5.04))
end

function reset_pixels!(camera::Camera)
    camera.pixels = zeros(camera.pixel_number[1], camera.pixel_number[2])
end

mutable struct Objective <: Detector
    # a 2D areal detector with a border defined by a implicit line function
    frame:: Frame
    NA:: Float64 # numerical aperture
    WD:: Float64 # working distance [um]
    collection:: Float64 # total light intensity collected by the object
    surface:: Surface
    function Objective(frame::Frame, NA::Float64, WD::Float64, collection::Float64=0.0, surface::Union{Surface, Nothing}=nothing)
        if surface === nothing
            radius = WD*tan(asin(NA))
            surface = make_circleplane(radius, frame.origin, frame.orientation)
        end
        new(frame, NA, WD, collection, surface)
    end
end

function reset_collection!(obj::Objective)
    obj.collection = 0.0
end

function photon_collection!(obj::Objective, morephotons::Float64)
    obj.collection += morephotons
end


mutable struct Celestial <: Detector
    # record the ray hitting the celestial sphere
    center:: Point
    radius:: Float64

    geometry:: Union{Volume, Surface} 

    function Celestial(center::Point, radius::Float64, geometry::Union{Volume, Surface, Nothing}=nothing)
        if geometry === nothing
            geometry = make_celestial(center, radius)
        end
        return new(center, radius, geometry)
    end
end

function add_child!(parent::TrajectoryNode, child::Ray)
    push!(parent.children, TrajectoryNode(child, []))
end



function find_leaf_nodes(node::TrajectoryNode, leaves::Vector{TrajectoryNode} = [])
    if isempty(node.children)
        push!(leaves, node)
    else
        for child in node.children
            find_leaf_nodes(child, leaves)
        end
    end
    # !! WARNING !! modifying the leaf_nodes vector will directly affect the trajectory 
    #               because leaf_nodes contains references to the actual TrajectoryNode objects in the tree.
    return leaves
end

function find_leaf_nodes(trajectory::Trajectory)
    leaves = TrajectoryNode[]
    find_leaf_nodes(trajectory.root, leaves)
    # !! WARNING !! modifying the leaf_nodes vector will directly affect the trajectory 
    #               because leaf_nodes contains references to the actual TrajectoryNode objects in the tree.
    return leaves
end

export TrajectoryNode, Trajectory, Detector, TrajectoryRecorder, Counter, Camera, Objective, Celestial, 
    find_leaf_nodes;

end # module end ----------------------------
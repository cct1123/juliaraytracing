
abstract type Detector end

mutable struct Counter <: Detector
    geometry:: Union{Volume, Surface} 
    NA:: Float64
    count:: Int64
end


mutable struct Camera <: Detector
    # a 2D areal detector with a border defined by a implicit line function
    geometry:: Union{Volume, Surface} 
    NA:: Float64
    bounds:: Vector{Tuple{Float64}}
    border:: Function
end


mutable struct TrajectoryNode
    ray::Ray
    children::Vector{TrajectoryNode}
end

mutable struct Trajectory
    root::TrajectoryNode
end

function add_child!(parent::TrajectoryNode, child::Ray)
    push!(parent.children, TrajectoryNode(child, []))
end


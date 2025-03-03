module ObjectModule # module begin ----------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    println("ONLY for DEVELPMENT: Script $(@__FILE__) is running directly!")
    include("./MathModule.jl")
    using .MathModule
    include("./GeometryModule.jl")
    using .GeometryModule: Volume, Surface
    include("./MaterialModule.jl")
    using .MaterialModule: Material
else
    using Main.GeometryModule: Volume, Surface
    using Main.MaterialModule: Material
end

struct Object
    geometry:: Union{Volume, Surface}
    material:: Material
end

export Object

end # module end --------------------------------------
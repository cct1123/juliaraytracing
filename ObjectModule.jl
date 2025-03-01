
abstract type Object end # object that can be hit by rays

struct Hiitable <: Object
    geometry:: Union{Volume, Surface}
    material:: Material
end


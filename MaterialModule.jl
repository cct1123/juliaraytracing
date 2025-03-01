#  Define the structure for material and objects ==========================================================================
abstract type Material end

struct Dielectric <: Material
    n:: Float64 # refractive index
    α:: Float64 # absorption coefficient, [um^-1]
    function Dielectric(n::Float64, α::Float64=0.0)
        return new(n, α)
    end
end
# define some default materials
AIR = Dielectric(1.0)

struct Metal <: Material
    # fuzz:: Float64
    # albedo:: Vec3
end
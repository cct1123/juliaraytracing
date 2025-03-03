
module MaterialModule # module begin----------------------
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


n_diamond = 2.4
diamond = Dielectric(n_diamond)
diamond_ib_heavy = Dielectric(n_diamond, 50.0E-4)  # assumed ~50/cm, https://doi.org/10.1098/rsta.2022.0314
diamond_ib_mild = Dielectric(n_diamond, 10.0E-4)  # assumed ~50/cm, https://doi.org/10.1098/rsta.2022.0314
diamond_iia = Dielectric(n_diamond, 0.1E-4)  # assumed ~50/cm, https://doi.org/10.1098/rsta.2022.0314
air = Dielectric(1.0)
n_PDMS = 1.4
pdms = Dielectric(n_PDMS)

end # module end ============   
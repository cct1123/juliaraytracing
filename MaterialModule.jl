
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
struct Metal <: Material
    σ:: Float64 # conductivity
end

struct SheetMetal <: Material
    σ:: Float64 # conductivity
    δ:: Float64 # skin depth
    thickness:: Float64
end

# some common materials
const air = Dielectric(1.0)
const AIR = Dielectric(1.0)
const n_diamond = 2.4
const diamond = Dielectric(n_diamond)
const diamond_ib_heavy = Dielectric(n_diamond, 50.0E-4)  # assumed ~50/cm, https://doi.org/10.1098/rsta.2022.0314
const diamond_ib_mild = Dielectric(n_diamond, 10.0E-4)  # assumed ~10/cm, https://doi.org/10.1098/rsta.2022.0314
const diamond_iia = Dielectric(n_diamond, 0.1E-4)  # assumed ~50/cm, https://doi.org/10.1098/rsta.2022.0314
const n_PDMS = 1.4
const pdms = Dielectric(n_PDMS)

export Material, Dielectric, Metal, SheetMetal, air, AIR, diamond, diamond_ib_heavy, diamond_ib_mild, diamond_iia, pdms

end # module end ============   
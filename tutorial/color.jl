include("vec3.jl")
const Color = Vec3
function write_color(out::IO, pixel_color::Color)
    r = pixel_color.x
    g = pixel_color.y
    b = pixel_color.z

    rbyte = Int(round(255.999 * r))
    gbyte = Int(round(255.999 * g))
    bbyte = Int(round(255.999 * b))

    println(out, "$rbyte $gbyte $bbyte")
end

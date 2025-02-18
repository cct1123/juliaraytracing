# interval.jl
struct Interval
    min::Float64
    max::Float64
end

# Default interval is empty
Interval() = Interval(Inf64, -Inf64)

# Check if the interval contains a value
function contains(interval::Interval, x::Float64)
    return interval.min <= x <= interval.max
end

# Check if the interval surrounds a value
function surrounds(interval::Interval, x::Float64)
    return interval.min < x < interval.max
end

# Return the size of the interval
function size(interval::Interval)
    return interval.max - interval.min
end

# Define empty and universe intervals
const empty_interval = Interval(Inf64, -Inf64)
const universe_interval = Interval(-Inf64, Inf64)

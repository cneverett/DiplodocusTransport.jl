function nextpow2_float(x::Real)
    x > 0 || throw(DomainError(x, "Input must be positive to find a power of 2."))
    
    s = significand(x)
    e = exponent(x)
    
    # If it's already an exact power of 2, significand(x) will be exactly 1.0
    return s == 1.0 ? ldexp(1.0, e) : ldexp(1.0, e + 1)
end
function flags = isInteger( x )

    flags = round(x) == x | isinf(x);
    
end
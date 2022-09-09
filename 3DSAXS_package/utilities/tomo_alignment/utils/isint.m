function val = isint(x)
    val = all( x(:) - floor(x(:)) == 0);
end
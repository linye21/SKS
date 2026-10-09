function done = timelimit(tspan, y, type, kf, kr, x) %#ok
done = 0; 
if toc > x(end) 
    done = 1; 
end
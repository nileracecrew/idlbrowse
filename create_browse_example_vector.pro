pro create_browse_example_vector
    compile_opt idl2, hidden
    
    nt = 1000
    nx = 50
    ny = 60
    
    xpos = findgen(nx)/nx*10. - 5.
    ypos = findgen(ny)/ny*10. - 5.

    a = 2.
    b0 = 3.
    dataset = fltarr(nt, nx, ny, 2)
    for y = 0, ny - 1 do begin
        for x = 0, nx - 1 do begin
            r = sqrt(xpos[x]^2 + ypos[y]^2)
            theta = atan(ypos[y], xpos[x])
            if r le a then begin
                dataset[*, x, y, 0] = b0 * r / a * sin(theta)
                dataset[*, x, y, 1] = -b0 * r / a * cos(theta) 
            endif else begin
                dataset[*, x, y, 0] = b0 * a / r * sin(theta)
                dataset[*, x, y, 1] = -b0 * a / r * cos(theta)
            endelse
        endfor
    endfor

    for t = 0, nt - 1 do begin
        dataset[t, *, *, *] *= cos(4.*2.*!pi * t / nt)
    endfor

    t = findgen(nt)/nt*2*!pi
    
    (scope_varfetch('dataset', level=-1, /enter)) = dataset
    (scope_varfetch('t', level=-1, /enter)) = t
    (scope_varfetch('x', level=-1, /enter)) = xpos
    (scope_varfetch('y', level=-1, /enter)) = ypos

end
pro create_browse_example
    compile_opt idl2, hidden
    
    nt = 1000
    nx = 50
    ny = 60
    
    xpos = findgen(nx)/nx*10. - 5.
    ypos = findgen(ny)/ny*10. - 5.
    bz = 13.0152	
    a = 6.
    r0 = 6.3802/bz*a
    r1 = 9.7610/bz*a
    m = 3
    n = 3
    r = fltarr(nx, ny)
    theta = fltarr(nx, ny)
    z = fltarr(nx, ny)
    for y = 0, ny - 1 do begin
        r[*, y] = sqrt(xpos^2 + ypos[y]^2)
        theta[*, y] = atan(ypos[y], xpos)
        z[*, y] = beselj(sqrt(xpos^2 + ypos[y]^2)/a*bz, m)
    endfor
    
    z[where(r gt a)] = 0.
    
    dataset = findgen(nt, nx, ny)
    zrot = fltarr(nx, ny)
    
    for t = 0., nt - 1 do begin
        zrot[where(r le r0)] = cos(n * (theta[where(r lt r0)] + 0.666667*!pi*t/nt))
        zrot[where(r gt r0 and r le r1)] = cos(n * (theta[where(r gt r0 and r le r1)] - 1.33333*!pi*t/nt))
        zrot[where(r gt r1)] = cos(n * (theta[where(r gt r1)] + 2.*!pi*t/nt))
        dataset[t, *, *] = z*zrot
    endfor
    
    t = findgen(nt)/nt*2*!pi
    
    (scope_varfetch('dataset', level=-1, /enter)) = dataset
    (scope_varfetch('t', level=-1, /enter)) = t
    (scope_varfetch('x', level=-1, /enter)) = xpos
    (scope_varfetch('y', level=-1, /enter)) = ypos

end
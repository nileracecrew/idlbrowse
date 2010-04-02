pro vvv, u_in, v_in, xpos_in, ypos_in, length = length, $
    veccolors=veccolors, colorbar=colorbar, bar_charsize=bar_charsize, bar_pad=bar_pad, bar_title=bar_title, lmin=lmin, lmax=lmax, no_reduce=no_reduce, overplot=overplot, max_dim=max_dim, no_arrow_resize=no_arrow_resize, clip=clip, noclip=noclip, _EXTRA=extra

    compile_opt idl2, hidden
    on_error, 2

    u = reform(u_in)
    v = reform(v_in)

    if n_elements(xpos_in) eq 0 then $
        xpos = findgen((size(u))[1]) $
    else $
        xpos = xpos_in

    if n_elements(ypos_in) eq 0 then $
        ypos = findgen((size(u))[2]) $
    else $
        ypos = ypos_in

    if !p.charsize eq 0.0 then begin
        charsize = 1.0 
    endif else begin
        charsize = !p.charsize
    endelse

    if n_elements(extra) ne 0 then begin
        if where(tag_names(extra) eq 'CHARSIZE') ne -1 then $
            charsize = extra.charsize
        if where(tag_names(extra) eq 'XRANGE') ne - 1 then $
            xrange = extra.xrange        
        if where(tag_names(extra) eq 'YRANGE') ne - 1 then $
            yrange = extra.yrange
    endif

    if keyword_set(overplot) then colorbar = 0

    if keyword_set(colorbar) then begin
        if n_elements(bar_charsize) eq 0 then bar_charsize = charsize
        if n_elements(bar_pad) eq 0 then bar_pad = bar_charsize
        old_region = !p.region
        !p.region=[0.0, 0.0, 0.87 - 0.05*bar_pad, 1.0]
    endif

    if n_elements(length) le 0 then length = 1.0

    if n_elements(max_dim) eq 0 then max_dim = 30

    nx = n_elements(xpos)
    ny = n_elements(ypos)
    xbin = (max(xpos) - min(xpos)) / (nx - 1.) 
    ybin = (max(ypos) - min(ypos)) / (ny - 1.) 

    ; find number of bins that will be in the plot area in each direction
    if n_elements(xrange) eq 0 then $
        nxbin = float(nx) $
    else $
        nxbin = abs(float(xrange[1] - xrange[0])/xbin)

    if n_elements(yrange) eq 0 then $
        nybin = float(ny) $
    else $
        nybin = abs(float(yrange[1] - yrange[0])/ybin)

    ; if there will be too many arrows, resize the dataset
    if ((nxbin gt max_dim) || (nybin gt max_dim)) && ~keyword_set(no_reduce) then begin
        nx = long(nx * (max_dim/nxbin)) < nx > 2
        ny = long(ny * (max_dim/nybin)) < ny > 2
        u = congrid(temporary(u), nx, ny, /interp, /minus_one)
        v = congrid(temporary(v), nx, ny, /interp, /minus_one)
        xpos = congrid(xpos, nx, /interp, /minus_one)
        ypos = congrid(ypos, ny, /interp, /minus_one)
        xbin = (max(xpos) - min(xpos)) / (nx - 1.) 
        ybin = (max(ypos) - min(ypos)) / (ny - 1.) 
        if n_elements(veccolors) gt 1 then begin
            veccolors=congrid(temporary(veccolors), nx, ny, /interp, /minus_one)
        endif
    endif 

    mag = sqrt(u^2 + v^2)
    
    if n_elements(veccolors) eq 0 then begin
        veccolors = bytscl(mag, min=lmin, max=lmax, top=253)+1
    endif else if n_elements(veccolors) eq 1 then begin
        veccolors = replicate(veccolors, n_elements(u))
    endif

    if ~keyword_set(overplot) then begin
        if n_elements(xrange) ne 0 then begin
            x0 = xrange[0]
            x1 = xrange[1]
        endif else begin
            x0 = min(xpos) - xbin
            x1 = max(xpos) + xbin
        endelse
        if n_elements(yrange) ne 0 then begin
            y0 = yrange[0]
            y1 = yrange[1]
        endif else begin
            y0 = min(ypos) - ybin
            y1 = max(ypos) + ybin
        endelse

        plot, [x0, x1], [y0, y1], /nodata, /xstyle, /ystyle, _EXTRA = extra
        if keyword_set(colorbar) then begin
            mag = sqrt(u^2. + v^2.)
            if n_elements(lmin) eq 0 then lmin = min(mag)
            if n_elements(lmax) eq 0 then lmax = max(mag)
            if n_elements(extra) ne 0 then begin
                if where(tag_names(extra) eq 'THICK') ne -1 then $
                    thick=extra.thick
                if where(tag_names(extra) eq 'XTHICK') ne -1 then $
                    xthick = extra.xthick
                if where(tag_names(extra) eq 'YTHICK') ne -1 then $
                    ythick = extra.ythick
                if where(tag_names(extra) eq 'CHARTHICK') ne -1 then $
                    chartick = extra.charthick
            endif
            browse_colorbar, /vertical, range=[lmin, lmax], format='(g0.2)', ncolors=254, charsize=bar_charsize, thick=thick, xthick=xthick, ythick=ythick, charthick=charthick, title=bar_title
            !p.region = old_region
        endif
    endif 

    ; scale the radii of the arrows so that they all fit into a xy-bin
    max_r = min([xbin, ybin])/2 * length
    if keyword_set(no_arrow_resize) then $
        r = replicate(max_r, n_elements(u)) $
    else begin
        r = mag/max(mag) > 0.15
        r *= max_r
    endelse
    
    theta = atan(v, u)

    for i = 0, n_elements(u) - 1 do begin
        x = xpos[i mod nx] 
        y = ypos[i / nx] 
        xp0 = r[i] * cos(theta[i])
        yp0 = r[i] * sin(theta[i])
        xp1 = r[i] * cos(theta[i] + 0.85*!pi)
        yp1 = r[i] * sin(theta[i] + 0.85*!pi)
        xp2 = xp1
        yp2 = yp1
        xp3 = r[i] * cos(theta[i] - 0.85*!pi)
        yp3 = r[i] * sin(theta[i] - 0.85*!pi)

        oldexcept = !except
        !except = 0
        polyfill, x+[xp0, xp1, xp2, xp3], y+[yp0, yp1, yp2, yp3], /data, COLOR = veccolors[i], clip=clip, noclip=0
        void = check_math()
        !except = oldexcept
    endfor
end

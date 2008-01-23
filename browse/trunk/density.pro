pro density, z, x, y, zmin=zmin, zmax=zmax, zlog=zlog, colorbar=colorbar, $
             nlevels=nlevels, bar_charsize=bar_charsize, $
             bar_pad=bar_pad, bar_title=bar_title, scientific=scientific, color=color, $
             _extra=extra
    compile_opt idl2, hidden
    on_error, 2

    z_error = 0
    z0 = reform(z)
    if n_elements(nlevels) eq 0 then $
        nlevels=60
    if n_elements(zmax) eq 0 then $
        zmax = max(z0)
    if n_elements(zmin) eq 0 then $
        zmin = min(z0)
    
    if !p.charsize eq 0.0 then begin
        charsize = 1.0 
    endif else begin
        charsize = !p.charsize
    endelse
    if n_elements(extra) ne 0 then begin
        if where(tag_names(extra) eq 'CHARSIZE') ne -1 then begin
            charsize = extra.charsize
        endif
    endif 

    if keyword_set(colorbar) then begin
        if n_elements(bar_charsize) eq 0 then bar_charsize = charsize
        if n_elements(bar_pad) eq 0 then bar_pad = bar_charsize
        old_region = !p.region
        !p.region=[0.0, 0.0, 0.87 - 0.05*bar_pad, 1.0]
    endif


    if keyword_set(zlog) then begin
        if zmin le 0 then begin
            if zmax le 0 then begin
                z_error = 1
                zmin = 1.
                zmax = 9.
            endif else begin
                zmin = min(z0[where(z0 gt 0)])
            endelse
        endif
        zmax = floor(alog10(zmax))+1
        zmin = floor(alog10(zmin))
        spacing = float(zmax - zmin) / nlevels
        levels=10.^(findgen(nlevels)*spacing + zmin)
    endif else begin       
        spacing = float(zmax - zmin) / nlevels
        levels = findgen(nlevels)*spacing + zmin
    endelse
    c_colors = bytscl(levels, top=254)

    if zmin ge zmax then begin
        z_error = 1
        zmax = zmin + 1
    endif

    ; if we have an error, then make a blank plot by undefining the contour levels and zeroing out z0.
    if z_error then begin
        void = temporary(levels)  
        void = temporary(c_colors)
        z0[*,*] = 0
    endif  

    if n_elements(x) eq 0 then $
        contour, z0, /fill, nlevels=nlevels, levels=levels, c_colors=c_colors, _extra=extra $
    else $
        contour, z0, x, y, /fill, nlevels=nlevels, levels=levels, c_colors=c_colors, _extra=extra

    if keyword_set(colorbar) then begin
        if keyword_set(zlog) then $
            format='(%"' + TeXtoIDL('10^{%i}') + '")' $
        else if keyword_set(scientific) then $
            format='(%"%0.2e")' $
        else $
            format='(g0.2)'
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
        colorbar, /vertical, range=[zmin,zmax], format=format, ncolors=254, charsize=bar_charsize, thick=thick, xthick=xthick, ythick=ythick, charthick=charthick, title=bar_title
        !p.region = old_region
    endif
end

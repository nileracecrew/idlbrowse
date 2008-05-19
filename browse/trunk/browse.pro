;------------------------------------------------------------------------------

pro open_z_buffer, x, y, scale
    set_plot, 'z'
    device, set_resolution=[x,y]*scale, decomposed=0, set_character_size=[6,10]
end

;------------------------------------------------------------------------------

pro close_z_buffer
    image=tvrd()
    set_plot, 'X'
    tv, image
end

;------------------------------------------------------------------------------

pro contour_plot, state, scale, no_z_buffer = no_z_buffer
    compile_opt idl2, hidden

    if ~keyword_set(no_z_buffer) then begin
        geometry = widget_info(state.wContour, /geometry)
        open_z_buffer, geometry.draw_xsize, geometry.draw_ysize, scale
    endif

    density_params = { $
        colorbar: 1, xstyle: 1, ystyle: 1, charsize: 1.4*scale, $
        charthick: 1.0*scale, bar_pad: 1.4, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.contour_axes[0].min, state.contour_axes[0].max], $
        yrange: [state.contour_axes[1].min, state.contour_axes[1].max], $
        zmin: state.contour_axes[2].min, zmax: state.contour_axes[2].max, $
        xlog: state.contour_axes[0].log, ylog: state.contour_axes[1].log, $
        zlog: state.contour_axes[2].log, isotropic: state.iso, $
        title: state.labels[0]+ ' = ' + string(state.t[state.ti], $
            format='(g0.5)'), $
        xtitle: state.labels[1], ytitle: state.labels[2], $
        bar_title: state.labels[3] $
    }
    if state.reduced then begin
        if state.zi eq 0 then $
            density, (*state.mag_red)[state.ti, *, *], *state.x_red, $
                *state.y_red, _extra = density_params $
        else $
            density, (*state.data_red)[state.ti, *, *, state.zi - 1], $
                *state.x_red, *state.y_red, _extra = density_params
    endif else begin
        if state.zi eq 0 then $
            density, (*state.mag)[state.ti, *, *], state.x, state.y, $
                _extra = density_params $
        else $
            density, (*state.data)[state.ti, *, *, state.zi - 1], state.x, $
                state.y, _extra = density_params
    endelse

    clip = [ state.contour_axes[0].min, state.contour_axes[1].min, $
        state.contour_axes[0].max, state.contour_axes[1].max ]
    plots, [state.x[state.xi], state.x[state.xi]], $
        [min(state.y), max(state.y)], /data, thick=2.0*scale, clip=clip, $
        noclip=0
    plots, [min(state.x), max(state.x)], $
        [state.y[state.yi], state.y[state.yi]], /data, thick=2.0*scale, $
        clip=clip, noclip=0
    if ~keyword_set(no_z_buffer) then $
        close_z_buffer
end

;------------------------------------------------------------------------------

pro surface_plot, state, scale, shade=shade, no_z_buffer=no_z_buffer
    compile_opt idl2, hidden

    if ~keyword_set(no_z_buffer) then begin
        geometry = widget_info(state.wContour, /geometry)
        open_z_buffer, geometry.draw_xsize, geometry.draw_ysize, scale
    endif

    surface_params = { $
        xstyle: 1, ystyle: 1, zstyle: 1, save: 1, charsize: 2.5*scale, $
        charthick: 1.0*scale, thick: 1.0*scale, t3d: 1, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        xlog: state.contour_axes[0].log, ylog: state.contour_axes[1].log, $
        zlog: state.contour_axes[2].log, isotropic: state.iso, $
        zrange: [state.contour_axes[2].min, state.contour_axes[2].max], $
        title: state.labels[0]+ ' = ' + string(state.t[state.ti], $
            format='(g0.5)'), $
        xtitle: state.labels[1], ytitle: state.labels[2] $
    }

    if state.reduced then begin
        if state.zi eq 0 then begin
            data = reform((*state.mag_red)[state.ti, *, *])
        endif else begin
            data = reform((*state.data_red)[state.ti, *, *, state.zi - 1])
        endelse
        shades = reform(bytscl(data, top=254, min=state.contour_axes[2].min, $
            max=state.contour_axes[2].max))
        xi_red = index(*state.x_red, state.x[state.xi])
        yi_red = index(*state.y_red, state.y[state.yi])
        if keyword_set(shade) then $
            shade_surf, data, *state.x_red, *state.y_red, shades=shades, $
                _extra=surface_params $
        else $
            surface, data, *state.x_red, *state.y_red, shades=shades, $
                _extra=surface_params
        plots, replicate(state.x[state.xi], n_elements(*state.y_red)), $
            *state.y_red, data[xi_red, *], /data, thick=4.*scale, color=255, $
            /t3d
        plots, *state.x_red, replicate(state.y[state.yi], $
            n_elements(*state.x_red)), data[*, yi_red], /data, $
            thick=4.*scale, color=255, /t3d
    endif else begin
        if state.zi eq 0 then begin
            data = reform((*state.mag)[state.ti, *, *] )
        endif else begin
            data = reform((*state.data)[state.ti, *, *, state.zi - 1])
        endelse
        shades = reform(bytscl(data, top=254, min=state.contour_axes[2].min, $
                                max=state.contour_axes[2].max))
        if keyword_set(shade) then $
            shade_surf, data, state.x, state.y, shades=shades, $
                _extra=surface_params $
        else $
            surface, data, state.x, state.y, shades=shades, $
                _extra=surface_params 
        plots, replicate(state.x[state.xi], n_elements(state.y)), state.y, $
            data[state.xi, *], /data, thick=4.*scale, color=255, /t3d
        plots, state.x, replicate(state.y[state.yi], n_elements(state.x)), $
            data[*, state.yi], /data, thick=4.*scale, color=255, /t3d
    endelse          
    if ~keyword_set(no_z_buffer) then $
        close_z_buffer
end

;------------------------------------------------------------------------------

pro shadesurf_plot, state, scale, no_z_buffer=no_z_buffer
    compile_opt idl2, hidden

    surface_plot, state, scale, /shade, no_z_buffer=no_z_buffer
end

;------------------------------------------------------------------------------

pro velovect_plot, state, scale, no_z_buffer=no_z_buffer
    compile_opt idl2, hidden
    
    if ~keyword_set(no_z_buffer) then begin
        geometry = widget_info(state.wContour, /geometry)
        open_z_buffer, geometry.draw_xsize, geometry.draw_ysize, scale
    endif
 
    velovect_params = { $
        colorbar: 1, xstyle: 1, ystyle: 1, charsize: 1.4*scale, $
        charthick: 1.0*scale, bar_pad: 1.4, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.contour_axes[0].min, state.contour_axes[0].max], $
        yrange: [state.contour_axes[1].min, state.contour_axes[1].max], $
        clip: [ state.contour_axes[0].min, state.contour_axes[1].min, $
                state.contour_axes[0].max, state.contour_axes[1].max ], $
        noclip: 0, $
        lmin: state.contour_axes[2].min, lmax: state.contour_axes[2].max, $
        xlog: state.contour_axes[0].log, ylog: state.contour_axes[1].log, $
        isotropic: state.iso, $
        length: state.max_mag[state.ti]/state.contour_axes[2].max*2., $
        title: state.labels[0]+ ' = ' + string(state.t[state.ti], $
            format='(g0.5)'), $
        xtitle: state.labels[1], ytitle: state.labels[2], $
        bar_title: state.labels[3]  $
    } 
    if state.reduced then begin
        velovectcolor, (*state.data_red)[state.ti, *, *, 0], $
            (*state.data_red)[state.ti, *, *, 1], *state.x_red, *state.y_red, $
            _extra=velovect_params
    endif else begin
        velovectcolor, (*state.data)[state.ti, *, *, 0], $
            (*state.data)[state.ti, *, *, 1], state.x, state.y, $
            _extra=velovect_params
    endelse

    plots, [state.x[state.xi], state.x[state.xi]], $
        [min(state.y), max(state.y)], /data, thick=2.0*scale, clip=clip, $
        noclip=0
    plots, [min(state.x), max(state.x)], $
        [state.y[state.yi], state.y[state.yi]], /data, $ 
        thick=2.0*scale, clip=clip, noclip=0
    if ~keyword_set(no_z_buffer) then $
        close_z_buffer
end

;------------------------------------------------------------------------------

pro xcut_plot, state, scale
    compile_opt idl2, hidden

    plot_params = { $
        xstyle: 1, ystyle: 1, charsize: 1.4*scale, charthick: 1.0*scale, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.cut_axes[0,0].min, state.cut_axes[0,0].max], $
        yrange: [state.cut_axes[0,1].min, state.cut_axes[0,1].max], $
        xlog: state.cut_axes[0,0].log, ylog: state.cut_axes[0,1].log, $
        title: state.labels[2] + ' = ' + string(state.y[state.yi], $
            format='(g0.5)'), $
        xtitle: state.labels[1], ytitle: state.labels[3] $
    }
    if state.zi eq 0 then $
        cut = (*state.mag)[state.ti, *, state.yi] $
    else $
        cut = (*state.data)[state.ti, *, state.yi, state.zi - 1]
    plot, state.x, cut, _extra=plot_params
end

;------------------------------------------------------------------------------

pro ycut_plot, state, scale
    compile_opt idl2, hidden

    plot_params = { $
        xstyle: 1, ystyle: 1, charsize: 1.4*scale, charthick: 1.0*scale, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.cut_axes[1,0].min, state.cut_axes[1,0].max], $
        yrange: [state.cut_axes[1,1].min, state.cut_axes[1,1].max], $
        xlog: state.cut_axes[1,0].log, ylog: state.cut_axes[1,1].log, $
        title: state.labels[1] + ' = ' + string(state.x[state.xi], $
            format='(g0.5)'), $
        xtitle: state.labels[2], ytitle: state.labels[3] $
    }
    if state.zi eq 0 then $
        cut = (*state.mag)[state.ti, state.xi, *] $
    else $
        cut = (*state.data)[state.ti, state.xi, *, state.zi - 1]
    plot, state.y, cut, _extra=plot_params 
end

;------------------------------------------------------------------------------

pro contour_redraw, state
    compile_opt idl2, hidden

    plot_state = save_plot_state()
    widget_control, state.wContour, get_value=win
    wset, win
    case widget_info(state.wContourType, /droplist_select) of
        0: contour_plot, state, 1.0
        1: surface_plot, state, 1.0
        2: shadesurf_plot, state, 1.0
        3: velovect_plot, state, 1.0
    endcase
    restore_plot_state, plot_state
end

;------------------------------------------------------------------------------

pro xcut_redraw, state
    compile_opt idl2, hidden

    plot_state = save_plot_state()
    widget_control, state.wXcut, get_value=win
    wset, win
    xcut_plot, state, 1.0
    restore_plot_state, plot_state
end

;------------------------------------------------------------------------------

pro ycut_redraw, state
    compile_opt idl2, hidden

    plot_state = save_plot_state()
    widget_control, state.wYcut, get_value=win
    wset, win
    ycut_plot, state, 1.0
    restore_plot_state, plot_state
end

;------------------------------------------------------------------------------

pro ContourMouseEvent, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy

    if (state.contour_type eq 1) || (state.contour_type eq 2) then begin
        widget_control, state.wContour, get_value=win
        plot_state = save_plot_state()
        wset, win
        xsize = !d.x_size
        ysize = !d.y_size
        restore_plot_state, plot_state
        widget_control, event.id, get_uvalue=rot_state, /no_copy
        if rot_state.rotating then begin
            x = 2.*event.x/xsize - 1
            y = 2.*event.y/ysize - 1 
            rot_state.xangle += -atan(y - rot_state.y0) * 180. / !pi
            rot_state.zangle += atan(x - rot_state.x0) * 180. / !pi 
            rot_state.zangle = rot_state.zangle mod 360.
            ; Keep the data z-axis in the positive device y-axis so that 
            ; left-right mouse motions don't become completely non-intuitive. 
            rot_state.xangle = (rot_state.xangle < 90.) > (-90.)
            scale3, az=rot_state.zangle, ax=rot_state.xangle
            contour_redraw, state
            rot_state.x0 = x
            rot_state.y0 = y
        endif
    
        if event.press then begin
            rot_state.rotating = 1
            rot_state.x0 = 2.*event.x/xsize - 1
            rot_state.y0 = 2.*event.y/ysize - 1
        endif
    
        if event.release then begin
            print, rot_state.xangle, rot_state.zangle
            rot_state.rotating = 0
        endif 
        widget_control, event.id, set_uvalue=rot_state, /no_copy
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro Tslider, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.ti = slider_change(state.wSlider[0], state.t, state.ts)
    contour_redraw, state
    xcut_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro TsliderValue, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.ti = slider_value_change(state.wSlider[0], state.t, state.ts)
    contour_redraw, state
    xcut_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro Xslider, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.yi = slider_change(state.wSlider[1], state.y, state.ys)
    contour_redraw, state
    xcut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro XsliderValue, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.yi = slider_value_change(state.wSlider[1], state.y, state.ys)
    contour_redraw, state
    xcut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro Yslider, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.xi = slider_change(state.wSlider[2], state.x, state.xs)
    contour_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro YsliderValue, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.xi = slider_value_change(state.wSlider[2], state.x, state.xs)
    contour_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

function slider_change, slider, q, qs
    compile_opt idl2, hidden

    widget_control, slider.id, get_value=q_index
    qi = qs[q_index]
    value = string(q[qi], format='(g0.5)')
    widget_control, slider.value_id, set_value=value
    return, qi
end

;------------------------------------------------------------------------------

function slider_value_change, slider, q, qs
    compile_opt idl2, hidden

    widget_control, slider.value_id, get_value=q_value
    qi = qs[index(q, float(q_value))]
    widget_control, slider.id, set_value=qi
    ; update the box with the closest value
    value = string(q[qi], format='(g0.5)')
    widget_control, slider.value_id, set_value=value
    return, qi
end

;------------------------------------------------------------------------------

pro change_slider_range, slider, q
    compile_opt idl2, hidden

    widget_control, slider.min_id, set_value=string(min(q), format='(g0.5)')
    widget_control, slider.max_id, set_value=string(max(q), format='(g0.5)')
    widget_control, slider.id, set_slider_max = n_elements(q)-1
end

;------------------------------------------------------------------------------

pro CutLogButtons, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    cut = strcmp(name, 'YCUT', 4)
    dir = strcmp(strmid(name, 4, 4), 'YLOG')
    state.cut_axes[cut, dir].log = event.select

    if cut eq 0 then $
        xcut_redraw, state $
    else $
        ycut_redraw, state

    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro CutRanges, event
    compile_opt idl2, hidden

    if (TAG_NAMES(event, /STRUCT) eq 'WIDGET_KBRD_FOCUS') $
        && (event.enter eq 1) then $
        return

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    widget_control, event.id, get_value=value
    value = float(value)
    widget_control, event.id, set_value=string(value, format='(g0.5)')
    cut = strcmp(name, 'YCUT', 4)
    dir = strcmp(strmid(name, 4, 1), 'Y')
    changed = 0
    case strmid(name, 5, 3) of
        'MIN': if state.cut_axes[cut, dir].min ne value then begin
                   changed = 1 
                   state.cut_axes[cut, dir].min = value
               endif            
        'MAX': if state.cut_axes[cut, dir].max ne value then begin
                   changed = 1 
                   state.cut_axes[cut, dir].max = value
               endif
    endcase
    if changed then begin
        if cut eq 0 then $
            xcut_redraw, state $
        else $
            ycut_redraw, state
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro ContourLogButtons, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    case strmid(name, 7, 4) of
        'XLOG' : state.contour_axes[0].log = event.select
        'YLOG' : state.contour_axes[1].log = event.select
        'ZLOG' : state.contour_axes[2].log = event.select
    endcase
    contour_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro ContourRanges, event
    compile_opt idl2, hidden

    if (TAG_NAMES(event, /STRUCT) eq 'WIDGET_KBRD_FOCUS') $
        && (event.enter eq 1) then $
        return

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    widget_control, event.id, get_value=value
    value = float(value)
    widget_control, event.id, set_value=string(value, format='(g0.5)')
    dir = strmid(name, 7, 1)
    case dir of
        'X': dir = 0
        'Y': dir = 1
        'Z': dir = 2
    endcase
    changed = 0
    case strmid(name, 8, 3) of 
        'MIN': if state.contour_axes[dir].min ne value then begin
                   changed = 1 
                   state.contour_axes[dir].min = value
               endif            
        'MAX': if state.contour_axes[dir].max ne value then begin
                   changed = 1 
                   state.contour_axes[dir].max = value
               endif
    endcase
    if changed then begin
        max_dim = 200
        nx = (size(*state.data, /dimensions))[1]
        ny = (size(*state.data, /dimensions))[2]
        if (nx gt max_dim) || (ny gt max_dim) then begin
            nt = (size(*state.data, /dimensions))[0]
            xirange = index(state.x, [state.contour_axes[0].min, $
                state.contour_axes[0].max])
            yirange = index(state.y, [state.contour_axes[1].min, $
                state.contour_axes[1].max])
            xirange = xirange[sort(xirange)]
            yirange = yirange[sort(yirange)]
            nx = (xirange[1] - xirange[0] + 1) > 2 
            ny = (yirange[1] - yirange[0] + 1) > 2 
            state.reduced = 1
            ptr_free, state.x_red
            ptr_free, state.y_red
            state.x_red = ptr_new(congrid(state.x[xirange[0]:xirange[1]], $
                nx < max_dim), /no_copy)
            state.y_red = ptr_new(congrid(state.y[yirange[0]:yirange[1]], $
                ny < max_dim), /no_copy)
            if (size(*state.data))[0] eq 4 then begin
                ndir = (size(*state.data, /dimensions))[3]
                data_red = fltarr(nt, nx < max_dim, ny < max_dim, ndir, /nozero)
                for i = 0, ndir - 1 do begin
                    data_red[*, *, *, i] = congrid((*state.data)[ $
                        *, xirange[0]:xirange[1], yirange[0]:yirange[1], i], $
                        nt, nx < max_dim, ny < max_dim)
                endfor 
                mag_red = congrid((*state.mag)[*, xirange[0]:xirange[1], $
                    yirange[0]:yirange[1]], nt, nx < max_dim, ny < max_dim)
            endif else begin
                data_red = congrid((*state.data)[*, xirange[0]:xirange[1], $
                    yirange[0]:yirange[1]], nt, nx < max_dim, ny < max_dim)
                mag_red = 0
            endelse
            ptr_free, state.data_red
            ptr_free, state.mag_red
            state.data_red = ptr_new(data_red, /no_copy)    
            state.mag_red = ptr_new(mag_red, /no_copy)    
        endif else begin
            state.reduced = 0
        endelse
        contour_redraw, state
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro dtValue, event
    compile_opt idl2, hidden
    widget_control, event.top, GET_UVALUE=state, /no_copy
    widget_control, event.id, get_value=value
    state.dt = fix(value) > 1
    widget_control, event.id, set_value=string(state.dt, format='(i0)')
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro PlayButton, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    if state.animating then begin
        state.animating = 0
        widget_control, event.id, set_value = 'Play'
    endif else begin
        state.animating = 1
        widget_control, event.id, set_value = 'Stop'
        widget_control, event.top, timer = 0.1
    endelse
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro ISOButton, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.iso = event.select
    contour_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end
 
;------------------------------------------------------------------------------

pro DimList, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    changed = where(state.wDims eq event.id)
    q = [0, 1, 2]
    for i = 0, 2 do begin
        p = widget_info(state.wDims[i], /droplist_select)
        q[p] = -1
        if (p eq event.index) && (state.wDims[i] ne event.id) then $
            duplicate = state.wDims[i]
    endfor
    if n_elements(duplicate) ne 0 then $
        widget_control, duplicate, set_droplist_select=where(q ne -1)
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro DimList4D, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    if event.index ne state.zi then begin
        state.zi = event.index
        if state.zi eq 0 then begin
            state.contour_axes[2].min = 0.
            state.contour_axes[2].max = state.ref_mag
            state.cut_axes[*, 1].min = 0.
            state.cut_axes[*, 1].max = state.ref_mag
        endif else begin
            if state.contour_type ne 3 then begin
                state.contour_axes[2].min = state.zmin
                state.contour_axes[2].max = state.zmax
            endif
            state.cut_axes[*, 1].min = state.zmin
            state.cut_axes[*, 1].max = state.zmax
        endelse
        widget = widget_info(state.wContourRanges, find_by_uname='CONTOURZMIN')
        widget_control, widget, set_value=string(state.contour_axes[2].min, $
            format='(g0.3)')
        widget = widget_info(state.wContourRanges, find_by_uname='CONTOURZMAX')
        widget_control, widget, set_value=string(state.contour_axes[2].max, $
            format='(g0.3)')
        uname = [ 'XCUT', 'YCUT' ]
        for cut = 0, n_elements(uname) - 1 do begin
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+'YMIN')
            widget_control, widget, set_value=string( $
                state.cut_axes[cut, 1].min, format='(g0.2)') 
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+'YMAX')
            widget_control, widget, set_value=string($ 
                state.cut_axes[cut, 1].max, format='(g0.2)') 
        endfor
        contour_redraw, state
        xcut_redraw, state
        ycut_redraw, state
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro DimSetButton, event
    compile_opt idl2, hidden
    
    widget_control, /hourglass
    widget_control, event.top, GET_UVALUE=state, /no_copy
    p = intarr(3)
    for i = 0, 2 do begin
        p[i] = widget_info(state.wDims[i], /droplist_select)
    endfor

    if array_equal(uniq(p), [0,1,2]) then begin
        newstate = init_state(*state.data, state.t, state.x, state.y, $
            state.labels, p=p, old_p=state.p) 
        newstate.wContour = state.wContour
        newstate.wContourType = state.wContourType
        newstate.wContourToggles = state.wContourToggles
        newstate.wContourRanges = state.wContourRanges
        newstate.wDims = state.wDims
        newstate.wXcut = state.wXCut
        newstate.wYcut = state.wYCut
        newstate.wCutToggles = state.wCutToggles
        newstate.wCutRanges = state.wCutRanges
        newstate.wSlider = state.wSlider

        reset_view, newstate
        widget_control, event.top, SET_UVALUE=newstate, /no_copy
    endif
    widget_control, hourglass=0

end

;------------------------------------------------------------------------------

pro ContourTypeList, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy

    if state.zi ne 0 then begin
        if (event.index eq 3) && (state.contour_type ne 3) then begin
            state.contour_axes[2].min = 0 
            state.contour_axes[2].max = state.ref_mag
            state.cut_axes[*, 1].min = 0.
            state.cut_axes[*, 1].max = state.ref_mag
        endif else if (event.index ne 3) && (state.contour_type eq 3) then begin
            state.contour_axes[2].min = state.zmin 
            state.contour_axes[2].max = state.zmax
            state.cut_axes[*, 1].min = state.zmin
            state.cut_axes[*, 1].max = state.zmax
        endif
    endif
    state.contour_type = event.index

    widget = widget_info(state.wContourRanges, find_by_uname='CONTOURZMIN')
    widget_control, widget, set_value=string(state.contour_axes[2].min, $
        format='(g0.3)')
    widget = widget_info(state.wContourRanges, find_by_uname='CONTOURZMAX')
    widget_control, widget, set_value=string(state.contour_axes[2].max, $
        format='(g0.3)')
    uname = [ 'XCUT', 'YCUT' ]
    for cut = 0, n_elements(uname) - 1 do begin
        widget = widget_info(state.wCutRanges[cut], $
            find_by_uname=uname[cut]+'YMIN')
        widget_control, widget, set_value=string( $
                                state.cut_axes[cut, 1].min, format='(g0.2)') 
        widget = widget_info(state.wCutRanges[cut], $
            find_by_uname=uname[cut]+'YMAX')
        widget_control, widget, set_value=string( $
                                state.cut_axes[cut, 1].max, format='(g0.2)') 
    endfor
    contour_redraw, state
    xcut_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro VAR_event, event
    compile_opt idl2, hidden

    if event.quit then begin
        if event.tag eq 'OK' then begin
            widget_control, event.id, get_value=form
            if form.varname eq '' then begin
                void = dialog_message('Please enter a variable name.', /error, $
                    dialog_parent=event.top)
                return
            endif
            widget_control, event.id, get_uvalue=wBase 
            widget_control, wBase, GET_UVALUE=state, /no_copy
            include_cuts = form.options[0]
            subscripts = form.options[1]

            if subscripts then begin
                if (size(*state.data))[0] eq 4 then begin
                    s = [ state.ti, state.xi, state.yi ]
                    if state.zi ne 0 then begin
                        s = [ s, state.zi - 1 ]
                    endif
                endif else begin
                    if include_cuts then begin
                        s = ([state.ti, state.xi, state.yi ])[state.p]
                    endif else begin
                        s = state.ti
                    endelse
                endelse
                (scope_varfetch(form.varname, /enter, level=1)) = s
            endif else begin
                if state.contour_type eq 3 then begin
                     slice = reform((*state.data)[state.ti, *, *, *])
                endif else begin
                    if state.zi eq 0 then begin
                        slice = reform((*state.mag)[state.ti, *, *])
                    endif else begin
                        slice = reform((*state.data)[state.ti, *, *, $
                                                        state.zi-1])
                    endelse
                endelse
                (scope_varfetch(form.varname, /enter, level=1)) = slice
                if include_cuts then begin
                    (scope_varfetch(form.varname + '_xcut', /enter, level=1)) $
                        = reform(slice[*, state.yi])
                    (scope_varfetch(form.varname + '_ycut', /enter, level=1)) $ 
                        = reform(slice[state.xi, *])  
                endif
            endelse
            widget_control, wBase, set_UVALUE=state, /no_copy
        endif
        widget_control, event.top, /destroy
    endif
end

;------------------------------------------------------------------------------

pro VARButton, event
    compile_opt idl2, hidden

    wVarBase = widget_base(group_leader = event.top, /modal)
    desc = [ $  
     '0, TEXT, , LABEL_LEFT=Enter variable name:, WIDTH=12, TAG=varname', $ 
     '1, BASE,, ROW', $ 
     '2, BUTTON, Include Cuts|Subscripts Only, NO_RELEASE, TAG=options', $ 
     '1, BASE,, ROW', $  
     '0, BUTTON, OK, QUIT, TAG=OK', $  
     '2, BUTTON, Cancel, QUIT, TAG=Cancel'] 
    wVarForm = cw_form(wVarBase, desc, /column)    
    widget_control, wVarForm, set_uvalue=event.top
    widget_control, wVarBase, /realize 
    xmanager, 'var', wVarBase
end
   
;------------------------------------------------------------------------------

pro PNGButton, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    filename = dialog_pickfile(filter='*.png', get_path=path, $
        dialog_parent=event.top)
    if filename ne '' then begin
        scale = 2.0
        geometry = widget_info(state.wContour, /geometry)
        png_xsize = geometry.draw_xsize*scale
        png_ysize = geometry.draw_ysize*scale
        geometry = widget_info(state.wXCut, /geometry)
        png_xcut_xsize = geometry.draw_xsize*scale
        png_xcut_ysize = geometry.draw_ysize*scale
        geometry = widget_info(state.wYCut, /geometry)
        png_ycut_xsize = geometry.draw_xsize*scale
        png_ycut_ysize = geometry.draw_ysize*scale
 
        basename = file_basename(filename, '.png', /fold_case) 
        window, /free, xsize=png_xsize, ysize=png_ysize, /pixmap
        case widget_info(state.wContourType, /droplist_select) of
            0: contour_plot, state, scale, /no_z_buffer
            1: surface_plot, state, scale, /no_z_buffer
            2: shadesurf_plot, state, scale, /no_z_buffer
            3: velovect_plot, state, scale, /no_z_buffer
        endcase
        write_png, path + basename + '.png', tvrd(/true)
        wdelete
        window, /free, xsize=png_xcut_xsize, ysize=png_xcut_ysize, /pixmap
        xcut_plot, state, scale
        write_png, path + basename + '_xcut.png', tvrd(/true)
        wdelete
        window, /free, xsize=png_ycut_xsize, ysize=png_ycut_ysize, /pixmap
        ycut_plot, state, scale
        write_png, path + basename + '_ycut.png', tvrd(/true)
    endif

    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro CTBUpdate, data=event
    compile_opt idl2, hidden
    widget_control, event.top, get_UVALUE=state, /no_copy
    contour_redraw, state
    xcut_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------

pro CTBButton, event
    compile_opt idl2, hidden

    xloadct, group=event.top, /use_current, updatecallback='CTBUpdate', $
        updatecbdata=event
end

;------------------------------------------------------------------------------

pro ResetButton, event
    compile_opt idl2, hidden

    widget_control, event.top, update=0
    widget_control, event.top, GET_UVALUE=state, /no_copy
    reset_view, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
    widget_control, event.top, update=1

end

;------------------------------------------------------------------------------

pro ExitButton, event
    compile_opt idl2, hidden

    widget_control, event.top, /destroy
end

;------------------------------------------------------------------------------

pro Browse_event, event
    compile_opt idl2, hidden

    case tag_names(event, /structure_name) of
        'WIDGET_TIMER': begin
            widget_control, event.top, GET_UVALUE=state, /no_copy
            if state.animating then begin
                widget_control, state.wSlider[0].id, get_value=t_index
                nt = n_elements(state.t)
                t_index += state.dt
                if t_index ge nt then $
                    t_index = t_index mod nt
                widget_control, state.wSlider[0].id, set_value=t_index
                widget_control, state.wSlider[0].id, send_event={id: 0L, $
                    top: 0L, handler: 0L, value: t_index, drag: 0 }
                widget_control, event.top, SET_UVALUE=state, /no_copy
                widget_control, event.top, timer = 0.1
            endif else begin
                widget_control, event.top, SET_UVALUE=state, /no_copy
            endelse
        end
        'WIDGET_BASE': begin
            widget_control, event.top, GET_UVALUE=state, /no_copy
            ; Apparently this is a bad idea in Windows.  See IDL help.
            widget_control, event.top, update=0
            sizes = browse_widget_sizes(event.x > 999, event.y > 600)
            widget_control, state.wContour, draw_xsize=sizes.contour_x, $
                draw_ysize=sizes.contour_y
            widget_control, state.wXCut, draw_xsize=sizes.cut_x, $
                draw_ysize=sizes.cut_y
            widget_control, state.wYCut, draw_xsize=sizes.cut_x, $
                draw_ysize=sizes.cut_y
            widget_control, state.wSlider[0].id, xsize=sizes.tslider
            widget_control, state.wSlider[1].id, ysize=sizes.xslider
            widget_control, state.wSlider[2].id, xsize=sizes.yslider
            uname = ['CONTOUR' + [ 'X', 'Y', 'Z' ] + 'MIN', $
                     'CONTOUR' + [ 'X', 'Y', 'Z' ] + 'MAX']
            for i = 0, n_elements(uname) - 1 do begin
                widget = widget_info(state.wContourRanges, find_by_uname=uname[i])
                widget_control, widget, xsize=sizes.contour_range
            endfor
        
            uname = [ 'XCUT' + [ 'XMIN', 'XMAX', 'YMIN', 'YMAX' ], $
                      'YCUT' + [ 'XMIN', 'XMAX', 'YMIN', 'YMAX' ]]
            for i = 0, n_elements(uname) - 1 do begin
                widget = widget_info(state.wCutRanges[i/4], find_by_uname=uname[i])
                widget_control, widget, xsize=sizes.cut_range
            endfor
            widget_control, event.top, update=1, /clear_events
            contour_redraw, state
            xcut_redraw, state
            ycut_redraw, state
            widget_control, event.top, SET_UVALUE=state, /no_copy
        end
        else: print, 'uncaught event'
    endcase
end

;------------------------------------------------------------------------------

pro BrowseCleanup, wBase
    compile_opt idl2, hidden

    widget_control, wBase, get_uvalue=state, /no_copy
    ptr_free, state.data
    ptr_free, state.data_red
    ptr_free, state.mag
    ptr_free, state.mag_red
end

;------------------------------------------------------------------------------

function save_plot_state
    compile_opt idl2, hidden

    return, { p: !p, x: !x, y: !y, z: !z, w: !d.window }
end

;------------------------------------------------------------------------------

pro restore_plot_state, plot_state
    compile_opt idl2, hidden

    empty
    wset, plot_state.w
    !p = plot_state.p
    !x = plot_state.x
    !y = plot_state.y
    !z = plot_state.z
end

;------------------------------------------------------------------------------

function init_state, data, t_in, x_in, y_in, labels, p=p, old_p=old_p
    compile_opt idl2, hidden

    if n_elements(old_p) eq 0 then old_p = [ 0, 1, 2 ]
    if n_elements(p) ne 0 then begin
        p2 = (sort(old_p))[p]
        data = transpose(temporary(data), p2)
        labels = labels[[p2, 3]]
        case p2[0] of 
            0: t = t_in
            1: t = x_in
            2: t = y_in
        endcase
        case p2[1] of 
            0: x = t_in
            1: x = x_in
            2: x = y_in
        endcase
        case p2[2] of 
            0: y = t_in
            1: y = x_in
            2: y = y_in
        endcase
    endif else begin
        p = [ 0, 1, 2 ]
        t = t_in
        x = x_in
        y = y_in
    endelse
    ts = sort(t)
    xs = sort(x)
    ys = sort(y)   
    nt = n_elements(t)
    nx = n_elements(x)
    ny = n_elements(y)
    max_mag = 0
    mag = 0 
    zi = 1
    if (size(data))[0] eq 4 then begin
        message, 'Calculating vector lengths...', /info
        mag = fltarr(nt, nx, ny, /nozero)
        max_mag = fltarr(nt)
        for i = 0, nt - 1 do begin
            mag[i, *, *] = sqrt(data[i, *, *, 0]^2 + data[i, *, *, 1]^2)
            max_mag[i] = max(mag[i, *, *])  
        endfor
        zi = 0
    endif
    ref_mag = max(max_mag)

    max_dim = 200
    if (nx gt max_dim) || (ny gt max_dim) then begin
        reduced = 1
        x_red = congrid(x, nx < max_dim)
        y_red = congrid(y, ny < max_dim)
        if (size(data))[0] eq 4 then begin
            ndir = (size(data, /dimensions))[3]
            data_red = fltarr(nt, nx < max_dim, ny < max_dim, ndir, /nozero)
            for i = 0, ndir - 1 do begin
                data_red[*, *, *, i] = congrid(data[*, *, *, i], nt, $
                    nx < max_dim, ny < max_dim)
            endfor 
            mag_red = congrid(mag, nt, nx < max_dim, ny < max_dim)
        endif else begin
            data_red = congrid(data, nt, nx < max_dim, ny < max_dim)
            mag_red = 0
            zi = 1
       endelse    
    endif else begin
        reduced = 0
        x_red = 0
        y_red = 0
        data_red = 0
        mag_red = 0
    endelse 

    data = ptr_new(data, /no_copy)
    data_red = ptr_new(data_red, /no_copy)
    mag = ptr_new(mag, /no_copy)
    mag_red = ptr_new(mag_red, /no_copy)
    x_red = ptr_new(x_red, /no_copy)
    y_red = ptr_new(y_red, /no_copy)
    contour_axes = replicate({ log: 0, min: 0.0, max: 0.0 }, 3)
    cut_axes = replicate({ log: 0, min: 0.0, max: 0.0 }, 2, 2)
    wSlider = replicate({ id: 0L, value_id: 0L, min_id: 0L, max_id: 0L }, 3)
    dt = nt / 120 > 1
    state = {  $
        wContour: 0L, $
        wContourType: 0L, $
        wContourToggles: 0L, $
        wContourRanges: 0L, $
        wDims: lonarr(3), $
        wXcut: 0L, $
        wYcut: 0L, $
        wCutToggles: lonarr(2), $
        wCutRanges: lonarr(2), $
        wSlider: wSlider, $
        labels: labels, $
        contour_type: 0, $
        zmin: min(*data), zmax: max(*data), $
        data: data, $
        mag: mag, $
        mag_red: mag_red, $
        x: x, y: y, t: t, $ 
        xs: xs, ys: ys, ts: ts, $
        xi: 0, yi: 0, ti: 0, zi: zi, $
        p: p, $
        ref_mag: ref_mag, $   
        max_mag: max_mag, $
        reduced: reduced, $
        data_red: data_red, $
        x_red: x_red, y_red: y_red, $
        cut_axes: cut_axes, $
        contour_axes: contour_axes, $
        iso: 0, $
        animating: 0, $
        dt: dt $
    }
    return, state
end

;------------------------------------------------------------------------------

pro reset_view, state, no_redraw=no_redraw
    compile_opt idl2, hidden

    state.contour_axes.log = intarr(3)
    state.contour_axes.min = [ min(state.x), min(state.y), state.zmin ]
    state.contour_axes.max = [ max(state.x), max(state.y), state.zmax ]

    rot_state = { xangle: 30., zangle: 30., rotating: 0, x0: 0.0, y0: 0.0 }
    widget_control, state.wContour, set_uvalue=rot_state

    state.cut_axes.log = intarr(2,2)
    state.cut_axes.min = [ [ min(state.x), min(state.y) ], $
        [ state.zmin, state.zmin ] ]
    state.cut_axes.max = [ [ max(state.x), max(state.y) ], $
        [ state.zmax, state.zmax ] ]
    state.iso = 0
    state.xi = state.xs[0]
    state.yi = state.ys[0]
    state.ti = state.ts[0]
    state.contour_type = widget_info(state.wContourType, /droplist_select)

    ; By default, scale the z-axis according to the vector magnitude if we have
    ; a velovect plot or 'MAG' is selected.
    if (state.contour_type eq 3) || (state.zi eq 0) then begin
        state.contour_axes[2].min = 0.
        state.contour_axes[2].max = state.ref_mag
        state.cut_axes[*, 1].min = 0.
        state.cut_axes[*, 1].max = state.ref_mag
    endif

    uname = ['XCUT', 'YCUT' ]
    xlog_enabled = [ min(state.x), min(state.y) ] ge 0.
    ylog_enabled = state.zmin ge 0.
    for i = 0, n_elements(uname) - 1 do begin
        widget = widget_info(state.wCutToggles[i], $
            find_by_uname=uname[i]+'XLOG')
        widget_control, widget, set_button=0, sensitive=xlog_enabled[i]
        widget = widget_info(state.wCutToggles[i], $
            find_by_uname=uname[i]+'YLOG')
        widget_control, widget, set_button=0, sensitive=ylog_enabled
    endfor

    uname = ['CONTOURXLOG', 'CONTOURYLOG', 'CONTOURZLOG', 'ISO' ]
    log_enabled = [ min(state.x), min(state.y), state.zmin, 1. ] ge 0.
    for i = 0, n_elements(uname) - 1 do begin
        widget = widget_info(state.wContourToggles, find_by_uname=uname[i])
        widget_control, widget, set_button=0, sensitive=log_enabled[i]
    endfor

    uname = 'CONTOUR' + [ 'X', 'Y', 'Z' ] + 'MIN'
    for i = 0, n_elements(uname) - 1 do begin
        widget = widget_info(state.wContourRanges, find_by_uname=uname[i])
        widget_control, widget, set_value=string(state.contour_axes[i].min, $
            format='(g0.3)')
    endfor

    uname = 'CONTOUR' + [ 'X', 'Y', 'Z' ] + 'MAX'
    for i = 0, n_elements(uname) - 1 do begin
        widget = widget_info(state.wContourRanges, find_by_uname=uname[i])
        widget_control, widget, set_value=string(state.contour_axes[i].max, $
            format='(g0.3)')
    endfor

    uname = [ 'XCUT', 'YCUT' ]
    for cut = 0, n_elements(uname) - 1 do begin
        dir = [ 'X', 'Y' ]
        for d = 0, 1 do begin
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+dir[d]+'MIN')
            widget_control, widget, set_value=string( $
                state.cut_axes[cut, d].min, format='(g0.2)') 
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+dir[d]+'MAX')
            widget_control, widget, set_value=string( $
                state.cut_axes[cut, d].max, format='(g0.2)') 
        endfor
    endfor
    widget_control, state.wSlider[0].id, set_value=0
    widget_control, state.wSlider[1].id, set_value=0
    widget_control, state.wSlider[2].id, set_value=0
    state.ti = slider_change(state.wSlider[0], state.t, state.ts)
    state.yi = slider_change(state.wSlider[1], state.y, state.ys)
    state.xi = slider_change(state.wSlider[2], state.x, state.xs)
    change_slider_range, state.wSlider[0], state.t
    change_slider_range, state.wSlider[1], state.y
    change_slider_range, state.wSlider[2], state.x
    scale3
    if ~keyword_set(no_redraw) then begin
        contour_redraw, state
        xcut_redraw, state
        ycut_redraw, state
    endif

end

;------------------------------------------------------------------------------

function browse_widget_sizes, x, y
    if n_elements(x) eq 0 then x = 1139
    if n_elements(y) eq 0 then y = 727

    contour_x = .5267*x
    contour_y = .6878*y
    if x lt 1024 then field_delta = -1 else field_delta = 0
    return, {browse_widget_sizes, $
        contour_x: contour_x, $
        contour_y: contour_y, $
        cut_x: .3687*x, $
        cut_y: .3851*y, $
        tslider: 0.69*contour_x, $
        xslider: 0.8*contour_y, $ 
        yslider: 0.75*contour_x, $
        slider_field: 7 + field_delta, $
        contour_range: 8 + field_delta, $
        cut_range: 7 + field_delta $
    }
end

;------------------------------------------------------------------------------


pro browse, data_in, t_in, x_in, y_in, $
        ttitle=ttitle, xtitle=xtitle, ytitle=ytitle, ztitle=ztitle, $
        p=p, f=f, h5select=h5select, help=help

    compile_opt idl2, hidden
    on_error, 2

    if !VERSION.OS_FAMILY ne 'unix' then $
        message, 'Browse currently only supports IDL on Unix platforms' + $
            '(including OS X).'

    if keyword_set(help) || (n_params() eq 0) then begin
        print
        print, 'Usage: browse, <dataset>, <t>, <x>, <y>, [ keywords ]'
        print, '    where <dataset> is a 3D or 4D array, or an h5variable object;'
        print, '    and <t>, <x>, and <y> are 1D arrays of axis labels.'
        print
        print, 'Keywords: '
        print, '    [txyz]title: A string specifying the plot titles for each dimension.'
        print, '                 Defaults are [''t'', ''x'', ''y'', ''data''].'
        print, '    p: 3-element array specifying how to permute the dataset dimensions.'
        print, '       For example, [2, 1, 0] swaps the first and third dimensions.'
        print, '    /f: Set this flag to title the t axis as f (for FFTs).'
        print, '    /h5select: Use the current selection on an h5variable dataset.'
        print, '    /help: Print this message.'
        print
        return
    endif

    restore_plot
    device, decomposed=0
    loadct, 39, /silent
    t3d, /reset
    scale3
    plot_state = save_plot_state()
    sizes = browse_widget_sizes()
    dataset_name = scope_varname(data_in, level=-1) 
    if dataset_name eq '' then dataset_name = 'data'

    widget_control, /hourglass
    
    if (size(data_in, /type) eq 11) && $
                                (obj_class(data_in) eq 'H5VARIABLE') then begin
        message, 'Restoring H5VAR dataset...', /info
        data = reform(data_in->r(select=h5select))
    endif else begin
        ; arg_present checks if data_in was passed by reference (true) or
        ; value (false)
        if arg_present(data_in) then begin
            data = reform(data_in)
        endif else begin
            data = reform(temporary(data_in))
        endelse
    endelse
    nonfinite = where(~finite(data), count)
    if count gt 0 then begin
        message, 'Data has non-finite elements such as NaN or Inf.'
    endif
    dims = size(data, /dimensions)
    ndims = n_elements(dims)
    if ndims eq 2 then begin
        data = reform(temporary(data), 1, dims[0], dims[1])
        if n_elements(t_in) eq 0 then t_in = findgen(1)
        if n_elements(x_in) eq 0 then x_in = findgen(dims[0])
        if n_elements(y_in) eq 0 then y_in = findgen(dims[1])
    endif else begin
        if n_elements(t_in) eq 0 then t_in = findgen(dims[0])
        if n_elements(x_in) eq 0 then x_in = findgen(dims[1])
        if n_elements(y_in) eq 0 then y_in = findgen(dims[2])
    endelse

    s = size(data, /dimensions)
    if (n_elements(t_in) lt s[0]) || (n_elements(x_in) lt s[1]) || $
                                        (n_elements(y_in) lt s[2]) then begin
        message, "T, X, and Y arrays must match the dimensions of the data" + $         
            "array."
    endif

    if keyword_set(f) then $
        labels = [ 'f', 'x', 'y', dataset_name ] $
    else $
        labels = [ 't', 'x', 'y', dataset_name ]

    if n_elements(ttitle) ne 0 then labels[0] = ttitle
    if n_elements(xtitle) ne 0 then labels[1] = xtitle
    if n_elements(ytitle) ne 0 then labels[2] = ytitle
    if n_elements(ztitle) ne 0 then labels[3] = ztitle
    
    state = init_state(data, t_in[0:s[0]-1], x_in[0:s[1]-1], y_in[0:s[2]-1], $
        labels, p=p)

    nt = n_elements(state.t)
    nx = n_elements(state.x)
    ny = n_elements(state.y)

    if dataset_name eq 'data' then $
        window_title = 'Data Browser' $
    else $
        window_title = 'Data Browser - ' + dataset_name
    wBase = widget_base(column=1, title=window_title,  $
        /tlb_size_events, kill_notify='BrowseCleanup')
    wToolbarBase = widget_base(wBase, /base_align_center, /align_left, /row)
    wPlotsBase = WIDGET_BASE(wBase, column=2, /BASE_ALIGN_center, /align_left)
    wContourCol = widget_base(wPlotsBase, /column)
    wCutCol = widget_base(wPlotsBase, /column)

    ; toolbar
    wAnimateBase = widget_base(wToolbarBase, /base_align_center, /row, frame=1)
    wDTLabel = widget_label(wAnimateBase, value='dt:')
    wDTValue = widget_text(wAnimateBase, xsize=5, /edit, event_pro='dtValue', $
        /kbrd_focus_events, value=string(state.dt, format='(i0)'))
    wPlayButton = widget_button(wAnimateBase, value='Play', $
        event_pro='PlayButton')
    ; dimension selector
    case ndims of
        2: begin
           end 
        3: begin    
                dimList = [ 'T', 'X', 'Y' ] + ': (' + string(s, format='(i0)') $
                    + ')'
                wDimsBase = widget_base(wToolbarBase, /base_align_center, $
                    /row, frame=1)
                if n_elements(p) eq 0 then p = [ 0, 1, 2 ]
                for d = 0, 2 do begin
                    state.wDims[d] = widget_droplist(wDimsBase, value=dimList, $
                        event_pro='DimList', title=string(d, format='(i0)')+':')
                    widget_control, state.wDims[d], set_droplist_select=p[d]
                endfor
                wDimsSetButton = widget_button(wDimsBase, value='Set Dims', $
                    event_pro='DimSetButton')
           end
        4: begin
                wDimsBase = widget_base(wToolbarBase, /base_align_center, $
                    /row, frame=1)
                dimList = [ 'MAG', string(indgen((size(*state.data, $
                    /dimensions))[3]), format='(i0)') ]
                state.wDims[0] = widget_droplist(wDimsBase, value = dimList, $
                    event_pro='DimList4D', title='DATA:')    
           end
    endcase

    wContourButtonsBase = widget_base(wToolbarBase, /base_align_center, /row, $
        frame=1)
    if ndims eq 4 then begin
        state.wContourType = widget_droplist(wContourButtonsBase, $
            value=['Contour', 'Surface', 'ShadeSurface', 'Velovect'], $
            event_pro='ContourTypeList')
        widget_control, state.wContourType, set_droplist_select=3
        state.contour_type = 3
    endif else begin
        state.wContourType = widget_droplist(wContourButtonsBase, $
            value=['Contour', 'Surface', 'ShadeSurface'], $
            event_pro='ContourTypeList')
        state.contour_type = 0
    endelse
    wCTBButton = widget_button(wContourButtonsBase, value='Color Table', $
        event_pro='CTBButton')

    wExportBase = widget_base(wToolbarBase, /base_align_center, /row, frame=1)
    wExportLabel = widget_label(wExportBase, value='Export to:')
    wPNGButton = widget_button(wExportBase, value='PNG', event_pro='PNGButton')
    wVARButton = widget_button(wExportBase, value='Variable', $
        event_pro='VARButton')

    wResetBase = widget_base(wToolbarBase, /row)
    wResetButton = widget_button(wResetBase, value='Reset View', $
        event_pro='ResetButton')
 
   ; T slider
    wTsliderBase = widget_base(wContourCol, /row, /align_left, $
        /base_align_center)
    wTSliderLabel = widget_label(wTSliderBase, value='Slice:')
    state.wSlider[0].value_id = widget_text(wTSliderBase, /edit, xsize=10, $
        event_pro='TsliderValue')
    state.wSlider[0].min_id = widget_label(wTsliderBase, value='12345678')
    if ndims eq 2 then begin
        state.wSlider[0].id = widget_slider(wTsliderBase, drag=1, max=1, $
            sensitive=0, event_pro='tslider', xsize=sizes.tslider, /suppress_value) 
        widget_control, state.wSlider[0].value_id, sensitive=0
    endif else $
        state.wSlider[0].id = widget_slider(wTsliderBase, drag=1, max=nt-1, $
            event_pro='tslider', xsize=sizes.tslider, /suppress_value)
    state.wSlider[0].max_id = widget_label(wTsliderBase, value='12345678')  
 
    ; contour column
    dir = ['X', 'Y', 'Z' ]
    wCutSliderBase = lonarr(2)
    wContourBase = widget_base(wContourCol, /row, /base_align_center)
    rot_state = { xangle: 30., zangle: 30., rotating: 0, x0: 0.0, y0: 0.0 }
    state.wContour = WIDGET_DRAW(wContourBase, XSIZE=sizes.contour_x, $         
        YSIZE=sizes.contour_y, /button_events, /motion_events, $
        event_pro='ContourMouseEvent', uvalue=rot_state)
    wCutSliderBase[0] = widget_base(wContourBase, /column, /align_top, $
        /base_align_center)
    wCutSliderBase[1] = widget_base(wContourCol, /row, /align_left, $
        /base_align_center)
    wContourOptions = widget_base(wContourCol, /row, /align_left, $
        /base_align_center)
    wContourToggles = widget_base(wContourOptions, /nonexclusive, $
        row=2, /base_align_center)          
    for i = 0, 2 do begin
        wLogButton = widget_button(wContourToggles, value=dir[i]+'LOG', $
            uname='CONTOUR'+dir[i]+'LOG', event_pro='ContourLogButtons')
        if i eq 1 then $
            wIsoButton = widget_button(wContourToggles, value='ISO ', $
                uname='ISO', event_pro='ISOButton') 
    endfor         
    wContourRanges = widget_base(wContourOptions, row=1, /base_align_center)
    wCurrentRange = widget_base(wContourRanges, /row, /base_align_center)
    for i = 0, 2 do begin
        wRangeLabel = widget_label(wCurrentRange, value=dir[i]+':')
        wRangeMin = widget_text(wCurrentRange, /edit, $
            xsize=sizes.contour_range, uname='CONTOUR'+dir[i]+'MIN', $
            event_pro='ContourRanges', /kbrd_focus_events)
        wRangeMax = widget_text(wCurrentRange, /edit, $
            xsize=sizes.contour_range, uname='CONTOUR'+dir[i]+'MAX', $
            event_pro='ContourRanges', /kbrd_focus_events)
    endfor
    state.wContourToggles = wContourToggles
    state.wContourRanges = wContourRanges

    ; insert x and y cut sliders
    cut = 0 
    state.wSlider[cut+1].value_id = widget_text(wCutSliderBase[cut], /edit, $
        xsize=sizes.slider_field, event_pro='XSliderValue')
    state.wSlider[cut+1].max_id = widget_label(wCutSliderBase[cut], $
        value='12345678')
    state.wSlider[cut+1].id = widget_slider(wCutSliderBase[cut], /drag, $
        /vertical, event_pro='XSlider', ysize=sizes.xslider, /suppress_value)
    state.wSlider[cut+1].min_id = widget_label(wCutSliderBase[cut], $
        value='12345678')
    
    cut = 1
    state.wSlider[cut+1].min_id = widget_label(wCutSliderBase[cut], $
        value='12345678')
    state.wSlider[cut+1].id = widget_slider(wCutSliderBase[cut], /drag, $
        event_pro='YSlider', xsize=sizes.yslider, /suppress_value)
    state.wSlider[cut+1].max_id = widget_label(wCutSliderBase[cut], $
        value='12345678')
    state.wSlider[cut+1].value_id = widget_text(wCutSliderBase[cut], /edit, $
        xsize=sizes.slider_field, event_pro='YSliderValue')
  
    ; cut column
    wCut = lonarr(2)
    dir = ['X', 'Y']
    naxis = [ nx, ny ]
    wCutToggles = lonarr(2)
    wCutRanges = lonarr(2)
    for cut = 0, 1 do begin
        wCut[cut] = widget_draw(wCutCol, xsize=sizes.cut_x, ysize=sizes.cut_y)
        wCutOptions = widget_base(wCutCol, /row, /base_align_center)
        wCutToggles[cut] = widget_base(wCutOptions, /nonexclusive, /row)
        wCutXlogButton = widget_button(wCutToggles[cut], value='XLOG', $
            uname=dir[cut]+'CUTXLOG', event_pro='CutLogButtons')
        wCutYlogButton = widget_button(wCutToggles[cut], value='YLOG', $
            uname=dir[cut]+'CUTYLOG', event_pro='CutLogButtons')
        wCutRanges[cut] = widget_base(wCutOptions, /row, /base_align_center)
        wXRangeLabel = widget_label(wCutRanges[cut], value='X:')
        wXRangeMin = widget_text(wCutRanges[cut], /edit, $
            xsize=sizes.cut_range, uname=dir[cut]+'CUTXMIN', $
            event_pro='CutRanges', /kbrd_focus_events)
        wXRangeMax = widget_text(wCutRanges[cut], /edit, $
            xsize=sizes.cut_range, uname=dir[cut]+'CUTXMAX', $
            event_pro='CutRanges', /kbrd_focus_events)
        wYRangeLabel = widget_label(wCutRanges[cut], value='Y:')
        wYRangeMin = widget_text(wCutRanges[cut], /edit, $
            xsize=sizes.cut_range, uname=dir[cut]+'CUTYMIN', $
            event_pro='CutRanges', /kbrd_focus_events)
        wYRangeMax = widget_text(wCutRanges[cut], /edit, $
            xsize=sizes.cut_range, uname=dir[cut]+'CUTYMAX', $
            event_pro='CutRanges', /kbrd_focus_events)
    endfor
    state.wXCut = wCut[0]
    state.wYCut = wCut[1]
    state.wCutToggles = wCutToggles
    state.wCutRanges = wCutRanges

    widget_control, wBase, /REALIZE
    widget_control, hourglass=0

    reset_view, state
    widget_control, state.wSlider[0].value_id, /input_focus
    widget_control, wBase, set_uvalue=state, /no_copy

    XMANAGER, 'Browse', wBase, /no_block
    restore_plot_state, plot_state

end

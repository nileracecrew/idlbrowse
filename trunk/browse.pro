; VERSION 0.9

; docformat = 'rst'

;------------------------------------------------------------------------------
;+
; Open a Z-buffer for doing buffered writes to the draw widgets.
;
; :Params:
;    x : in, required, type=integer
;       Horizontal size of Z buffer in pixels
;    y : in, required, type=integer
;       Vertical size of Z buffer in pixels
;
;-
pro open_z_buffer, x, y
    set_plot, 'Z'
    device, set_resolution=[x,y], decomposed=0, set_character_size=[6,10]
end

;------------------------------------------------------------------------------
;+
; Write the contents of the Z-buffer to the 'X' device.
;-
pro close_z_buffer
    image=tvrd()
    case strupcase(!version.os_family) of
      'WINDOWS': set_plot, 'WIN'
      'MAC': set_plot, 'MAC'
      else: set_plot, 'X'
    endcase
    tv, image
end

;------------------------------------------------------------------------------
;+
; Create a filled contour or level plot of the current slice.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;
; :Keywords:
;    no_z_buffer : in, optional, type=boolean
;       Set this keyword to disable buffered writes to the graphics window.
;    fill: in, optional, type=boolean
;       If set, makes level plot instead of a filled contour plot.
;-
pro contour_plot, state, scale, no_z_buffer=no_z_buffer, no_fill=no_fill
    compile_opt idl2, hidden

    if ~keyword_set(no_z_buffer) then begin
        geometry = widget_info(state.wMainPlot, /geometry)
        open_z_buffer, geometry.draw_xsize*scale, geometry.draw_ysize*scale
    endif

    if keyword_set(no_fill) then nlevels = 25 else nlevels = 60

    density_params = { $
        fill: ~keyword_set(no_fill), nlevels: nlevels, $ 
        colorbar: 1, xstyle: 1, ystyle: 1, charsize: 1.4*scale, $
        charthick: 1.0*scale, bar_pad: 1.8, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.mainplot_axes[0].min, state.mainplot_axes[0].max], $
        yrange: [state.mainplot_axes[1].min, state.mainplot_axes[1].max], $
        zmin: state.mainplot_axes[2].min, zmax: state.mainplot_axes[2].max, $
        xlog: state.mainplot_axes[0].log, ylog: state.mainplot_axes[1].log, $
        zlog: state.mainplot_axes[2].log, isotropic: state.iso, $
        title: state.titles[0]+ ' = ' + string(state.t[state.ti], $
            format='(g0.5)'), $
        xtitle: state.titles[1], ytitle: state.titles[2], $
        bar_title: state.titles[3] $
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

    clip = [ state.mainplot_axes[0].min, state.mainplot_axes[1].min, $
        state.mainplot_axes[0].max, state.mainplot_axes[1].max ]
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
;+
; Create a level plot of the current slice.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;
; :Keywords:
;    no_z_buffer : in, optional, type=boolean
;       Set this keyword to disable buffered writes to the graphics window.
;-
pro level_plot, state, scale, no_z_buffer=no_z_buffer
    compile_opt idl2, hidden

    contour_plot, state, scale, no_z_buffer=no_z_buffer, /no_fill
end

;------------------------------------------------------------------------------
;+
; Create a wireframe or shaded surface plot of the current slice.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;
; :Keywords:
;    shade : in, optional, type=boolean
;       Set this keyword to use SHADE_SURF instead of SURFACE.
;    no_z_buffer : in, optional, type=boolean
;       Set this keyword to disable buffered writes to the graphics window.
;
;-
pro surface_plot, state, scale, shade=shade, no_z_buffer=no_z_buffer
    compile_opt idl2, hidden

    if ~keyword_set(no_z_buffer) then begin
        geometry = widget_info(state.wMainPlot, /geometry)
        open_z_buffer, geometry.draw_xsize*scale, geometry.draw_ysize*scale
    endif

    surface_params = { $
        xstyle: 1, ystyle: 1, zstyle: 1, save: 1, charsize: 2.5*scale, $
        charthick: 1.0*scale, thick: 1.0*scale, t3d: 1, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        xlog: state.mainplot_axes[0].log, ylog: state.mainplot_axes[1].log, $
        zlog: state.mainplot_axes[2].log, isotropic: state.iso, $
        zrange: [state.mainplot_axes[2].min, state.mainplot_axes[2].max], $
        title: state.titles[0]+ ' = ' + string(state.t[state.ti], $
            format='(g0.5)'), $
        xtitle: state.titles[1], ytitle: state.titles[2] $
    }

    if state.reduced then begin
        if state.zi eq 0 then begin
            data = reform((*state.mag_red)[state.ti, *, *])
        endif else begin
            data = reform((*state.data_red)[state.ti, *, *, state.zi - 1])
        endelse
        shades = reform(bytscl(data, top=254, min=state.mainplot_axes[2].min, $
            max=state.mainplot_axes[2].max))
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
        shades = reform(bytscl(data, top=254, min=state.mainplot_axes[2].min, $
                                max=state.mainplot_axes[2].max))
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
;+
; Create a shaded surface plot of the current slice.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;
; :Keywords:
;    no_z_buffer : in, optional, type=boolean
;       Set this keyword to disable buffered writes to the graphics window.
;
;-
pro shadesurf_plot, state, scale, no_z_buffer=no_z_buffer
    compile_opt idl2, hidden

    surface_plot, state, scale, /shade, no_z_buffer=no_z_buffer
end

;------------------------------------------------------------------------------
;+
; Create a vector field arrow plot of the current slice.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;
; :Keywords:
;    no_z_buffer : in, optional, type=boolean
;       Set this keyword to disable buffered writes to the graphics window.
;
;-
pro velovect_plot, state, scale, no_z_buffer=no_z_buffer
    compile_opt idl2, hidden
    
    if ~keyword_set(no_z_buffer) then begin
        geometry = widget_info(state.wMainPlot, /geometry)
        open_z_buffer, geometry.draw_xsize*scale, geometry.draw_ysize*scale
    endif
 
    velovect_params = { $
        colorbar: 1, xstyle: 1, ystyle: 1, charsize: 1.4*scale, $
        charthick: 1.0*scale, bar_pad: 1.4, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.mainplot_axes[0].min, state.mainplot_axes[0].max], $
        yrange: [state.mainplot_axes[1].min, state.mainplot_axes[1].max], $
        clip: [ state.mainplot_axes[0].min, state.mainplot_axes[1].min, $
                state.mainplot_axes[0].max, state.mainplot_axes[1].max ], $
        noclip: 0, $
        lmin: state.mainplot_axes[2].min, lmax: state.mainplot_axes[2].max, $
        xlog: state.mainplot_axes[0].log, ylog: state.mainplot_axes[1].log, $
        isotropic: state.iso, $
        title: state.titles[0]+ ' = ' + string(state.t[state.ti], $
            format='(g0.5)'), $
        xtitle: state.titles[1], ytitle: state.titles[2], $
        bar_title: state.titles[3],  $
        no_arrow_resize: 1 $
    } 
    if state.reduced then begin
        vvv, (*state.data_red)[state.ti, *, *, 0], $
            (*state.data_red)[state.ti, *, *, 1], *state.x_red, *state.y_red, $
            _extra=velovect_params
    endif else begin
        vvv, (*state.data)[state.ti, *, *, 0], $
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
;+
; Create a filled contour or level plot of the current slice with an 
; overplotted vector field.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;
; :Keywords:
;    no_z_buffer : in, optional, type=boolean
;       Set this keyword to disable buffered writes to the graphics window.
;    fill: in, optional, type=boolean
;       If set, makes level plot instead of a filled contour plot.
;-
pro velovect_contour_plot, state, scale, no_z_buffer=no_z_buffer, no_fill=no_fill
    compile_opt idl2, hidden

    if ~keyword_set(no_z_buffer) then begin
        geometry = widget_info(state.wMainPlot, /geometry)
        open_z_buffer, geometry.draw_xsize*scale, geometry.draw_ysize*scale
    endif

    if keyword_set(no_fill) then nlevels = 25 else nlevels = 60

    density_params = { $
        fill: ~keyword_set(no_fill), nlevels: nlevels, $ 
        colorbar: 1, xstyle: 1, ystyle: 1, charsize: 1.4*scale, $
        charthick: 1.0*scale, bar_pad: 1.8, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.mainplot_axes[0].min, state.mainplot_axes[0].max], $
        yrange: [state.mainplot_axes[1].min, state.mainplot_axes[1].max], $
        zmin: state.mainplot_axes[2].min, zmax: state.mainplot_axes[2].max, $
        xlog: state.mainplot_axes[0].log, ylog: state.mainplot_axes[1].log, $
        zlog: state.mainplot_axes[2].log, isotropic: state.iso, $
        title: state.titles[0]+ ' = ' + string(state.t[state.ti], $
            format='(g0.5)'), $
        xtitle: state.titles[1], ytitle: state.titles[2], $
        bar_title: state.titles[3] $
    }
    if state.reduced then begin
        xbin = (max(*state.x_red) - min(*state.x_red)) / (n_elements(*state.x_red) - 1.) * 2
        ybin = (max(*state.y_red) - min(*state.y_red)) / (n_elements(*state.y_red) - 1.) * 2
        density_params.xrange += [ -xbin, xbin ]
        density_params.yrange += [ -ybin, ybin ]
        if state.zi eq 0 then $
            density, (*state.mag_red)[state.ti, *, *], *state.x_red, $
                *state.y_red, _extra = density_params $
        else $
            density, (*state.data_red)[state.ti, *, *, state.zi - 1], $
                *state.x_red, *state.y_red, _extra = density_params
        vvv, (*state.data_red)[state.ti, *, *, 0], (*state.data_red)[state.ti, *, *, 1], *state.x_red, *state.y_red, length=1.0*state.max_mag[state.ti]/state.ref_mag, veccolor=255, /overplot
    endif else begin
        xbin = (max(state.x) - min(state.x)) / (n_elements(state.x) - 1.) * 2
        ybin = (max(state.y) - min(state.y)) / (n_elements(state.y) - 1.) * 2
        ;density_params.xrange += [ -xbin, xbin ]
       ;xs density_params.yrange += [ -ybin, ybin ]
        if state.zi eq 0 then $
            density, (*state.mag)[state.ti, *, *], state.x, state.y, $
                _extra = density_params $
        else $
            density, (*state.data)[state.ti, *, *, state.zi - 1], state.x, $
                state.y, _extra = density_params
        vvv, (*state.data)[state.ti, *, *, 0], (*state.data)[state.ti, *, *, 1], state.x, state.y, length=1.0*state.max_mag[state.ti]/state.ref_mag, veccolors=255, /overplot
    endelse

    clip = [ state.mainplot_axes[0].min, state.mainplot_axes[1].min, $
        state.mainplot_axes[0].max, state.mainplot_axes[1].max ]
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
;+
; Draws the x-cut plot for the current slice and y selection.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;-
pro xcut_plot, state, scale
    compile_opt idl2, hidden

    plot_params = { $
        xstyle: 1, ystyle: 1, charsize: 1.4*scale, charthick: 1.0*scale, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.cut_axes[0,0].min, state.cut_axes[0,0].max], $
        yrange: [state.cut_axes[0,1].min, state.cut_axes[0,1].max], $
        xlog: state.cut_axes[0,0].log, ylog: state.cut_axes[0,1].log, $
        title: state.titles[2] + ' = ' + string(state.y[state.yi], $
            format='(g0.5)'), $
        xtitle: state.titles[1], ytitle: state.titles[3], xmargin: [12,3] $
    }
    if state.zi eq 0 then $
        cut = (*state.mag)[state.ti, *, state.yi] $
    else $
        cut = (*state.data)[state.ti, *, state.yi, state.zi - 1]
    plot, state.x, cut, _extra=plot_params
end

;------------------------------------------------------------------------------
;+
; Draws the y-cut plot for the current slice and x selection.  
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;    scale : in, required, type=float
;       Size of a plot pixel relative to a device pixel. Usually set to 1.
;-
pro ycut_plot, state, scale
    compile_opt idl2, hidden

    plot_params = { $
        xstyle: 1, ystyle: 1, charsize: 1.4*scale, charthick: 1.0*scale, $
        xthick: 1.0*scale, ythick: 1.0*scale, zthick: 1.0*scale, $
        thick: 1.0*scale, $
        xrange: [state.cut_axes[1,0].min, state.cut_axes[1,0].max], $
        yrange: [state.cut_axes[1,1].min, state.cut_axes[1,1].max], $
        xlog: state.cut_axes[1,0].log, ylog: state.cut_axes[1,1].log, $
        title: state.titles[1] + ' = ' + string(state.x[state.xi], $
            format='(g0.5)'), $
        xtitle: state.titles[2], ytitle: state.titles[3], xmargin: [12,3] $
    }
    if state.zi eq 0 then $
        cut = (*state.mag)[state.ti, state.xi, *] $
    else $
        cut = (*state.data)[state.ti, state.xi, *, state.zi - 1]
    plot, state.y, cut, _extra=plot_params 
end

;------------------------------------------------------------------------------
;+
; Draws the 2D slice plot in the main plot area. 
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;-
pro mainplot_redraw, state, scale=scale
    compile_opt idl2, hidden

    if n_elements(scale) eq 0 then scale = 1.0

    plot_state = current_plot_state()
    scale3, az=state.rotation.zangle, ax=state.rotation.xangle
    widget_control, state.wMainPlot, get_value=win
    wset, win
    case widget_info(state.wMainPlotType, /droplist_select) of
        0: contour_plot, state, 1.0
        1: level_plot, state, 1.0
        2: surface_plot, state, 1.0
        3: shadesurf_plot, state, 1.0
        4: velovect_plot, state, 1.0
        5: velovect_contour_plot, state, 1.0
    endcase
    restore_plot_state, plot_state
end

;------------------------------------------------------------------------------
;+
; Draws the x-cut plot in its plot area. 
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;-
pro xcut_redraw, state
    compile_opt idl2, hidden

    plot_state = current_plot_state()
    widget_control, state.wXcut, get_value=win
    wset, win
    xcut_plot, state, 1.0
    restore_plot_state, plot_state
end

;------------------------------------------------------------------------------
;+
; Draws the y-cut plot in its plot area. 
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;-
pro ycut_redraw, state
    compile_opt idl2, hidden

    plot_state = current_plot_state()
    widget_control, state.wYcut, get_value=win
    wset, win
    ycut_plot, state, 1.0
    restore_plot_state, plot_state
end

;------------------------------------------------------------------------------
;+
; Handle events from the slice selection slider.
;
; :Params:
;    event : in, required, type=struct
;       Slider event structure
;-
pro Tslider, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.ti = slider_index(state.wSlider[0], state.t, state.ts)
    mainplot_redraw, state
    xcut_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the slice selection text field.
;
; :Params:
;    event : in, required, type=struct
;       Field event structure
;-
pro TsliderValue, event
    compile_opt idl2, hidden

    ; don't do anything if entering the text field 
    if (TAG_NAMES(event, /STRUCT) eq 'WIDGET_KBRD_FOCUS') $
        && (event.enter eq 1) then $
        return

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.ti = slider_value_index(state.wSlider[0], state.t, state.ts)
    mainplot_redraw, state
    xcut_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the x-cut selection slider.
;
; :Params:
;    event : in, required, type=struct
;       Slider event structure
;-
pro Xslider, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.yi = slider_index(state.wSlider[1], state.y, state.ys)
    mainplot_redraw, state
    xcut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the x-cut selection text field.
;
; :Params:
;    event : in, required, type=struct
;       Field event structure
;-
pro XsliderValue, event
    compile_opt idl2, hidden

    ; don't do anything if entering the text field 
    if (TAG_NAMES(event, /STRUCT) eq 'WIDGET_KBRD_FOCUS') $
        && (event.enter eq 1) then $
        return

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.yi = slider_value_index(state.wSlider[1], state.y, state.ys)
    mainplot_redraw, state
    xcut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the y-cut selection slider.
;
; :Params:
;    event : in, required, type=struct
;       Slider event structure
;-
pro Yslider, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.xi = slider_index(state.wSlider[2], state.x, state.xs)
    mainplot_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the y-cut selection text field.
;
; :Params:
;    event : in, required, type=struct
;       Field event structure
;-
pro YsliderValue, event
    compile_opt idl2, hidden

    ; don't do anything if entering the text field 
    if (TAG_NAMES(event, /STRUCT) eq 'WIDGET_KBRD_FOCUS') $
        && (event.enter eq 1) then $
        return

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.xi = slider_value_index(state.wSlider[2], state.x, state.xs)
    mainplot_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Returns the index in the axis coordinate array corresponding to the current
; slider position. Also updates the slider's text field with the current 
; coordinate. 
;
; :Returns: integer
;
; :Params:
;    slider: in, required, type=struct
;       structure with widget info for the slider and its text field
;    q: in, required, type=fltarr()
;       array of axis coordinates 
;    qs: in, required, type=intarr()
;       result of sort(q)
;-
function slider_index, slider, q, qs
    compile_opt idl2, hidden

    widget_control, slider.id, get_value=q_index
    qi = qs[q_index]
    value = string(q[qi], format='(g0.5)')
    ; update the text field with the new value
    widget_control, slider.value_id, set_value=value
    return, qi
end

;------------------------------------------------------------------------------
;+
; Returns the index in the axis coordinate array corresponding to the current
; value of the slider text field.  Also moves the slider to match this index.
;
; :Returns: integer
;
; :Params:
;    slider: in, required, type=struct
;       structure with widget info for the slider and its text field
;    q: in, required, type=fltarr()
;       array of axis coordinates 
;    qs: in, required, type=intarr()
;       result of sort(q)
;-
function slider_value_index, slider, q, qs
    compile_opt idl2, hidden

    widget_control, slider.value_id, get_value=q_value
    ; update the slider text field with the value nearest to the entered value
    qi = qs[index(q, float(q_value))]
    widget_control, slider.id, set_value=qi
    value = string(q[qi], format='(g0.5)')
    widget_control, slider.value_id, set_value=value
    return, qi
end

;------------------------------------------------------------------------------
;+
; Change the minimum and maximum labels for a slider.
;
; :Params:
;    slider: in, required, type=struct
;       structure with widget info for the slider and its text field
;    q: in, required, type=fltarr()
;       array of axis coordinates 
;-
pro change_slider_range, slider, q
    compile_opt idl2, hidden

    widget_control, slider.min_id, set_value=string(min(q), format='(g0.5)')
    widget_control, slider.max_id, set_value=string(max(q), format='(g0.5)')
    widget_control, slider.id, set_slider_max = n_elements(q)-1
end

;------------------------------------------------------------------------------
;+
; Change a plot range field.
;
; :Params:
;    id: in, required, type=integer
;       widget id of the text field
;    value: in, required, type=float
;       value to set
;-
pro change_range_field, id, value
    compile_opt idl2, hidden

    widget_control, id, set_value=string(value, format='(g0.3)')
end

;------------------------------------------------------------------------------
;+
; Handle events from the XLOG/YLOG toggles for the cut plots.
;
; :Params:
;    event: in, required, type=struct
;       button event structure
;-
pro CutLogButtons, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    ; strcmp returns 0 for 'XCUT', 1 for 'YCUT'
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
;+
; Handle events from the xrange/yrange fields for the cut plots.
;
; :Params:
;    event: in, required, type=struct
;       text field event structure
;-
pro CutRanges, event
    compile_opt idl2, hidden

    ; don't do anything if entering the text field 
    if (TAG_NAMES(event, /STRUCT) eq 'WIDGET_KBRD_FOCUS') $
        && (event.enter eq 1) then $
        return

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    widget_control, event.id, get_value=value
    value = float(value)
    change_range_field, event.id, value
    ; strcmp returns 0 for 'XCUT', 1 for 'YCUT'
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
;+
; Handle mouse events in the main plot area.  At present, this only handles 
; rotating the surface plots.
;
; :Params:
;    event : in, required, type=struct
;       Mouse event structure
;-
pro MainPlotMouseEvent, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy

    if (state.mainplot_type eq 2) || (state.mainplot_type eq 3) then begin
        widget_control, state.wMainPlot, get_value=win
        plot_state = current_plot_state()
        wset, win
        xsize = !d.x_size
        ysize = !d.y_size
        restore_plot_state, plot_state
        
        if state.rotation.rotating then begin
            x = 2.*event.x/xsize - 1
            y = 2.*event.y/ysize - 1 
            state.rotation.xangle += -atan(y - state.rotation.y0) * 180. / !pi
            state.rotation.zangle += atan(x - state.rotation.x0) * 180. / !pi 
            state.rotation.zangle = state.rotation.zangle mod 360.
            ; Keep the data z-axis in the positive device y-axis so that 
            ; left-right mouse motions don't become completely non-intuitive. 
            state.rotation.xangle = (state.rotation.xangle < 90.) > (-90.)
            mainplot_redraw, state
            state.rotation.x0 = x
            state.rotation.y0 = y
        endif
    
        if event.press then begin
            state.rotation.rotating = 1
            state.rotation.x0 = 2.*event.x/xsize - 1
            state.rotation.y0 = 2.*event.y/ysize - 1
        endif
    
        if event.release then begin
            state.rotation.rotating = 0
        endif 
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the XLOG/YLOG/ZLOG toggles for the main plot.
;
; :Params:
;    event: in, required, type=struct
;       button event structure
;-
pro MainPlotLogButtons, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    case strmid(name, 7, 4) of
        'XLOG' : state.mainplot_axes[0].log = event.select
        'YLOG' : state.mainplot_axes[1].log = event.select
        'ZLOG' : state.mainplot_axes[2].log = event.select
    endcase
    mainplot_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the [xyz]range fields for the main plot.
;
; :Params:
;    event: in, required, type=struct
;       text field event structure
;-
pro MainPlotRanges, event
    compile_opt idl2, hidden

    ; don't do anything if entering the text field 
    if (TAG_NAMES(event, /STRUCT) eq 'WIDGET_KBRD_FOCUS') $
        && (event.enter eq 1) then $
        return

    widget_control, event.top, GET_UVALUE=state, /no_copy
    name = widget_info(event.id, /uname)
    widget_control, event.id, get_value=value
    value = float(value)
    change_range_field, event.id, value
    dir = strmid(name, 7, 1)
    case dir of
        'X': dir = 0
        'Y': dir = 1
        'Z': dir = 2
    endcase
    changed = 0
    case strmid(name, 8, 3) of 
        'MIN': if state.mainplot_axes[dir].min ne value then begin
                   changed = 1 
                   state.mainplot_axes[dir].min = value
               endif            
        'MAX': if state.mainplot_axes[dir].max ne value then begin
                   changed = 1 
                   state.mainplot_axes[dir].max = value
               endif
    endcase
    if changed then begin
        ; if necessary, re-decimate the data based on the new [xy]ranges
        max_dim = 200
        nx = (size(*state.data, /dimensions))[1]
        ny = (size(*state.data, /dimensions))[2]
        if (nx gt max_dim) || (ny gt max_dim) then begin
            nt = (size(*state.data, /dimensions))[0]
            ; find the indices for the min and max coordinate values
            xirange = index(state.x, [state.mainplot_axes[0].min, $
                state.mainplot_axes[0].max])
            yirange = index(state.y, [state.mainplot_axes[1].min, $
                state.mainplot_axes[1].max])
            xirange = xirange[sort(xirange)]
            yirange = yirange[sort(yirange)]
            if yirange[0] eq yirange[1] then $
                yirange = [0, 1]
            ; need at least 2 points else the dimension will be reformed away
            nx = (xirange[1] - xirange[0] + 1) > 2 
            ny = (yirange[1] - yirange[0] + 1) > 2 
            state.reduced = 1
            ptr_free, state.x_red
            ptr_free, state.y_red
            state.x_red = ptr_new(congrid(state.x[xirange[0]:xirange[1]], $
                nx < max_dim, /minus_one, /interp), /no_copy)
            state.y_red = ptr_new(congrid(state.y[yirange[0]:yirange[1]], $
                ny < max_dim, /minus_one, /interp), /no_copy)
            if (size(*state.data))[0] eq 4 then begin
                ndir = (size(*state.data, /dimensions))[3]
                data_red = fltarr(nt, nx < max_dim, ny < max_dim, ndir, $
                    /nozero)
                for i = 0, ndir - 1 do begin
                    data_red[*, *, *, i] = congrid((*state.data)[ $
                        *, xirange[0]:xirange[1], yirange[0]:yirange[1], i], $
                        nt, nx < max_dim, ny < max_dim, /minus_one, /interp)
                endfor 
                mag_red = congrid((*state.mag)[*, xirange[0]:xirange[1], $
                    yirange[0]:yirange[1]], nt, nx < max_dim, ny < max_dim, $
                    /minus_one, /interp)
            endif else begin
                data_red = congrid((*state.data)[*, xirange[0]:xirange[1], $
                    yirange[0]:yirange[1]], nt, nx < max_dim, ny < max_dim, $
                    /minus_one, /interp)
                mag_red = 0
            endelse

            ptr_free, state.data_red
            ptr_free, state.mag_red
            state.data_red = ptr_new(data_red, /no_copy)    
            state.mag_red = ptr_new(mag_red, /no_copy)    
        endif else begin
            state.reduced = 0
        endelse
        mainplot_redraw, state
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the animation "dt" text field.
;
; :Params:
;    event: in, required, type=struct
;       text field event structure
;-
pro dtValue, event
    compile_opt idl2, hidden
    widget_control, event.top, GET_UVALUE=state, /no_copy
    widget_control, event.id, get_value=value
    state.dt = fix(value) > 1
    widget_control, event.id, set_value=string(state.dt, format='(i0)')
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Start and stop animation.
;
; :Params:
;    event: in, required, type=struct
;       button event structure
;-
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
;+
; Set/unset the "isotropic" flag for the main plot window.
;
; :Params:
;    event: in, required, type=struct
;       button event structure
;-
pro ISOButton, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    state.iso = event.select
    mainplot_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end
 
;------------------------------------------------------------------------------
;+
; Handle events from the dimension order droplists.  The dataset is not actually
; changed until DimSetButton is clicked.
;
; :Params:
;    event: in, required, type=struct
;       droplist event structure
;-
pro DimList, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    selected_dims = intarr(3)
    ; loop through each droplist to find the one that has the same value as the
    ; one that was just changed.
    for i = 0, 2 do begin
        p = widget_info(state.wDims[i], /droplist_select)
        selected_dims[p] = 1
        ; if this droplist has the same value as the one that just changed, and
        ; it is not the droplist that just changed, then it is the duplicate.
        if (p eq event.index) && (state.wDims[i] ne event.id) then $
            duplicate = state.wDims[i]
    endfor
    ; switch the duplicate to the selected droplist's old value.
    if n_elements(duplicate) ne 0 then $
        widget_control, duplicate, set_droplist_select=where(selected_dims eq 0)
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the dimension droplist for 4D datasets.  It selects which
; vector direction (x, y, or vector magnitude) is used for contour/surface
; plots and cuts.  The axis ranges for the plots are reset to min and max of
; the selected direction as appropriate.
;
; :Params:
;    event: in, required, type=struct
;       droplist event structure
;-
pro DimList4D, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    if event.index ne state.zi then begin
        ; reset the zranges to something appropriate for the new selection
        if event.index eq 0 then begin
            state.mainplot_axes[2].min = 0.
            state.mainplot_axes[2].max = state.ref_mag
            state.cut_axes[*, 1].min = 0.
            state.cut_axes[*, 1].max = state.ref_mag
        endif else begin
            ; only change zrange if old state.zi was 0
            if state.zi eq 0 then begin
                state.mainplot_axes[2].min = state.zmin
                state.mainplot_axes[2].max = state.zmax
            endif
            state.cut_axes[*, 1].min = state.zmin
            state.cut_axes[*, 1].max = state.zmax
        endelse
        state.zi = event.index
        widget = widget_info(state.wMainPlotRanges, find_by_uname='CONTOURZMIN')
        change_range_field, widget, state.mainplot_axes[2].min
        widget = widget_info(state.wMainPlotRanges, find_by_uname='CONTOURZMAX')
        change_range_field, widget, state.mainplot_axes[2].max

        uname = [ 'XCUT', 'YCUT' ]
        for cut = 0, n_elements(uname) - 1 do begin
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+'YMIN')
            change_range_field, widget, state.cut_axes[cut, 1].min
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+'YMAX')
            change_range_field, widget, state.cut_axes[cut, 1].max
        endfor
        mainplot_redraw, state
        xcut_redraw, state
        ycut_redraw, state
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the 'Set Dims' button.  Reinitializes the state structure 
; with the new dimension permutation.
;
; :Params:
;    event: in, required, type=struct
;       button event structure
;-
pro DimSetButton, event
    compile_opt idl2, hidden
    
    widget_control, /hourglass
    widget_control, event.top, GET_UVALUE=state, /no_copy
    p = intarr(3)
    for i = 0, 2 do begin
        p[i] = widget_info(state.wDims[i], /droplist_select)
    endfor

    ; make sure the selected permutation is valid, and then update the state.
    if array_equal(uniq(p), [0,1,2]) then begin
        newstate = init_state(*state.data, state.t, state.x, state.y, $
            state.titles, p=p, old_p=state.p) 
        free_browse_pointers, state
        newstate.wMainPlot = state.wMainPlot
        newstate.wMainPlotType = state.wMainPlotType
        newstate.wMainPlotToggles = state.wMainPlotToggles
        newstate.wMainPlotRanges = state.wMainPlotRanges
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
;+
; Handle events from the droplist that selects the main plot type (contour,
; surface, etc.)
; of plot for the main plot.
;
; :Params:
;    event: in, required, type=struct
;       droplist event structure
;-
pro MainPlotTypeList, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy

    if state.mainplot_type ne event.index then begin
        changed = 0
        if state.zi ne 0 then begin
            ; change zranges to vector magnitude if going to velovect
            if (event.index eq 4) then begin
                changed = 1
                state.mainplot_axes[2].min = 0 
                state.mainplot_axes[2].max = state.ref_mag
                state.cut_axes[*, 1].min = 0.
                state.cut_axes[*, 1].max = state.ref_mag
            ; change zranges to zmin/zmax if coming from velovect
            endif else if (state.mainplot_type eq 4) then begin
                changed = 1
                state.mainplot_axes[2].min = state.zmin 
                state.mainplot_axes[2].max = state.zmax
                state.cut_axes[*, 1].min = state.zmin
                state.cut_axes[*, 1].max = state.zmax
            endif
        endif
    
        ; update any xyzrange text fields if there were changes
        if changed then begin 
            widget = widget_info(state.wMainPlotRanges, $
                find_by_uname='CONTOURZMIN')
            change_range_field, widget, state.mainplot_axes[2].min
            widget = widget_info(state.wMainPlotRanges, $
                find_by_uname='CONTOURZMAX')
            change_range_field, widget, state.mainplot_axes[2].max
            uname = [ 'XCUT', 'YCUT' ]
            for cut = 0, n_elements(uname) - 1 do begin
                widget = widget_info(state.wCutRanges[cut], $
                    find_by_uname=uname[cut]+'YMIN')
                change_range_field, widget, state.cut_axes[cut, 1].min
                widget = widget_info(state.wCutRanges[cut], $
                    find_by_uname=uname[cut]+'YMAX')
                change_range_field, widget, state.cut_axes[cut, 1].max
            endfor
        endif

        state.mainplot_type = event.index
        mainplot_redraw, state
        xcut_redraw, state
        ycut_redraw, state
    endif
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the 'Export to variable' form.  See VARButton below.
;
; :Params:
;   event: in, required, type=struct
;       form event structure
;-
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
                if state.mainplot_type eq 4 then begin
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
;+
; Handle events from the 'Export to variable' button. The actual export is 
; performed by VAR_event.
;
; :Params:
;   event: in, required, type=struct
;       button event structure
;-
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
;+
; Handle events from the 'Export to PNG' button. 
;
; :Params:
;   event: in, required, type=struct
;       button event structure
;-
pro PNGButton, event
    compile_opt idl2, hidden

    widget_control, event.top, GET_UVALUE=state, /no_copy
    filename = dialog_pickfile(filter='*.png', get_path=path, $
        dialog_parent=event.top)
    if filename ne '' then begin
        scale = 2.0
        geometry = widget_info(state.wMainPlot, /geometry)
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
        scale3, az=state.rotation.zangle, ax=state.rotation.xangle
        case widget_info(state.wMainPlotType, /droplist_select) of
            0: contour_plot, state, scale, /no_z_buffer
            1: level_plot, state, scale, /no_z_buffer
            2: surface_plot, state, scale, /no_z_buffer
            3: shadesurf_plot, state, scale, /no_z_buffer
            4: velovect_plot, state, scale, /no_z_buffer
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
;+
; Redraws the plots after a color table update from xloadct.
;
; :Keywords:
;    data: in, required, type=struct
;       color table event structure
;-
pro CTBUpdate, data=event
    compile_opt idl2, hidden
    widget_control, event.top, get_UVALUE=state, /no_copy
    mainplot_redraw, state
    xcut_redraw, state
    ycut_redraw, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
end

;------------------------------------------------------------------------------
;+
; Handle events from the 'Color Table' button. 
;
; :Params:
;   event: in, required, type=struct
;       button event structure
;-
pro CTBButton, event
    compile_opt idl2, hidden

    xloadct, group=event.top, /use_current, updatecallback='CTBUpdate', $
        updatecbdata=event
end

;------------------------------------------------------------------------------
;+
; Handle events from the 'Reset View' button.
;
; :Params:
;    event: in, required, type=struct
;       button event structure
pro ResetButton, event
    compile_opt idl2, hidden

    widget_control, event.top, update=0
    widget_control, event.top, GET_UVALUE=state, /no_copy
    reset_view, state
    widget_control, event.top, SET_UVALUE=state, /no_copy
    widget_control, event.top, update=1

end

;------------------------------------------------------------------------------
;+
; Handle events from the top level base.  Currently handles the animation timer
; and resize events.
;
; :Params:
;   event: in, required, type=struct
;       event structure
;-
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
            widget_control, state.wMainPlot, draw_xsize=sizes.mainplot_x, $
                draw_ysize=sizes.mainplot_y
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
                widget = widget_info(state.wMainPlotRanges, $
                    find_by_uname=uname[i])
                widget_control, widget, xsize=sizes.mainplot_range
            endfor
        
            uname = [ 'XCUT' + [ 'XMIN', 'XMAX', 'YMIN', 'YMAX' ], $
                      'YCUT' + [ 'XMIN', 'XMAX', 'YMIN', 'YMAX' ]]
            for i = 0, n_elements(uname) - 1 do begin
                widget = widget_info(state.wCutRanges[i/4], $
                    find_by_uname=uname[i])
                widget_control, widget, xsize=sizes.cut_range
            endfor
            widget_control, event.top, update=1, /clear_events
            mainplot_redraw, state
            xcut_redraw, state
            ycut_redraw, state
            widget_control, event.top, SET_UVALUE=state, /no_copy
        end
        else: print, 'uncaught event'
    endcase
end

;------------------------------------------------------------------------------
;+
; Deallocate data pointers in state structure.
;
; :Params:
;    state : in, required, type=struct
;       Browse state structure
;-
pro free_browse_pointers, state
    compile_opt idl2, hidden

    ptr_free, state.data
    ptr_free, state.data_red
    ptr_free, state.mag
    ptr_free, state.mag_red
    ptr_free, state.x_red
    ptr_free, state.y_red
end

;------------------------------------------------------------------------------
;+
; Clean up procedure.
;
; :Params:
;    wBase: in, required, type=integer
;       widget id of the top level base.
;-
pro BrowseCleanup, wBase
    compile_opt idl2, hidden

    widget_control, wBase, get_uvalue=state, /no_copy
    free_browse_pointers, state
end

;------------------------------------------------------------------------------
;+
; Returns the current plot state. Used to save the plot state so that plot 
; windows created by the user won't be affected.
;   
;  :Returns: plot state structure
;-
function current_plot_state
    compile_opt idl2, hidden

    return, { p: !p, x: !x, y: !y, z: !z, w: !d.window }
end

;------------------------------------------------------------------------------
;+
; Resets the plot parameters to those given in plot_state.  See
; current_plot_state for more details.
;
; :Params:
;   plot_state: in, required, type=struct
;      plot state structure
;-
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
;+
; Initializes the state structure for the application. 
;
; :Params:
;   data: in, required, type=fltarr
;       dataset to browse
;   t_in: in, required, type=fltarr
;       t-axis coordinates
;   x_in: in, required, type=fltarr
;       x-axis coordinates
;   y_in: in, required, type=fltarr
;       y-axis coordinates
;   titles: in, required, type=strarr
;       axis titles
;   
; :Keywords:
;   p: in, optional, type=intarr(3)
;       dimension permutation to apply to the dataset; relative to the order
;       when the application was first started.
;   old_p: in, optional, type=intarr(3)
;       the current dimension permutation of the dataset, as compared to its
;       state when the application was first started.
;-
function init_state, data, t_in, x_in, y_in, titles, p=p, old_p=old_p
    compile_opt idl2, hidden

    if n_elements(old_p) eq 0 then old_p = [ 0, 1, 2 ]
    if n_elements(p) ne 0 then begin
        ; find the permutation relative to the current state of the dataset
        p2 = (sort(old_p))[p]
        data = transpose(temporary(data), p2)
        ; reorder only t,x, and y titles
        titles = titles[[p2, 3]]
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

    ; calculate vector lengths if 4D dataset
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

    ; create a decimated dataset if nx or ny is large
    max_dim = 200
    if (nx gt max_dim) || (ny gt max_dim) then begin
        reduced = 1
        x_red = congrid(x, nx < max_dim, /minus_one, /interp)
        y_red = congrid(y, ny < max_dim, /minus_one, /interp)
        if (size(data))[0] eq 4 then begin
            ndir = (size(data, /dimensions))[3]
            data_red = fltarr(nt, nx < max_dim, ny < max_dim, ndir, /nozero)
            for i = 0, ndir - 1 do begin
                data_red[*, *, *, i] = congrid(data[*, *, *, i], nt, $
                    nx < max_dim, ny < max_dim, /minus_one, /interp)
            endfor 
            mag_red = congrid(mag, nt, nx < max_dim, ny < max_dim, /minus_one, /interp)
        endif else begin
            data_red = congrid(data, nt, nx < max_dim, ny < max_dim, /minus_one, /interp)
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

    ; plot parameter structure
    axis_structure = { $
        log: 0, $               ; true if axis is log, false if axis is linear
        min: 0.0, max: 0.0 $    ; axis range
    }
    ; x,y,z axes for main plot
    mainplot_axes = replicate(axis_structure, 3)
    ; x, y axes for each cut plot.  first dimension selects the cut, second
    ; selects the axis.
    cut_axes = replicate(axis_structure, 2, 2)

    ; id is the slider id, value_id is the associated text field id,
    ; min/max_id are the widget_label ids that show min/max for the slider
    wSlider = replicate({ id: 0L, value_id: 0L, min_id: 0L, max_id: 0L }, 3)

    ; set dt such that playing the entire dataset takes approx 120 seconds.
    dt = nt / 120 > 1

    ; default surface plot view
    rotation = { xangle: 30., zangle: 30., rotating: 0, x0: 0.0, y0: 0.0 }

    state = {  $
        wMainPlot: 0L, $                       ; main plot id
        wMainPlotType: 0L, $                   ; main plot type list id
        wMainPlotToggles: 0L, $                ; main plot toggles base id
        wMainPlotRanges: 0L, $                 ; main plot ranges base id
        wDims: lonarr(3), $                    ; dim. selector droplist ids
        wXcut: 0L, $                           ; xcut draw widget id
        wYcut: 0L, $                           ; ycut draw widget id
        wCutToggles: lonarr(2), $              ; x/y cut toggles base ids
        wCutRanges: lonarr(2), $               ; x/y cut ranges base ids
        wSlider: wSlider, $                    ; slider id structures
        titles: titles, $                      ; strarr of axis titles
        mainplot_type: 0, $                    ; contour/surface/velovect/etc
        zmin: min(*data), zmax: max(*data), $  ; min/max of dataset
        data: data, $                          ; ptr to dataset
        mag: mag, $                            ; ptr to vector magnitudes
        mag_red: mag_red, $                    ; ptr to reduced dataset
        x: x, y: y, t: t, $                    ; axis coordinates
        xs: xs, ys: ys, ts: ts, $              ; sort(t/x/y)
        xi: 0, yi: 0, ti: 0, zi: zi, $         ; current index in each dim.
        p: p, $                                ; current dimension permutation
        max_mag: max_mag, $                    ; max vector mag at each timestep
        ref_mag: ref_mag, $                    ; max(max_mag)
        reduced: reduced, $                    ; use reduced dataset flag
        data_red: data_red, $                  ; ptr to reduced dataset
        x_red: x_red, y_red: y_red, $          ; axis coords for reduced dataset
        cut_axes: cut_axes, $                  ; plot params for cuts
        mainplot_axes: mainplot_axes, $        ; plot params for main plot
        iso: 0, $                              ; isotropic flag
        animating: 0, $                        ; animation on/off flag
        dt: dt, $                              ; dt between animation frames
        rotation: rotation $                   ; surface plot rotation state
    }
    return, state
end

;------------------------------------------------------------------------------
;+ 
; Resets the plot parameters, sliders, rotation state, etc. back to their
; defaults.  Does not change the main plot type or the dimension order.
;
; :Params:
;   state: in, required, type=struct
;       application state structure
;
; :Keywords:
;   no_redraw: in, optional, type=boolean
;       if set, do not redraw the plots
;-
pro reset_view, state, no_redraw=no_redraw
    compile_opt idl2, hidden

    state.mainplot_axes.log = intarr(3)
    state.mainplot_axes.min = [ min(state.x), min(state.y), state.zmin ]
    state.mainplot_axes.max = [ max(state.x), max(state.y), state.zmax ]

    state.rotation = { xangle: 30., zangle: 30., rotating: 0, x0: 0.0, y0: 0.0 }

    ; set axes to linear, ranges to min/max of the dataset.
    state.cut_axes.log = intarr(2,2)
    state.cut_axes.min = [ [ min(state.x), min(state.y) ], $
        [ state.zmin, state.zmin ] ]
    state.cut_axes.max = [ [ max(state.x), max(state.y) ], $
        [ state.zmax, state.zmax ] ]
    state.iso = 0
    state.xi = state.xs[0]
    state.yi = state.ys[0]
    state.ti = state.ts[0]
    state.mainplot_type = widget_info(state.wMainPlotType, /droplist_select)

    ; By default, scale the z-axis according to the vector magnitude if we have
    ; a velovect plot and 'MAG' is selected.
    if (state.mainplot_type eq 4) && (state.zi eq 0) then begin
        state.mainplot_axes[2].min = 0.
        state.mainplot_axes[2].max = state.ref_mag
        state.cut_axes[*, 1].min = 0.
        state.cut_axes[*, 1].max = state.ref_mag
    endif

    ; update the text fields
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
        widget = widget_info(state.wMainPlotToggles, find_by_uname=uname[i])
        widget_control, widget, set_button=0, sensitive=log_enabled[i]
    endfor

    uname = 'CONTOUR' + [ 'X', 'Y', 'Z' ] + 'MIN'
    for i = 0, n_elements(uname) - 1 do begin
        widget = widget_info(state.wMainPlotRanges, find_by_uname=uname[i])
        change_range_field, widget, state.mainplot_axes[i].min
    endfor

    uname = 'CONTOUR' + [ 'X', 'Y', 'Z' ] + 'MAX'
    for i = 0, n_elements(uname) - 1 do begin
        widget = widget_info(state.wMainPlotRanges, find_by_uname=uname[i])
        change_range_field, widget, state.mainplot_axes[i].max
    endfor

    uname = [ 'XCUT', 'YCUT' ]
    for cut = 0, n_elements(uname) - 1 do begin
        dir = [ 'X', 'Y' ]
        for d = 0, 1 do begin
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+dir[d]+'MIN')
            change_range_field, widget, state.cut_axes[cut, d].min
            widget = widget_info(state.wCutRanges[cut], $
                find_by_uname=uname[cut]+dir[d]+'MAX')
            change_range_field, widget, state.cut_axes[cut, d].max
        endfor
    endfor
    widget_control, state.wSlider[0].id, set_value=0
    widget_control, state.wSlider[1].id, set_value=0
    widget_control, state.wSlider[2].id, set_value=0
    state.ti = slider_index(state.wSlider[0], state.t, state.ts)
    state.yi = slider_index(state.wSlider[1], state.y, state.ys)
    state.xi = slider_index(state.wSlider[2], state.x, state.xs)
    change_slider_range, state.wSlider[0], state.t
    change_slider_range, state.wSlider[1], state.y
    change_slider_range, state.wSlider[2], state.x
    scale3
    if ~keyword_set(no_redraw) then begin
        mainplot_redraw, state
        xcut_redraw, state
        ycut_redraw, state
    endif

end

;------------------------------------------------------------------------------
;+
; Calculate sizes for the draw widgets and sliders based on the window size
;
;  :Returns: structure of sizes
;         
;  :Params:
;     x: in, optional, type=integer
;       horizontal size of window in pixels
;     y: in, optional, type=integer
;       vertical size of window in pixels
;-
function browse_widget_sizes, x, y
    ; default sizes. 
    if n_elements(x) eq 0 then x = 1139
    if n_elements(y) eq 0 then y = 727

    mainplot_x = .5267*x
    mainplot_y = .6878*y
    ; if the window width is very small, then shrink the text fields
    if x lt 1024 then field_delta = -1 else field_delta = 0

    return, {browse_widget_sizes, $
        mainplot_x: mainplot_x, $
        mainplot_y: mainplot_y, $
        cut_x: .3687*x, $
        cut_y: .3851*y, $
        tslider: 0.69*mainplot_x, $
        xslider: 0.8*mainplot_y, $ 
        yslider: 0.75*mainplot_x, $
        slider_field: 7 + field_delta, $
        mainplot_range: 8 + field_delta, $
        cut_range: 7 + field_delta $
    }
end

;------------------------------------------------------------------------------
;+
; Main application procedure. Parses arguments, creates the widgets and realizes
; the application window.
;
;  :Params:
;     data_in: in, required, type=fltarr or h5variable
;     t_in: in, optional, type=fltarr  
;        t axis coordinates
;     x_in: in, optional, type=fltarr  
;        x axis coordinates
;     y_in: in, optional, type=fltarr  
;        y axis coordinates
;  
;  :Keywords:
;     ttitle: in, optional, type=string
;       t axis title
;     xtitle: in, optional, type=string
;       x axis title
;     ytitle: in, optional, type=string
;       y axis title
;     data_title: in, optional, type=string
;       dataset title
;     p: in, optional, type=intarr(3)
;       array specifying how to permute the dataset dimensions
;     f: in, optional, type=boolean
;       set this flag to title the t axis as f (e.g., for FFTs)
;     h5select: in, optional, type=boolean
;       use the current selection on an h5variable dataset
;     help: in, optional, type=boolean
;       print a help message
;-
pro browse, data_in, t_in, x_in, y_in, scalar = scalar, $
        ttitle=ttitle, xtitle=xtitle, ytitle=ytitle, data_title=ztitle, $
        p=p, f=f, h5select=h5select, help=help

    compile_opt idl2, hidden
    on_error, 2

    if keyword_set(help) || (n_params() eq 0) then begin
        print
        print, 'Usage: browse, DATA, T, X, Y, [ named keywords ]'
        print, '    where DATA is a 3D or 4D array, or an h5variable ' $
            + 'object;'
        print, '    and T, X, and Y are 1D arrays of coordinates for ' $
            + 'each axis (optional).'
        print
        print, 'Named keywords: '
        print, '    SCALAR: Optional 3D/4D scalar dataset to be concatenated' $
            +  ' with the primary '
        print, '       dataset. The first three dimensions of SCALAR and ' $
            + 'DATA must match.'
        print, '    [TXY]TITLE: A string specifying the plot titles for ' $
            + 'each dimension.'
        print, '       Default titles are [''t'', ''x'', ''y''].'
        print, '    DATA_TITLE: A string specifying the title of the dataset.'
        print, '       Default is ''data''.'
        print, '    P: 3-element array specifying how to permute the ' $
            + 'dataset dimensions.'
        print, '       For example, [2, 1, 0] swaps the first and third ' $
            + 'dimensions.'
        print, '    /F: Set this flag to title the T axis as ''f'' (e.g., for ' $
            + 'FFTs).'
        print, '    /H5SELECT: Use the current selection on an h5variable ' $
            + 'dataset.'
        print, '    /HELP: Print this message.'
        print
        return
    endif

    case strupcase(!version.os_family) of
      'WINDOWS': set_plot, 'WIN'
      'MAC': set_plot, 'MAC'
      else: set_plot, 'X'
    endcase
    device, decomposed=0
    loadct, 39, /silent
    plot_state = current_plot_state()

    if n_elements(ztitle) eq 0 then begin
        ; if possible, use the IDL variable name of the dataset as the title
        ztitle = scope_varname(data_in, level=-1) 
        if ztitle eq '' then begin 
            ztitle = 'data' 
            window_title = 'Data Browser'
        endif else begin
            window_title = 'Data Browser - ' + ztitle
        endelse
    endif else begin
        window_title = 'Data Browser - ' + ztitle
    endelse

    widget_control, /hourglass

    if n_elements(scalar) eq 0 then begin
        nscalar = 0
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
    endif else begin
        ; flag if h5variable objects need to be restored
        h5variable = 0

        ; First, calculate the dimensions of the scalar dataset and create the appropriate IDL code to read it in.
        if (size(scalar, /type) eq 11) && $
                (obj_class(scalar) eq 'H5VARIABLE') then begin
            h5variable = 1
            scalar_dims = scalar->dims()
            scalar_string = 'reform(scalar->r(select=scalarselect), product(scalar_dims), /overwrite)'
        endif else if arg_present(scalar) then begin
            scalar_dims = size(scalar, /dimensions)
            scalar_string = 'reform(scalar, product(scalar_dims))'
        endif else if ~arg_present(scalar) then begin
            scalar_dims = size(scalar, /dimensions)
            scalar_string = 'reform(temporary(scalar), product(scalar_dims))'
        endif

        ; Do the same for the vector dataset.
        if (size(data_in, /type) eq 11) && $
                (obj_class(data_in) eq 'H5VARIABLE') then begin
            h5variable = 1
            data_dims = data_in->dims()
            data_string = 'reform(data_in->r(select=h5select), product(data_dims), /overwrite)'
        endif else if arg_present(data_in) then begin
            data_dims = size(data_in, /dimensions)
            data_string = 'reform(data_in, product(data_dims))'
        endif else if ~arg_present(data_in) then begin
            data_dims = size(data_in, /dimensions)
            data_string = 'reform(temporary(data_in), product(data_dims))'
        endif
        data_dims = data_dims[where(data_dims ne 1)]
        scalar_dims = scalar_dims[where(scalar_dims ne 1)]
        if n_elements(scalar_dims) eq 3 then $
            nscalar = 1 $
        else $
            nscalar = scalar_dims[3]


        if n_elements(data_dims) le 3 then $
            message, 'SCALAR keyword only allowed for vector (4D) datasets.'
        if ~array_equal(data_dims[0:2], scalar_dims[0:2]) then $
            message, 'SCALAR array dimensions do not match DATA array dimensions.'
  
        if h5variable then $
            message, 'Restoring H5VAR dataset...', /info

        success = execute('data = [' + data_string + ', ' + scalar_string + ']')
        if ~success then $
            message, 'Error combining DATA and SCALAR variables.'  
        data_dims[3] += nscalar
        data = reform(data, data_dims, /overwrite)
    endelse

    nonfinite = where(~finite(data), count)
    if count gt 0 then begin
        message, 'Data has non-finite elements such as NaN or Inf.'
    endif

    ; create axis coordinates if needed
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
        titles = [ 'f', 'x', 'y', ztitle ] $
    else $
        titles = [ 't', 'x', 'y', ztitle ]

    if n_elements(ttitle) ne 0 then titles[0] = ttitle
    if n_elements(xtitle) ne 0 then titles[1] = xtitle
    if n_elements(ytitle) ne 0 then titles[2] = ytitle
    if n_elements(ztitle) ne 0 then titles[3] = ztitle
    
    state = init_state(data, t_in[0:s[0]-1], x_in[0:s[1]-1], y_in[0:s[2]-1], $
        titles, p=p)

    sizes = browse_widget_sizes()
    nt = n_elements(state.t)
    nx = n_elements(state.x)
    ny = n_elements(state.y)

    ; create the main widget bases
    wBase = widget_base(column=1, title=window_title,  $
        /tlb_size_events, kill_notify='BrowseCleanup')
    wToolbarBase = widget_base(wBase, /base_align_center, /align_left, /row)
    wPlotsBase = WIDGET_BASE(wBase, column=2, /BASE_ALIGN_center, /align_left)
    wMainPlotCol = widget_base(wPlotsBase, /column)
    wCutCol = widget_base(wPlotsBase, /column)

    ; populate the toolbar
    wAnimateBase = widget_base(wToolbarBase, /base_align_center, /row, frame=1)
    wDTLabel = widget_label(wAnimateBase, value='dt:')
    wDTValue = widget_text(wAnimateBase, xsize=5, /edit, event_pro='dtValue', $
        /kbrd_focus_events, value=string(state.dt, format='(i0)'))
    wPlayButton = widget_button(wAnimateBase, value='Play', $
        event_pro='PlayButton')
    ; create the dimension selector
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
                if nscalar eq 0 then $
                    dimList = [ 'MAG', string(indgen((size(*state.data, $
                    /dimensions))[3]), format='(i0)') ] $
                else if nscalar eq 1 then $
                    dimList = [ 'MAG', string(indgen(data_dims[3]-nscalar), $
                    format='(i0)'), 'S'] $
                else $
                    dimList = [ 'MAG', string(indgen(data_dims[3]-nscalar), $
                    format='(i0)'), 'S' + string(indgen(nscalar), $
                    format='(i0)') ] 
                    
                state.wDims[0] = widget_droplist(wDimsBase, value = dimList, $
                    event_pro='DimList4D', title='DATA:')    
           end
        else: message, 'DATA must have 2, 3, or 4 dimensions.'
    endcase
    ; create the plot type droplist
    wMainPlotButtonsBase = widget_base(wToolbarBase, /base_align_center, $
        /row, frame=1)
    if ndims eq 4 then begin
        state.wMainPlotType = widget_droplist(wMainPlotButtonsBase, $
            value=['Contour', 'Level', 'Surface', 'ShadeSurface', 'Vector', 'Vector+Contour'], $
            event_pro='MainPlotTypeList')
        widget_control, state.wMainPlotType, set_droplist_select=4
        state.mainplot_type = 4
    endif else begin
        state.wMainPlotType = widget_droplist(wMainPlotButtonsBase, $
            value=['Contour', 'Level', 'Surface', 'ShadeSurface'], $
            event_pro='MainPlotTypeList')
        state.mainplot_type = 0
    endelse
    ; color table button
    wCTBButton = widget_button(wMainPlotButtonsBase, value='Color Table', $
        event_pro='CTBButton')
    ; export options buttons
    wExportBase = widget_base(wToolbarBase, /base_align_center, /row, frame=1)
    wExportLabel = widget_label(wExportBase, value='Export to:')
    wPNGButton = widget_button(wExportBase, value='PNG', event_pro='PNGButton')
    wVARButton = widget_button(wExportBase, value='Variable', $
        event_pro='VARButton')
    ; reset button
    wResetButton = widget_button(wToolbarBase, value='Reset View', $
        event_pro='ResetButton')
 
    ; T slider
    wTsliderBase = widget_base(wMainPlotCol, /row, /align_left, $
        /base_align_center)
    wTSliderLabel = widget_label(wTSliderBase, value='Slice:')
    state.wSlider[0].value_id = widget_text(wTSliderBase, /edit, $
        xsize=sizes.slider_field, /kbrd_focus_events, $
        event_pro='TsliderValue')
    state.wSlider[0].min_id = widget_label(wTsliderBase, value='12345678')
    if ndims eq 2 then begin
        state.wSlider[0].id = widget_slider(wTsliderBase, drag=1, max=1, $
            sensitive=0, event_pro='tslider', xsize=sizes.tslider, $
            /suppress_value) 
        widget_control, state.wSlider[0].value_id, sensitive=0
    endif else $
        state.wSlider[0].id = widget_slider(wTsliderBase, drag=1, max=nt-1, $
            event_pro='tslider', xsize=sizes.tslider, /suppress_value)
    state.wSlider[0].max_id = widget_label(wTsliderBase, value='12345678')  
 
    ; main plot column
    dir = ['X', 'Y', 'Z' ]
    wCutSliderBase = lonarr(2)
    wMainPlotBase = widget_base(wMainPlotCol, /row, /base_align_center)
    ; main plot draw widget
    state.wMainPlot = WIDGET_DRAW(wMainPlotBase, XSIZE=sizes.mainplot_x, $         
        YSIZE=sizes.mainplot_y, /button_events, /motion_events, $
        event_pro='MainPlotMouseEvent')
    ; bases for cut sliders
    wCutSliderBase[0] = widget_base(wMainPlotBase, /column, /align_top, $
        /base_align_center)
    wCutSliderBase[1] = widget_base(wMainPlotCol, /row, /align_left, $
        /base_align_center)
    ; main plot toggles and ranges
    wMainPlotOptions = widget_base(wMainPlotCol, /row, /align_left, $
        /base_align_center)
    wMainPlotToggles = widget_base(wMainPlotOptions, /nonexclusive, $
        row=2, /base_align_center)          
    for i = 0, 2 do begin
        wLogButton = widget_button(wMainPlotToggles, value=dir[i]+'LOG', $
            uname='CONTOUR'+dir[i]+'LOG', event_pro='MainPlotLogButtons')
        if i eq 1 then $
            wIsoButton = widget_button(wMainPlotToggles, value='ISO ', $
                uname='ISO', event_pro='ISOButton') 
    endfor         
    wMainPlotRanges = widget_base(wMainPlotOptions, row=1, /base_align_center)
    wCurrentRange = widget_base(wMainPlotRanges, /row, /base_align_center)
    for i = 0, 2 do begin
        wRangeLabel = widget_label(wCurrentRange, value=dir[i]+':')
        wRangeMin = widget_text(wCurrentRange, /edit, $
            xsize=sizes.mainplot_range, uname='CONTOUR'+dir[i]+'MIN', $
            event_pro='MainPlotRanges', /kbrd_focus_events)
        wRangeMax = widget_text(wCurrentRange, /edit, $
            xsize=sizes.mainplot_range, uname='CONTOUR'+dir[i]+'MAX', $
            event_pro='MainPlotRanges', /kbrd_focus_events)
    endfor
    state.wMainPlotToggles = wMainPlotToggles
    state.wMainPlotRanges = wMainPlotRanges

    ; insert x and y cut sliders
    cut = 0 
    state.wSlider[cut+1].value_id = widget_text(wCutSliderBase[cut], /edit, $
        xsize=sizes.slider_field, /kbrd_focus_events, event_pro='XSliderValue')
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
        xsize=sizes.slider_field, /kbrd_focus_events, event_pro='YSliderValue')
  
    ; cut plots column
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
            event_pro='CutRanges')
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

pro restore_plot

    thisOS = !VERSION.OS_FAMILY
    thisOS = STRMID(thisOS, 0, 3)
    thisOS = STRUPCASE(thisOS)
    CASE thisOS of
        'MAC': SET_PLOT, thisOS
        'WIN': SET_PLOT, thisOS
        ELSE: SET_PLOT, 'X'
    ENDCASE

    !p.multi = 0
    !p.font = -1
end

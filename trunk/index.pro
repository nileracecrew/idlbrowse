function index, array, value
    nv = n_elements(value)
    subscript = lonarr(nv)
    for i = 0, nv - 1 do begin
        void = min(abs(array - value[i]), s)
        subscript[i] = s
    endfor
    return, reform(subscript)
end
function regrid, array, lon, lat, lon2, lat2, cubic=cubic, missing=missing
; Return array gridded onto new coordinates lon2,lat2.
; Otherwise, works like INTERPOLATE.

; Modified 10/04 to accomodate lons going either from -180 to +180 or
; 0-360, but not both in the same grid!
; lons and lats much each be in ascending order.

nxorig = n_elements(lon)
nyorig = n_elements(lat)
iy = interpol(findgen(nyorig), lat, lat2)

newarr = fltarr(n_elements(lon2),n_elements(lat2))

if min(lon) lt 0 and max(lon2) gt max(lon) then begin
    wh2 = where(lon2 le max(lon), cnt)
    if cnt gt 0 then begin
        ix = interpol(findgen(nxorig), lon, lon2[wh2])
        newarr[wh2,*] = $
          interpolate(array, ix, iy, /grid, cubic=cubic, missing=missing)
    endif
    wh2 = where(lon2 gt max(lon))
    wh = where(lon le 0)
    whs = where(lon gt 0)
    sh_lon = [lon[whs], lon[wh]+360.]
    sh_array = [array[whs,*],array[wh,*]]
    ix = interpol(findgen(nxorig), sh_lon, lon2[wh2])
    newarr[wh2,*] = $
      interpolate(sh_array, ix, iy, /grid, cubic=cubic, missing=missing)
endif else begin
    if min(lon2) lt 0 and max(lon) ge 180 then begin
        wh2 = where(lon2 ge 0, cnt)
        if cnt gt 0 then begin
            ix = interpol(findgen(nxorig), lon, lon2[wh2])
            newarr[wh2,*] = $
              interpolate(array, ix, iy, /grid, cubic=cubic, missing=missing)
        endif
        wh2 = where(lon2 lt 0)
        wh = where(lon ge 180)
        whs = where(lon lt 180)
        sh_lon = [lon[wh]-360., lon[whs]]
        sh_array = [array[wh,*],array[whs,*]]
        ix = interpol(findgen(nxorig), sh_lon, lon2[wh2])
        newarr[wh2,*] = $
          interpolate(sh_array, ix, iy, /grid, cubic=cubic, missing=missing)
    endif else begin ; no shifts necessary
        ix = interpol(findgen(nxorig), lon, lon2)
        newarr = $
          interpolate(array, ix, iy, /grid, cubic=cubic, missing=missing)
    endelse
endelse


return, newarr
end

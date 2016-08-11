function get_all_w, wap, lon, lat, maxlat, D=D, plev=plev, normval=normval, $
                    hus=hus, hdivprof=hdivprof, hnodivprof=hnodivprof, $
                    sflux=sflux, sfall=sfall, cl=cl, ocean=ocean

; get total mass flux profile for ascent+divergent columns, 
;  normalized by wap_hi in ascent columns

; Use cos(lat) weighting.  Keep only locations with good data at 925 hPa.
;  For profiles, below 925 hPa, use all good data available.

; wap_lo = avg at 850 and 700 (levels 2,3)
; wap_hi = avg from 600-400 (levels 4-6)

; NORMVAL is the normalization value used
; HUS is an input 3-D+time tracer array
; HDIVPROF is the mean profile of HUS in ascent+div columns
; HNODIVPROF is the mean in the remainder of the cool-ocean-region
; SFLUX is any 2-D+time field (e.g. surface flux)
; SFALL is the total  (not mean) of SFLUX over the whole region
; CL is the index of the lowest must-be-good level (default 2, or closest 
;   to 980 hPa if plev is supplied)

; if PLEV not supplied, assume standard levels.

l1 = (lon+360) mod 360
whl = where(l1 gt 200 or l1 lt 45)  ; cool-ocean region
wht = where(abs(lat) le maxlat)        ;
nmon = get_index(wap,4)

if not keyword_set(cl) then cl = 2  ; lowest must-be-good level
bl1 = 2  &  bl2 = 3   ; lo levels
tl1 = 4  &  tl2 = 5  &  tl3 = 6  ; hi levels

nlev = get_index(wap,3)
if keyword_set(plev) then begin  ; find nearest ones in given set
    bl1 = round(interpol(indgen(nlev), plev, 850.))
    bl2 = round(interpol(indgen(nlev), plev, 700.))
    tl1 = round(interpol(indgen(nlev), plev, 600.))
    tl2 = round(interpol(indgen(nlev), plev, 500.))
    tl3 = round(interpol(indgen(nlev), plev, 400.))
    cl = round(interpol(indgen(nlev), plev, 980.))
endif

; Calculate weighting, based on land and latitude
restore, '1x1_landmask.save'
land = regrid(float(l), glon, glat, lon, lat)
wgt = float(reform(abs(wap[whl,min(wht):max(wht),cl,*]) lt 100)) ; good values
if keyword_set(ocean) then $
  for i=0,nmon-1 do wgt[*,*,i] = wgt[*,*,i] * (1-(land[whl,*])[*,wht]) 
for j=0,n_elements(wht)-1 do wgt[*,j,*]=wgt[*,j,*]*cos(lat[wht[j]]/180*!PI)

w1 = (wap[whl,min(wht):max(wht),bl1,*]+wap[whl,min(wht):max(wht),bl2,*])/2
w2 = (wap[whl,min(wht):max(wht),tl1,*]+wap[whl,min(wht):max(wht),tl2,*]+wap[whl,min(wht):max(wht),tl3,*])/3

wh = where(w2 lt 0 and w1 lt 0 and wgt gt 0)   ; ascent
wh2 = where(w2 lt 0 and w1 lt w2 and wgt gt 0) ; div+asc

normval = total((w2*wgt)[wh])
D = total(((w1-w2)*wgt)[wh2]) / normval

if keyword_set(sflux) then begin
    s1 = sflux[whl,min(wht):max(wht),*]
    sfall = total((s1*wgt)[where(abs(w1) lt 100)])
endif

profile = fltarr(nlev) - 999.
if keyword_set(hus) then begin
    hdivprof = profile & hnodivprof = profile
    whnodiv = where((w2 ge 0 or w1 ge w2) and wgt gt 0)
endif

for i=nlev-1,0,-1 do begin
    w = wap[whl,min(wht):max(wht),i,*]
    if i lt cl then begin
        wh2 = where(w2 lt 0 and w1 lt w2 and wgt gt 0 $
                          and abs(wap[whl,min(wht):max(wht),i,*]) lt 100,cnt)
        whnodiv = where((w2 ge 0 or w1 ge w2) and wgt gt 0 $
                                and abs(wap[whl,min(wht):max(wht),i,*]) lt 100)
        if cnt eq 0 then return, profile
    endif
    profile[i] = total((w*wgt)[wh2]) / normval
    if keyword_set(hus) then begin
        h = hus[whl,min(wht):max(wht),i,*]
        hdivprof[i] = mean(h[wh2], weight=wgt[wh2])
        hnodivprof[i] = mean(h[whnodiv], weight=wgt[whnodiv])
    endif
endfor

return, profile

end

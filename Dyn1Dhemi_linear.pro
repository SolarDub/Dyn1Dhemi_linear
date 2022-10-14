pro Dyn1Dhemi_linear

;
;  Dynamo model to produce and sustain poloidal & toroidal fields
;  using a finite difference method
;

;set_plot,'win'

close, /ALL

Sym = 0                       ; Solution symmetry (0 = antisymmetric B (sym A); 1 = symmetric B (antisym A))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Spatial (latitudinal) grid ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

NumPoints = 51                  ; Number of points in grid
NminusOne = NumPoints-1

; Define spatial grid (zero to pi/2)
x    = (0.5*!PI/(NminusOne)) * findgen(NumPoints)   ; Radians
xdeg = 180.*x/!PI                                   ; Degrees

latplot = 60                                        ; Latitude at which to plot field/potential over time
minlat=min(abs(xdeg-latplot),ilat)                  ; Find index at which chosen latitude is positioned

deltaXi   = 1.0/(NminusOne)                         ; Grid resolution
deltaXiSq = deltaXi*deltaXi

;;;;;;;;;;;;;
; Time grid ;
;;;;;;;;;;;;;

Beta         = 0.25
deltaTau     = Beta*deltaXiSq     ; Time-step

dTaudXi      = deltaTau/deltaXi
dTaudXisq    = deltaTau/deltaXiSq

endTau       = 2.                                       ; Time extent of simulation
NumTimeSteps = long(endTau*round(1/deltaTau))           ; Number of time-steps
Tau          = findgen(NumTimeSteps)*deltaTau           ; Time grid

print, 'Time-step = ', deltaTau
print, 'Number of time-steps = ', NumTimeSteps

;;;;;;;;;;;;;;;;;;
; Dynamo Numbers ;
;;;;;;;;;;;;;;;;;;

dynNumCrit = 390.        ; antisymmetric B simulation (linear case)

dynNumRatio = 1.
dynNum      = dynNumRatio*dynNumCrit  ; Dynamo number

dynRatio = 0.1   ; Ratio of cAlpha/cOmega - tends to control amplitude of a

; One of these needs to be negative to produce equatorward flow?
cAlpha = sqrt(dynNum*dynRatio)      ; Magnetic Reynolds Number (alpha-effect) positive due to sense of Coriolis force
cOmega = -1.*sqrt(dynNum/dynRatio)  ; Magnetic Reynolds Number (shear) negative to make Alpha(Grad(Omega)) X ePhi equatorward ?
print,'cAlpha = ', cAlpha, '  cOmega = ', cOmega


;;;;;;;;;;;;;;;;;
; Set up arrays ;
;;;;;;;;;;;;;;;;;

; Set up magnetic field parameter arrays
a = fltarr(NumPoints,NumTimeSteps)  ; Dimensionless vector potential along x-grid for each timestep
b = fltarr(NumPoints,NumTimeSteps)  ; Dimensionless magnetic field along x-grid for each timestep
f = fltarr(NumPoints)               ; Alpha-effect cut-off factor

; Set up arrays describing time-dependent field characteristics
aTau       = fltarr(NumTimeSteps)
bTau       = fltarr(NumTimeSteps)
aPeak      = fltarr(NumTimeSteps)
bPeakNorth = fltarr(NumTimeSteps)
bPeakSouth = fltarr(NumTimeSteps)
MagEn      = fltarr(NumTimeSteps)


; Initial conditions of field (a is already initialized to zero)
if (Sym eq 0) then begin
    b[*,0] = 0.5*sin(2.*x)         ; Symmetric (Stix eq. 11) - div by 2 to regulate field amplitudes
    b[NminusOne,0] = 0.0           ; Force value due to calculation rounding errors
endif else begin
    b[*,0] = sin(x)                ; Antisymmetric (Stix eq. 12)
endelse

;ld = sin(x)*sin(x)*cos(x)
ld = cos(x)                      ; Latitudinal dependence of helicity (See Stix full paragraph before eq. 3)
ld[NminusOne] = 0.0              ; Force value due to calculation rounding errors

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop over timesteps                                                            ;
; n is the present time index and n+1 is the next time index, holding the result ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for n=long(0),NumTimeSteps-2 do begin

    i = 0

    a[i,n+1]=0.0
    b[i,n+1]=0.0

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;  Loop over spatial grid ;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;

;   Find values from north pole to grid-point before equator

    for i=1,NumPoints-2 do begin

;       Poloidal part of B field (dA/dx)
        dAdx = (a[i+1,n]-a[i-1,n])*0.5*dTaudXi ; (Forward diff.)

        f[i] = 1.    ; Keeps equations linear

     ;  Use 2nd order finite differencing to calculate next a & b (Stix Eq 4 & 3, respectively)
        a[i,n+1] = a[i,n]                                                             $ ; Previous value
                 + cAlpha * ld[i] * f[i] * b[i,n] * deltaTau                          $ ; c_alph*f*b           - Alpha effect
                 + (a[i+1,n] - 2.*a[i,n] + a[i-1,n]) * dTaudXisq                        ; (d2A/dx2)            - Dissipative term

        b[i,n+1] = b[i,n]                                                             $ ; Previous value
                 + cOmega * dAdx                                                      $ ;  c_omega*(dA/dx)     - Shearing
                 + (b[i+1,n] - 2.*b[i,n] + b[i-1,n]) * dTaudXisq                        ; (d2B/dx2)            - Dissipative term

    endfor


; Find values at equator

    i=Numpoints-1          ; i.e. last grid point

    if (Sym eq 0) then begin   ; Antisymmetric (Stix Eq 9)
        a[i,n+1] = (4.*a[i-1,n+1] - a[i-2,n+1])/3.   ; using a one sided backward difference (2nd order) (with b = 0)
        b[i,n+1] = 0.0
    endif

    if (Sym eq 1) then begin   ; Symmetric (Stix Eq 10)
        a[i,n+1] = 0.0
        b[i,n+1] = (4.*b[i-1,n+1] - b[i-2,n+1])/3.   ; using a one sided backward difference (2nd order) (with a = 0)
    endif


    MagEn(n)=total((b[*,n]*b[*,n])/(8.*!PI))     ; Magnetic energy density, u = (B^2)/(2mu0)

    aTau[n]=a[ilat,n]                            ; Vector Potential at chosen latitude at given time
    bTau[n]=b[ilat,n]                            ; Magnetic Field at chosen latitude at given time

    peak=max(abs(b[0:NminusOne,n]),index)        ; Find spatial grid index at which Magnetic Field is a maximum
    bPeakNorth[n]=xdeg(index)                    ; Latitude at which Magnetic Field is a maximum

endfor


LOADCT,13
window,0
plot, Tau, MagEn $
    , TITLE='Magnetic Energy over spatial range', XTITLE='Time', YTITLE='Magnetic Energy' $
    , BACKGROUND= 255+255*256,COLOR=0    ; White background, black text

window,2
plot, Tau, bTau $
    , LINE=0, TITLE='Magnetic Field near '+string(latplot,format='(i2.2)')+' degrees latitude' $
    , XTITLE='Time', YTITLE='Magnetic Field' $
    , BACKGROUND=255+255*256, COLOR=0    ; White background, black text

window,3
plot, Tau, aTau $
    , LINE=0, TITLE='Vector Potential near '+string(latplot,format='(i2.2)')+' degrees latitude' $
    , XTITLE='Time', YTITLE='Vector Potential' $
    , BACKGROUND=255+255*256, COLOR=0    ; White background, black text


window,1
;contour,transpose(b), Tau, xdeg, C_COLORS=[50,150,250] $
contour,transpose(b), Tau, xdeg, C_COLORS=[125+125*50,125+255*256,255+255*256] $
    , TITLE='Northern Hemisphere Time-Latitude Butterfly Diagram', XTITLE='Tau', YTITLE='Latitude', /YSTYLE, YRANGE=[0,90] $
    , LEVELS=[-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7]
oplot, Tau, bPeakNorth, LINE=1, COLOR=120*256

window,4
plot, Tau, bPeakNorth $
    , TITLE='Position of Maximum Peak', XTITLE='Times', YTITLE='Latitude Of Peak', LINE=2, /YSTYLE, YRANGE=[0,90] $
    , BACKGROUND=255+255*256, COLOR=0    ; White background, black text

print,'Program end.'
end
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
COMMON SHARE5,time_h,dt_h,mass_h,ener_h,sx_h,sy_h,sz_h,bx_h,by_h,bz_h
COMMON SHAREERR,error

;-----------------------------------------------------------------------
; Procedure READBIN: Reads ATHENA binary dumps
;
PRO readbin,filename
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
;
; Read number of zones and variables
;
ndata=LONARR(2)
openr,1,filename
readu,1,ndata
nx=ndata[0]
nvar=ndata[1]
;
; Read (gamma-1) and isothermal sound speed
;
eos=fltarr(2)
readu,1,eos
gamm1=eos[0]
isocs=eos[1]
;
; Read grid coordinates
;
x=fltarr(nx)
readu,1,x
;
; Read data.
; nvar=4 means isothermal hydro.  nvar=5 means adiabatic hydro
; nvar=7 means isothermal MHD.    nvar=8 means adiabatic mhd
;
d =fltarr(nx)
e =fltarr(nx)
vx=fltarr(nx)
vy=fltarr(nx)
vz=fltarr(nx)
bx=fltarr(nx)
by=fltarr(nx)
bz=fltarr(nx)
readu,1,d
readu,1,vx
readu,1,vy
readu,1,vz
IF nvar eq 5 OR nvar eq 8 THEN readu,1,e
IF nvar eq 7 OR nvar eq 8 THEN BEGIN
  readu,1,by
  readu,1,bz
  readu,1,bx
ENDIF
vx=vx/d
vy=vy/d
vz=vz/d
IF gamm1 NE 0 THEN p=gamm1*(e-0.5*d*(vx^2+vy^2+vz^2)-0.5*(bx^2+by^2+bz^2))
IF gamm1 EQ 0 THEN p=isocs*isocs*d
;
close,1
END
;-----------------------------------------------------------------------
PRO anim_plot,nfiles,filename
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
;
; Load data into floating point array
;
nxwin = !D.X_VSIZE
nywin = !D.Y_VSIZE
readbin,filename
data=FLTARR(nx,nfiles)
data(*,0) = d(*)
FOR i=1,nfiles-1 DO BEGIN
  number=STRMID(filename,3,3)
  inumber=FIX(number)+1
  IF inumber LT 10 THEN BEGIN
    number=STRING(inumber,FORMAT='(I1)')
    STRPUT,filename,number,5
  ENDIF
  IF inumber GE 10 AND inumber LT 100 THEN BEGIN
    number=STRING(inumber,FORMAT='(I2)')
    STRPUT,filename,number,4
  ENDIF
  IF inumber GE 100 AND inumber LT 1000 THEN BEGIN
    number=STRING(inumber,FORMAT='(I3)')
    STRPUT,filename,number,3
  ENDIF
  IF inumber ge 1000 THEN BEGIN
    print,'ERROR in ANIM_PLOT, Invalid Filename'
    RETURN
  ENDIF
  readbin,filename
  data(*,i) = d(*)
ENDFOR
;
; Make plots, save as images
;
;  1D ANIMATION
dmax = MAX(data)
dmin = MIN(data)
images=BYTARR(nxwin,nywin,nfiles)
FOR i=0,nfiles-1 DO BEGIN
  PLOT,x,data[*,i],PSYM=6,YRANGE=[dmin,dmax]
  PLOT,x,data[*,i],YRANGE=[dmin,dmax]
plot,x,data[*,i],psym=6,symsize=.4,YRANGE=[dmin,dmax]
images[*,*,i]=TVRD(0,0,nxwin,nywin)
ENDFOR
;
;  2D ANIMATION
;dmax = MAX(data)
;dmin = MIN(data)
;images=BYTARR(nx,ny,nfiles)
;images[*,*,*] = BYTSCL(data[*,*,0,*],MIN=dmin,MAX=dmax)
;
;  3D ANIMATION OF FIELD LOOP
;dmax = MAX(data)
;dmin = MIN(data)
;images=BYTARR(nx,ny,nfiles)
;FOR i=0,nfiles-1 DO BEGIN
;  iz = ((14+i) MOD 29)
;  images[*,*,i] = BYTSCL(data[*,*,iz,i],MIN=dmin,MAX=dmax)
;ENDFOR
;
; Animate images
;
base=widget_base(title='Animation Widget')
animate = cw_animate(base,nxwin,nywin,nfiles)
widget_control, /realize, base
for i=0,nfiles-1 do cw_animate_load, animate, frame=i, image=images(*,*,i)
cw_animate_getp, animate, pixmap_vect
cw_animate_run, animate
xmanager, 'cw_animate demo', base, event_handler='ehandler'

END

pro ehandler, ev
widget_control, /destroy, ev.top
end
;-----------------------------------------------------------------------
PRO read_hist,nlines,filename
COMMON SHARE5,time_h,dt_h,mass_h,ener_h,sx_h,sy_h,sz_h,bx_h,by_h,bz_h

histdat=fltarr(10,nlines)
openr,3,filename
readf,3,histdat
close,3

time_h=histdat[0,*]
dt_h=histdat[1,*]
mass_h=histdat[2,*]
ener_h=histdat[3,*]
sx_h=histdat[4,*]
sy_h=histdat[5,*]
sz_h=histdat[6,*]
bx_h=histdat[7,*]
by_h=histdat[8,*]
bz_h=histdat[9,*]

END

;-----------------------------------------------------------------------
PRO six_plot
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz

!P.MULTI=[0,2,3,0,0]
dmin=MIN(d)
dmax=MAX(d)
;plot,x,d,psym=6,symsize=.4,YTITLE='D',YRANGE=[dmin,dmax]
plot,x,d,YTITLE='D',YRANGE=[dmin,dmax]
dmin=MIN(p)
dmax=MAX(p)
;plot,x,p,psym=6,symsize=.4,YTITLE='P',YRANGE=[dmin,dmax]
plot,x,p,YTITLE='P',YRANGE=[dmin,dmax]
dmin=MIN(vx)
dmax=MAX(vx)
;plot,x,vx,psym=6,symsize=.4,YTITLE='Vx',YRANGE=[dmin,dmax]
plot,x,vx,YTITLE='Vx',YRANGE=[dmin,dmax]
dmin=MIN(vy)
dmax=MAX(vy)
;plot,x,vy,psym=6,symsize=.4,YTITLE='Vy',YRANGE=[dmin,dmax]
plot,x,vy,YTITLE='Vy',YRANGE=[dmin,dmax]
dmin=MIN(by)
dmax=MAX(by)
;plot,x,by,psym=6,symsize=.4,YTITLE='By',YRANGE=[dmin,dmax]
plot,x,by,YTITLE='By',YRANGE=[dmin,dmax]
dmin=MIN(p/d)
dmax=MAX(p/d)
;plot,x,p/d,psym=6,symsize=.4,YTITLE='T',YRANGE=[dmin,dmax]
plot,x,p/d,YTITLE='T',YRANGE=[dmin,dmax]
!P.MULTI=0
END

;-----------------------------------------------------------------------
PRO four_plot
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz

!P.MULTI=[0,2,2,0,0]
dmin=MIN(d)
dmax=MAX(d)
plot,x,d,YTITLE='D',psym=6,symsize=.4,YRANGE=[dmin,dmax]
;plot,x,d,YTITLE='D',YRANGE=[dmin,dmax]
dmin=MIN(p)
dmax=MAX(p)
plot,x,p,YTITLE='P',psym=6,symsize=.4,YRANGE=[dmin,dmax]
;plot,x,p,YTITLE='P',YRANGE=[dmin,dmax]
dmin=MIN(vx)
dmax=MAX(vx)
plot,x,vx,YTITLE='Vx',psym=6,symsize=.4,YRANGE=[dmin,dmax]
;plot,x,vx,YTITLE='Vx',YRANGE=[dmin,dmax]
dmin=MIN(p/d)
dmax=MAX(p/d)
plot,x,p/d,YTITLE='T',psym=6,symsize=.4,YRANGE=[dmin,dmax]
;plot,x,p/d,YTITLE='T',YRANGE=[dmin,dmax]
!P.MULTI=0
END
;-----------------------------------------------------------------------
PRO nine_plot,flag
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz

!P.MULTI=[0,3,3,0,0]
dmin=MIN(d)
dmax=MAX(d)
IF (flag EQ 0) THEN plot,x,d,psym=6,symsize=.4,YTITLE='D',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,d,YTITLE='D',YRANGE=[dmin,dmax]
dmin=MIN(p)
dmax=MAX(p)
IF (flag EQ 0) THEN plot,x,p,psym=6,symsize=.4,YTITLE='P',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,p,YTITLE='P',YRANGE=[dmin,dmax]
dmin=MIN(e)
dmax=MAX(e)
IF (flag EQ 0) THEN plot,x,e,psym=6,symsize=.4,YTITLE='E',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,e,YTITLE='E',YRANGE=[dmin,dmax]
dmin=MIN(vx)
dmax=MAX(vx)
IF (flag EQ 0) THEN plot,x,vx,psym=6,symsize=.4,YTITLE='Vx',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vx,YTITLE='Vx',YRANGE=[dmin,dmax]
dmin=MIN(vy)
dmax=MAX(vy)
IF (flag EQ 0) THEN plot,x,vy,psym=6,symsize=.4,YTITLE='Vy',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vy,YTITLE='Vy',YRANGE=[dmin,dmax]
dmin=MIN(vz)
dmax=MAX(vz)
IF (flag EQ 0) THEN plot,x,vz,psym=6,symsize=.4,YTITLE='Vz',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,vz,YTITLE='Vz',YRANGE=[dmin,dmax]
dmin=MIN(by)
dmax=MAX(by)
IF (flag EQ 0) THEN plot,x,by,psym=6,symsize=.4,YTITLE='By',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,by,YTITLE='By',YRANGE=[dmin,dmax]
dmin=MIN(bz)
dmax=MAX(bz)
IF (flag EQ 0) THEN plot,x,bz,psym=6,symsize=.4,YTITLE='Bz',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,bz,YTITLE='Bz',YRANGE=[dmin,dmax]
phi = atan(bz/by)
dmin=MIN(phi)
dmax=MAX(phi)
IF (flag EQ 0) THEN plot,x,phi,psym=6,symsize=.4,YTITLE='PHI',YRANGE=[dmin,dmax]
IF (flag EQ 1) THEN plot,x,phi,YTITLE='PHI',YRANGE=[dmin,dmax]
!P.MULTI=0
END
;-----------------------------------------------------------------------
; Procedure FLINES: plots 2D field lines
;
PRO flines,nlev
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
dims = SIZE(d,/dimensions)
vecpot=fltarr(dims[0],dims[1])
dx = x[1]-x[0]
dy = y[1]-y[0]
vecpot[0,0] = 0.0
FOR J=1,dims[1]-1 DO vecpot[0,J] = vecpot[0,j-1] + dy*bx[0,j]
FOR I=1,dims[0]-1 DO vecpot[I,*] = vecpot[i-1,*] - dx*by[i,*]
contour,vecpot,nlevels=nlev
END
;-----------------------------------------------------------------------
PRO shu_osher_plot
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
readbin,'bin047aa_shuosher_800'
plot,x,d,YTITLE='D'
readbin,'bin047aa_shuosher_200'
oplot,x,d,psym=6,symsize=.4
END
;-----------------------------------------------------------------------
PRO sod_plot
COMMON SHARE1,nx,nvar,x
COMMON SHARE3,gamm1,isocs
COMMON SHARE4,d,e,p,vx,vy,vz,bx,by,bz
time=0.25
vs = 1.7522
vc = 0.92745
vf = 0.07027
vh = 1.1832
xs = 0.5+vs*time
xc = 0.5+vc*time
xf = 0.5-vf*time
xh = 0.5-vh*time
dsod=FLTARR(500)
vsod=FLTARR(500)
esod=FLTARR(500)
psod=FLTARR(500)
xsod=FLTARR(500)
xsod[0] = 0.0
FOR I=1,499 DO xsod[I]=xsod[I-1] + 1./500.
FOR I=0,499 DO BEGIN
  IF xsod[I] GT xs THEN BEGIN
    dsod[I] = 0.125
    psod[I] = 0.1
    vsod[I] = 0.0
    esod[I] = 2.0
  ENDIF
  IF xsod[I] GT xc AND xsod[I] LT xs THEN BEGIN
    dsod[I] = 0.26557
    psod[I] = 0.30313
    vsod[I] = 0.92745
    esod[I] = 2.8535
  ENDIF
  IF xsod[I] GT xf AND xsod[I] LT xc THEN BEGIN
    dsod[I] = 0.42632
    psod[I] = 0.30313
    vsod[I] = 0.92745
    esod[I] = 1.7776
  ENDIF
  IF xsod[I] GT xh AND xsod[I] LT xf THEN BEGIN
    vsod[I] = 0.92745*(xsod[I]-xh)/(xf-xh)
    dsod[I] = 0.42632*(1.0+0.20046*(0.92745-vsod[I]))^5
    psod[I] = 0.30313*(1.0+0.20046*(0.92745-vsod[I]))^7
    esod[I] = psod[I]/(0.4*dsod[I])
  ENDIF
  IF xsod[I] LT xh THEN BEGIN
    dsod[I] = 1.0
    psod[I] = 1.0
    vsod[I] = 0.0
    esod[I] = 2.5
  ENDIF
ENDFOR
;
readbin,'bin025aa'
!P.MULTI=[0,2,2,0,0]
dmin=MIN(d)
dmax=MAX(d)
plot,x,d,YTITLE='D',psym=6,symsize=.4,YRANGE=[dmin,dmax]
oplot,xsod,dsod
dmin=MIN(p)
dmax=MAX(p)
plot,x,p,YTITLE='P',psym=6,symsize=.4,YRANGE=[dmin,dmax]
oplot,xsod,psod
dmin=MIN(vx)
dmax=MAX(vx)
plot,x,vx,YTITLE='Vx',psym=6,symsize=.4,YRANGE=[dmin,dmax]
oplot,xsod,vsod
dmin=MIN(p/d)
dmax=MAX(p/d)
plot,x,p/d,YTITLE='T',psym=6,symsize=.4,YRANGE=[dmin,dmax]
oplot,xsod,psod/dsod
!P.MULTI=0
END
;-----------------------------------------------------------------------
; Procedure  ERRORPLOT
;
PRO err_plot
COMMON SHAREERR,error
;
; PLot errors in fast wave
;
!P.MULTI=[0,3,1,0,0]
error = FLTARR(8,7)
openr,1,'err_fw2.dat'
readf,1,error
close,1
plot,error[0,*],error[1,*],/XLOG,/YLOG,XRANGE=[10.,2000.0],XSTYLE=1,YRANGE=[1.e-12,1.e-7],XTITLE='NZONES',YTITLE='L1 ERROR IN DENSITY'
oplot,error[0,*],error[1,*],psym=6
openr,1,'err_fw3.dat'
readf,1,error
close,1
oplot,error[0,*],error[1,*]
oplot,error[0,*],error[1,*],psym=6
;
; PLot errors in Alfven wave
;
openr,1,'err_aw3.dat'
readf,1,error
close,1
plot,error[0,*],error[3,*],/XLOG,/YLOG,XRANGE=[10.,2000.0],XSTYLE=1,YRANGE=[1.e-12,1.e-6],XTITLE='NZONES',YTITLE='L1 ERROR IN VY'
oplot,error[0,*],error[3,*],psym=6
openr,1,'err_aw2.dat'
readf,1,error
close,1
oplot,error[0,*],error[3,*]
oplot,error[0,*],error[3,*],psym=6
;
; PLot errors in slow wave
;
openr,1,'err_sw3.dat'
readf,1,error
close,1
plot,error[0,*],error[1,*],/XLOG,/YLOG,XRANGE=[10.,2000.0],XSTYLE=1,YRANGE=[1.e-11,1.e-6],XTITLE='NZONES',YTITLE='L1 ERROR IN DENSITY'
oplot,error[0,*],error[1,*],psym=6
openr,1,'err_sw2.dat'
readf,1,error
close,1
oplot,error[0,*],error[1,*]
oplot,error[0,*],error[1,*],psym=6
!P.MULTI=0
END


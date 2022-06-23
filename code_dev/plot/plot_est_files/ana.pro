; *******************************************
; ~/imaka/ana/imaka/imaka_cov/ana.pro
; *******************************************
;
; Journal File for mchun@opakapaka.ifa.hawaii.edu
; Working directory: /Users/mchun/students/REU/2020/Eden McEwen/datapipeline
; Date: Fri Jul 17 16:30:27 2020

;; det_thres is related to the strength of the layer.  A strong layer 
;; jumps out of acor while a weak layer does not (in sigma above background)
;;
;; vsd_thres is related to how much we're willing to allow the derived 
;; velocity to vary along the radial path.  This is a fractional threshold
;; Namely, Vsd is normalized by Vspd and checked that its smaller than Vsd_thres (fractional).

PRO find_layers, acovsm, outfile, $
	outfilepre=outfilepre, $	; text to add to beginning of each line (aocb filename, datetime, etc.)
	fsamp=fsamp, $ 			; sampling rate for this acov
	det_thres=det_thres, $		; detection threshold (minimium)
	vsd_thres=vsd_thres, $		; maximum variation in velocity along radial line
	layers=layers, $		; array that returns the found layers if any.  layers[nlayers,4] => layerDect, Vdir, Vspd, Vsd
	debug=debug, $	
	silent=silent			; can silence output graphics

	IF ( N_ELEMENTS(det_thres) EQ 0 ) THEN det_thres = 3.0
	IF ( N_ELEMENTS(vsd_thres) EQ 0 ) THEN vsd_thres = 0.3	;; fractional
	IF ( N_ELEMENTS(fsamp) EQ 0 ) THEN fsamp=180.
	IF ( N_ELEMENTS(outfilepre) EQ 0 ) THEN outfilepre=''

	IF ( KEYWORD_SET(debug) ) THEN BEGIN
		PRINT, 'find_layers: det_thres = ',STRING(det_thres,FORMAT='(F5.1)')
		PRINT, 'find_layers: vsd_thres = ',STRING(vsd_thres,FORMAT='(F5.1)')
	ENDIF
	x0 = 0
	y0 = 0
	r0 = 7
	rmax = 6	;; how far out from the center to look for peaks...edge is noisy
	dsub = 2.2/8.   ;; meters
	dmag = 4	;; magnification on the dsplay of the cov
	nmindetections = 5	;; minimum number of detections to cluster.  Same as nlayers in clustering?
	nmaxlayers = 5		;; maximum cluters to look for

	nt = (SIZE(acovsm))[3]
	ntcalc = nt ;; max dt to analyze here.	

	makexy, -r0,r0,1, -r0,r0,1, x, y
	r = SQRT(x^2+y^2)
	t = ATAN(y,x)

	; So calculate peaks above some thres=0 + sigma 
	; Do this for all annuli
	; Do this for each dt.  
	; The associated speed is r/dt
	; For each speed, look in each direction and ask if there is a detection in that 
	; direction in the appropriate dt 
	; For each pixel in the annulus, calculate its detection level in sigma over/under 
	; the mean in that annulus 
	
	IF ( KEYWORD_SET(debug) ) THEN writefits, 'covmap.fits', acovsm

	det = FLTARR(15*15,nt)
	meanarray = FLTARR(15*15,nt)
	stddevarray = FLTARR(15*15,nt)
	FOR dt=0,ntcalc-1 DO BEGIN

    		;; im is the "current" time slice.  
    		im = REFORM(acovsm[*,*,dt])
    		;; Consider data in annuli of r=1pix to r=7pix
    		FOR dr=1,rmax DO BEGIN 
        		; Pick out pixels in an annulus
        		idx = WHERE( dr GE r-0.5 AND dr LE r+0.5 )
        		; Calculate the mean and stddev within this annulus.
        		; Sigma-clip pixels in annulus
        		MEANCLIP, im[idx], m, sd, CLIPSIG=3	;; default clipsig=3, MAXITER=5, 
			; Eden saves the mean and stddev arrays - might help to compare...
			; Send Eden printout (longer times) and send cubes like she carries along
			;print, dt, dr, m, sd
			meanarray[idx,dt] = m
			stddevarray[idx,dt] = sd
        		;; At each location assign a detection level = deviation over mean in sd
        		det[idx,dt] = (im[idx]-m)/sd ;; detection is in units of sigma over mean
    		ENDFOR

    		IF (KEYWORD_SET(silent) EQ 0 AND EVEN(dt)) THEN BEGIN
;    			TVSCL, REBIN(im,dmag*15,dmag*15, /SAMPLE), dt*dmag*15/2, 0
;			TVSCL, BYTSCL(REBIN(REFORM(det[*,dt],15,15),dmag*15,dmag*15, $
;				/SAMPLE),MIN=0,MAX=4), dt*dmag*15/2,dmag*15
		ENDIF
	ENDFOR
	det = REFORM(det,15,15,nt)
	meanarray = REFORM(meanarray,15,15,nt)
	stddevarray = REFORM(stddevarray,15,15,nt)
	;;	 
	IF ( KEYWORD_SET(debug) ) THEN WRITEFITS, "detmap.fits", det
	IF ( KEYWORD_SET(debug) ) THEN WRITEFITS, "mnmap.fits", meanarray
	IF ( KEYWORD_SET(debug) ) THEN WRITEFITS, "sdmap.fits", stddevarray

	;; 2018-08-13 3') for each point in det(x,y,dy) above a threshold, calculate vsi, vdi, deti'
	vspdi = [-1.]
	vdiri = [-1.]
	vdeti = [-1.]
	vxposi = [-1.]
	vyposi = [-1.]
	vxedgei = [-1.]
	vyedgei = [-1.]
	vti = [-1.]
	vrposi = [-1.]
	FOR i=0,14 DO FOR j=0,14 DO BEGIN
		FOR k=1,ntcalc-1 DO BEGIN
			IF ( det[i,j,k] GE det_thres ) THEN BEGIN
				dri = SQRT( (x[i,j]-x0)^2 + (y[i,j]-y0)^2 ) * dsub
				dti = k / fsamp
				thetai = ATAN(y[i,j], x[i,j])
				vspdi = [vspdi,dri/dti]
				vdiri = [vdiri, thetai]
				vdeti = [vdeti, det[i,j,k]]
				vxposi = [vxposi, x[i,j]]
				vyposi = [vyposi, y[i,j]]
				vxedgei = [vxedgei, r0*cos(thetai)]	;; edge intersect of this direction
				vyedgei = [vyedgei, r0*sin(thetai)] ;; use for clustering
				vti = [vti, dti]
				vrposi = [vrposi, dri]
			ENDIF
		ENDFOR
	ENDFOR
	IF ( N_ELEMENTS(vspdi) GT nmindetections ) THEN BEGIN

junk = WHERE(vrposi GE 2 )
IF ( junk[0] NE -1 ) THEN STOP
		vspdi = vspdi[1:*]
		vdiri = vdiri[1:*]
		vdeti = vdeti[1:*]
		vxposi = vxposi[1:*]
		vyposi = vyposi[1:*]
		vxedgei = vxedgei[1:*]
		vyedgei = vyedgei[1:*]
		vti = vti[1:*]
		vrposi = vrposi[1:*]
		;array = TRANSPOSE([ [vspdi], [vdiri] ])  
		array = TRANSPOSE([ [vspdi], [vxedgei], [vyedgei] ])  

		IF ( KEYWORD_SET(debug) ) THEN BEGIN
			OPENW, lun,  '/tmp/ana.cluster.txt', /GET_LUN
			FOR i=0,N_ELEMENTS(vspdi)-1 DO PRINTF, lun, vspdi[i], vdiri[i], vxedgei[i], vyedgei[i], FORMAT='(4F8.4)'
			FREE_LUN, lun
		ENDIF

		w = CLUST_WTS(array, N_CLUSTERS=nmaxlayers)
		res = CLUSTER(array, w, N_CLUSTERS=nmaxlayers)
		r2 = res(SORT(res))
		r2 = r2[UNIQ(r2)]
		nlayers = N_ELEMENTS(r2)
		IF ( KEYWORD_SET(SILENT) EQ 0 ) THEN BEGIN
			PLOT, [-180,180], [0, 40], pSYM=3
			CLUT = (FINDGEN(nlayers)+1)/FLOAT(nlayers) * 250.
			LOADCT, 39, /SILENT
		ENDIF
		layers = FLTARR(nlayers,4)	;; layer detection, Vdir, Vspd, Vsd
		OPENW, 1, outfile, /APPEND
		FOR i=0,nlayers-1 DO BEGIN
		
			idx = WHERE( res EQ r2[i] )
			vsi = vspdi[idx]
			vdi = vdiri[idx]/!DTOR
			vxi = vxedgei[idx]
			vyi = vyedgei[idx]

			layerDect = AVG(vdeti[idx])
			layerVsd = STDDEV(vsi)
			layerVspd = AVG(vsi)
			layerVdir = ATAN(AVG(vyi),AVG(vxi))/!DTOR
			layers[i,*] = [layerDect, layerVdir, layerVspd, layerVdir]

    			IF ( layerDect GE det_thres AND layerVsd/layerVspd LE vsd_thres ) THEN BEGIN
        		IF ( KEYWORD_SET(SILENT) EQ 0 ) THEN BEGIN
				OPLOT, vdi, vsi, PSYM=4, SYMSIZE=0.3, COLOR=CLUT[i]
				OPLOT, [1,1]*layerVdir, [1,1]*layerVspd, PSYM=1, SYMSIZE=5, THICK=2, COLOR=CLUT[i]
				PRINT, layerDect, layerVdir, layerVspd, layerVsd, FORMAT = '(" LayerDect=", F8.2, " Vdir=",F8.1," Vspd=", F8.1," Vsd=", F8.2)'
			ENDIF
        		PRINTF, 1, outfilepre, layerDect, layerVdir, layerVspd, layerVsd, FORMAT = '(A,"  ", F8.2, F8.1, F8.1, F8.2)'
    			ENDIF

		ENDFOR
		CLOSE, 1
	ENDIF

	;; END OF CLUSTERING APPROACH

;STOP
RETURN
END

;------------------------------------------------------------------------------

FUNCTION readacov, filename, h, style=style, dtmax=dtmax

        CASE style OF
                1: BEGIN ;; MC code
                        acov = readfits(filename, h)
                        nt = 1.*(SIZE(acov))[3]
                        nwfs = SXPAR(h,'NWFS')
                        acov = TOTAL(acov,4)/nwfs
                        END
                ELSE: BEGIN ;; Eden style
                        h = headfits(filename) 
                        nwfs = SXPAR(h,'NWFS')
                        nt = SXPAR(h,'TMAX')
                        dt = nt < dtmax
                        FOR i=0,nwfs-1 DO BEGIN
                                d = mrdfits(filename, i+1)
                                IF ( i EQ 0 ) THEN acov = TOTAL(d[*,*,0:dt-1,i,*],5)/2. ELSE acov = acov + TOTAL(d[*,*,0:dt-1,i,*],5)/2.
                        ENDFOR
                        FOR i=0,dt-1 DO BEGIN
                                acov[*,*,i] = ROTATE(REFORM(acov[*,*,i]),2)
                        ENDFOR
                        END
        ENDCASE 

RETURN, acov
END

; IDL Version 8.4 (linux x86_64 m64)
; Journal File for mchun@ehu.ifa.hawaii.edu
; Working directory: /home/mchun/eden/inject
; Date: Fri Jun  3 09:04:14 2022
 
FUNCTION readlayers, inputfile
        READCOL, inputfile, fname, layDect, layDir, laySpd, layVsd, FORMAT='(A,F,F,F,F)'
        n = N_ELEMENTS(fname)
        vindir = FLOAT( STRMID(fname,STRPOS(fname[0],'vdir_')+5,1) )
        vinspd = FLOAT( STRMID(fname,STRPOS(fname[0],'vspd_')+5,5) )
        vinr0 = FLOAT( STRMID(fname,STRPOS(fname[0],'r0_')+3,4) )
        hlp, fname, layspd, laydir, vinspd, vindir, vinr0

        layer = {layer_elements, fname: '', vindir: 0.0, vinspd:0.0, vinr0:0.0, $
                        layDect: 0.0, layDir:0.0, laySpd:0.0, layVsd:0.0 }
        layers = REPLICATE({layer_elements}, n)

        layers.fname = fname
        layers.vindir = vindir
        layers.vinspd = vinspd
        layers.vinr0 = vinr0
        layers.layDect = layDect
        layers.layDir = layDir
        layers.laySpd = laySpd
        layers.layVsd = layVsd

RETURN, layers
END


PRO file2layers, fname, dtmax, nwfsmax, $
	det_thres=det_thres, $
	vsd_thres=vsd_thres, $
	outfile=outfile, $
	xcov=xcov, $
	silent=silent, debug=debug

;STOP
    	h = headfits(fname, /SILENT)
    	datetime = SXPAR(h, 'DATETIME')
	nwfs = SXPAR(h,'NWFS')
	IF ( nwfs GT nwfsmax ) THEN nwfs=nwfsmax

	fname2 = STRSPLIT(fname,'/',/EXTRACT)
	fname2 = fname2[N_ELEMENTS(fname2)-1]	;; pull of last bit
	dts = STRSPLIT(datetime,'-',/EXTRACT)
	dts = dts[0]+' '+dts[1]			;; reformat to remove '-'
	PRINT, fname2, dts, nwfs, FORMAT='(A," ",A," NWFS=",I2)'
	outfilepre = STRING(fname2, dts, FORMAT='(A," ",A," ")')

	h = headfits(fname, exten=1, /SILENT)
	nt = SXPAR(h,'NAXIS3')

	scov = FLTARR(15,15,nt,nwfs)
	fsamps = [-1]
	IF ( KEYWORD_SET(xcov) ) THEN BEGIN	;; Do xcov
		FOR j=0,nwfs-1 DO BEGIN 	
    			;.. loop thru all WFSs
    			d = mrdfits(fname, j+1, h, /SILENT)
			FOR k=j+1,nwfs-1 DO BEGIN
				fsamps = [fsamps, SXPAR(h,'FSAMPLE')]
    				si = (d[*,*,*,k,0] + d[*,*,*,k,1])/2.
    				scov[*,*,*,j] += si 
			ENDFOR
		ENDFOR
	ENDIF ELSE BEGIN 	;; Do autocov
		FOR j=0,nwfs-1 DO BEGIN
    			;.. loop thru all WFSs
    			d = mrdfits(fname, j+1, h, /SILENT)
			fsamps = [fsamps, SXPAR(h,'FSAMPLE')]
    			si = (d[*,*,*,j,0] + d[*,*,*,j,1])/2.
    			scov[*,*,*,j] = si
		ENDFOR
	ENDELSE
	fsamps = fsamps[1:*]
	cnt = FLOAT(N_ELEMENTS(fsamps))
	cov = TOTAL(scov,4)/cnt
	fsampsd = STDDEV(fsamps)
	fsamp = AVG(fsamps)

	;; Crop image to dtmax if dtmax < nt
	IF ( nt LT dtmax ) THEN dtmax=nt
	nt = dtmax
	;hlp, cov, dtmax, nt
	cov = cov[*,*,0:nt-1]
	;hlp, cov, dtmax, nt
	
	;; Subtract a mean image...
	mcov = total(cov,3)/FLOAT(nt)
	covsm = cov
	FOR j=0,nt-1 DO covsm[*,*,j] -= mcov

	IF ( fsampsd LE 1 ) THEN BEGIN
		find_layers, covsm, outfile, fsamp=fsamp, det_thres=det_thres, vsd_thres=vsd_thres, outfilepre=outfilepre, /silent
		IF ( KEYWORD_SET(debug) ) THEN BEGIN
			PRINT, fname2+'.cluster.txt'
			SPAWN, 'cp /tmp/ana.cluster.txt ./'+fname2+'.cluster.txt'
		ENDIF

	ENDIF ELSE BEGIN
		PRINT, 'skipping sampling rates not equal'
	ENDELSE
RETURN
END

PRO night2layers, ddir , $
	det_thres=det_thres, $
	vsd_thres=vsd_thres, $
	dtmax=dtmax, $
	nwfsmax=nwfsmax, $
	outfile=outfile, $
	xcov=xcov, $
	aocbgrade=aocbgrade, $
	silent=silent, debug=debug

	IF ( N_ELEMENTS(det_thres) EQ 0 ) THEN det_thres=3.0
	IF ( N_ELEMENTS(vsd_thres) EQ 0 ) THEN vsd_thres=0.3
	IF ( N_ELEMENTS(dtmax) EQ 0 ) THEN dtmax=200
	IF ( N_ELEMENTS(nwfsmax) EQ 0 ) THEN nwfsmax=3
	IF ( N_ELEMENTS(xcov) EQ 0 ) THEN xcov=0
	IF ( N_ELEMENTS(debug) EQ 0 ) THEN debug=0
	IF ( N_ELEMENTS(outfile) EQ 0 ) THEN outfile='/tmp/ANA.mc.out'

	SPAWN, 'ls '+ddir+'/*stt.fits', list
	n = N_ELEMENTS(list)

	IF ( n EQ 0 ) THEN BEGIN
		PRINT, 'night2layers:  No files found: '+ddir
		RETURN
	ENDIF

	;; filter for only files that have ratings of 1 in eden's conglomeration spreadsheet
	IF ( N_ELEMENTS(aocbgrade) NE 0 ) THEN BEGIN
		;res = read_csv('/home/emcewen/data/aocb_rating_20201111.csv', header=hdr)
		res = read_csv('/home/emcewen/code_dev/csv/corr_classified.csv', header=hdr)
		;; HARDCODED BUG - format of res file is hardcoded here.
		nres = N_ELEMENTS(res.field01)
		class = res.field02
		dates = res.field03
		files = res.field04
		grades = res.field09
		rating = INTARR(n)-1
		FOR i=0,n-1 DO BEGIN
			fi = STRSPLIT(list[i], '/', /EXTRACT)
			fi = fi[N_ELEMENTS(fi)-1]            
			jdx = [-1]
			FOR j=0,nres-1 DO BEGIN
				k = WHERE( STRPOS(fi,files[j]) NE -1 AND STRPOS(fi,STRTRIM(STRING(dates[j]),2)) NE -1 )
				IF ( k NE -1 ) THEN jdx = [jdx,j]
			ENDFOR
			IF ( N_ELEMENTS(jdx) GT 1 ) THEN rating[i] = grades[jdx[1]]	;; WHAT IF THERE IS MORE THAN ONE!
		ENDFOR
		idx = WHERE(rating EQ AOCBGRADE)
		IF ( idx[0] NE -1 ) THEN BEGIN
			list = list[idx]
			n = N_ELEMENTS(list)
		ENDIF ELSE BEGIN
			PRINT, 'night2layers:  No files found with rating = '+STRING(AOCBGRADE)
			RETURN
		ENDELSE
	ENDIF

	PRINT, 'night2layers:  Processing '+STRING(n)+' files.'

	IF ( KEYWORD_SET(SILENT) EQ 0 ) THEN BEGIN
		WINDOW, XSI=30*4*15, YSI=450+4*15*2,/FREE, TITLE = ddir
		!Y.MARGIN = [35,2]   
		!P.MULTI = [4,4,1]   
		!P.CHARSIZE=2
	ENDIF

	FOR i=0,n-1 DO BEGIN
		fname = list[i]
		file2layers, fname, dtmax, nwfsmax, det_thres=det_thres, vsd_thres=vsd_thres, outfile=outfile, xcov=0, silent=silent, debug=debug
	ENDFOR
	
	IF ( KEYWORD_SET(SILENT) EQ 0 ) THEN BEGIN
		!P.MULTI = 0
		!Y.MARGIN = [4,2]   
	ENDIF

END

PRO do_anatest

	dtA = 3.00	;; Good 
	vtA = 0.3	;; Good, note fractional
	;dtA = 0.0	;; save just about everything!
	;vtA = 100	;; save evrything...
	night2layers, "/data/emcewen/out/20180302/fits/", det_thres=dtA, vsd_thres=vtA, outfile='./ANA.acov.test.txt', /SILENT
	night2layers, "/data/emcewen/out/20180302/fits/", det_thres=dtA, vsd_thres=vtA, /xcov, outfile='./ANA.xcov.test.txt', /SILENT
	
	
STOP
RETURN
END


PRO do_ana, det=det, vsd=vsd
; det = 2.5, vsd = 0.3 (Good old way)
; det = 0.0, vsd = 100 (save everything)
; det=4.0, vsd=0.3 (2020-08-13)
; det=3.0, vsd=0.3 (2022-06-07)

	IF (N_ELEMENTS(det) EQ 0) THEN det=3.5	
	IF (N_ELEMENTS(vsd) EQ 0) THEN vsd=0.3	

	SPAWN, 'date +%Y%m%d', ymd
	aocbgrade = 1	; RATING_MC = 1`
	outfile1 = './ANA.xcov.'+ymd+'.txt'
	outfile2 = './ANA.xcov.'+ymd+'.cfht.txt'
	outfile3 = './ANA.xcov.'+ymd+'.gfs.txt'
	outfile4 = './ANA.acov.'+ymd+'.txt'
	outfile5 = './ANA.acov.'+ymd+'.cfht.txt'
	outfile6 = './ANA.acov.'+ymd+'.gfs.txt'
	outfileA = './ANA.acov.'+ymd+'.GL.txt'

	ddir     = '/data/emcewen/out/'
	dirs = ['20180301' ,'20180302' ,'20180525' ,'20180526' ,'20180527' ,'20180528' ,'20180531', $
		'20180601' ,'20180602' ,'20180603' ,'20180604' ,'20180821' ,'20180822', $
		'20181218' ,'20181219' ,'20181220' ,'20181221' ,'20181222' ,'20181223' ,'20181224', $
		'20190226' ,$
		'20200117', '20200119', '20200121' ,'20200122' ,'20200129', $
		'20210429' ,'20210430' ,'20210501' ,'20210502' ,'20210503' ] ;, $
;		'20210828' ,'20210829' ,'20210830', $
;		'20220120' ,'20220121' ,'20220122' ,'20220123' ,'20220124' ,'20220125']
	print, dirs
	n = N_ELEMENTS(dirs)
STOP

	IF ( 1 ) THEN BEGIN
		SPAWN, 'rm '+outfile1
		SPAWN, 'rm '+outfile4
		FOR i=0,n-1 DO BEGIN
			night2layers, ddir+dirs[i]+"/fits/", det_thres=det, vsd_thres=vsd, /xcov, outfile=outfile1, aocbgrade=grade, /SILENT 
			night2layers, ddir+dirs[i]+"/fits/", det_thres=det, vsd_thres=vsd, outfile=outfile4, aocbgrade=grade, /SILENT 
		ENDFOR
		SPAWN, 'rm '+outfile2
		cfht2aocb, outfile1, outfile2 
		SPAWN, 'rm '+outfile5
		cfht2aocb, outfile4, outfile5 
	ENDIF

	SPAWN, 'rm '+outfile3
	gfs2aocb, outfile1, outfile3 
	SPAWN, 'rm '+outfile6
	gfs2aocb, outfile4, outfile6 
	markgl, outfile1, outfile4, outfileA
STOP
END

PRO mkcfht

	READCOL, 'cfht/cfht-wx.2018.dat', y, m, d, hr, mn, vs, vd, temp, rh, p, x
	vs = vs * 0.5144	;; convert to m/s
	js = ymds2js(y,m,d,(hr+10.)*3600+mn*60.)
	cfht_js = [js]
	cfht_vs = [vs]
	cfht_vd = [vd]

	READCOL, 'cfht/cfht-wx.2019.dat', y, m, d, hr, mn, vs, vd, temp, rh, p, x
	vs = vs * 0.5144	;; convert to m/s
	js = ymds2js(y,m,d,(hr+10.)*3600+mn*60.)
	cfht_js = [cfht_js, js]
	cfht_vs = [cfht_vs, vs]
	cfht_vd = [cfht_vd, vd]
	
	READCOL, 'cfht/cfht-wx.2020.dat', y, m, d, hr, mn, vs, vd, temp, rh, p, x
	vs = vs * 0.5144	;; convert to m/s
	js = ymds2js(y,m,d,(hr+10.)*3600+mn*60.)
	cfht_js = [cfht_js, js]
	cfht_vs = [cfht_vs, vs]
	cfht_vd = [cfht_vd, vd]

	SAVE, filename='cfht-wx.dat', cfht_js, cfht_vs, cfht_vd	
RETURN
END

PRO cfht2aocb, fana, fout
	
	RESTORE, 'cfht-wx.dat'	;; cfht_js, cfht_vs, cfht_vd

	READCOL, fana, $
		fname, ymd, hms, det, vdir, vspd, sigma_vspd, $
		FORMAT='(A,A,A,F,F,F,F)'
	n = N_ELEMENTS(ymd)
	c_js = FLTARR(n)-99
	c_vs = FLTARR(n)-99
	c_vd = FLTARR(n)-99
	
	FOR i=0L,n-1 DO BEGIN
		y = FIX(STRMID(ymd[i],0,4))
		m = FIX(STRMID(ymd[i],4,2))
		d = FIX(STRMID(ymd[i],6,2))
		hr = FLOAT(STRMID(hms[i],0,2))
		mn = FLOAT(STRMID(hms[i],2,2))
		sc = FLOAT(STRMID(hms[i],4,2))
		js = ymds2js(y,m,d,hr*3600.+mn*60.+sc)

		idx = WHERE ( ABS(cfht_js-js) LE 5*60. )
		IF ( idx[0] NE -1 ) THEN BEGIN
			c_js[i] = AVG(cfht_js[idx])
			c_vs[i] = AVG(cfht_vs[idx])
			c_vd[i] = AVG(cfht_vd[idx])
		ENDIF
	ENDFOR

	OPENW, 1, fout
	FOR i=0,n-1 DO PRINTF, 1, c_js[i], c_vs[i], c_vd[i], FORMAT='(G16.8, F12.4, F12.4)'
	CLOSE, 1

RETURN
END

PRO gfs2aocb, fana, fout
	
	RESTORE, 'gfs.windsUT12.dat'	;;
	;;JS              DOUBLE    = Array[5878]
	;;WDIR            FLOAT     = Array[10, 5878]
	;;WSPD            FLOAT     = Array[10, 5878]
	gfs_js = js
	gfs_vs = wspd
	gfs_vd = wdir

	READCOL, fana, $
		fname, ymd, hms, det, vdir, vspd, sigma_vspd, $
		FORMAT='(A,A,A,F,F,F,F)'
	n = N_ELEMENTS(ymd)
	g_js = FLTARR(10,n)-99
	g_vs = FLTARR(10,n)-99
	g_vd = FLTARR(10,n)-99
	
	FOR i=0L,n-1 DO BEGIN
		y = FIX(STRMID(ymd[i],0,4))
		m = FIX(STRMID(ymd[i],4,2))
		d = FIX(STRMID(ymd[i],6,2))
		hr = FLOAT(STRMID(hms[i],0,2))
		mn = FLOAT(STRMID(hms[i],2,2))
		sc = FLOAT(STRMID(hms[i],4,2))
		js = ymds2js(y,m,d,hr*3600.+mn*60.+sc)

		idx = WHERE ( ABS( gfs_js-js) LE 12*60*60. )	;; 12 hrs...
		IF ( idx[0] NE -1 ) THEN BEGIN
			;; GFS has one point per night so should have one element in idx...
			g_js[*,i] = gfs_js[idx[0]]
			g_vs[*,i] = gfs_vs[*,idx[0]] * 0.5144	;; Convert from knots to m/s
			g_vd[*,i] = gfs_vd[*,idx[0]]
		ENDIF
	ENDFOR

	OPENW, 1, fout
	FOR i=0,n-1 DO PRINTF, 1, g_js[i], g_vs[*,i], g_vd[*,i], FORMAT='(G16.8, 10F12.4, 10F12.4)'
	CLOSE, 1

RETURN
END

PRO markgl, fxcor, facor, fout
;; identify GL in acor maps
; Date: Thu Aug 13 11:01:10 2020
	vmargin = 4.0 ; m/s
	dmargin = 30.0 ; deg

	READCOL, fxcor, xfname, xymd, xhms, xdet, xvdir, xvspd, xvsd, FORMAT='(A,A,A,F,F,F,F)'
	READCOL, facor, afname, aymd, ahms, adet, avdir, avspd, avsd, FORMAT='(A,A,A,F,F,F,F)'

	na = N_ELEMENTS(afname)
	gl = INTARR(na) ;; Is this a GL (e.g. in xcor?) default is no
	dv = FLTARR(na)-99
	dd = FLTARR(na)-99

	nd = N_ELEMENTS(xfname)
	FOR di = 0,nd-1 DO BEGIN
		idx = WHERE( afname EQ xfname[di] )
		print, xfname[di]
		dvi = ABS(xvspd[di]-avspd[idx])
		ddi = ABS(xvdir[di]-avdir[idx])
		dv[idx] = dvi
		dd[idx] = ddi
		jdx = WHERE( dvi LE vmargin AND ddi LE dmargin )
		IF ( jdx[0] NE -1 ) THEN BEGIN
			FORPRINT, ABS(xvspd[di] - avspd[idx[jdx]]), ABS(xvdir[di]-avdir[idx[jdx]])
			GL[idx[jdx]] = 1
		ENDIF
	ENDFOR
	HELP, WHERE( GL EQ 1 )

	OPENW, 1, fout
	;FOR i=0,na-1 DO PRINTF, 1, afname[i], aymd[i], ahms[i], adet[i], avdir[i], avspd[i], avsd[i], GL[i]
	FOR i=0,na-1 DO PRINTF, 1, afname[i], GL[i]
	CLOSE, 1
END

; PLOTS
;
; IDL Version 8.5.1 (linux x86_64 m64)
; Journal File for mchun@ehu.ifa.hawaii.edu
; Working directory: /home/mchun/eden
; Date: Mon Jul 27 11:37:13 2020
 

PRO do_plots

	f1 = 'ANA.acov.20220607.txt'
	f2 = 'ANA.acov.20220607.cfht.txt'
	f3 = 'ANA.acov.20220607.gfs.txt'
	fA = 'ANA.acov.20220607.GL.txt'

	;; Read in and analyze the ana output files
	READCOL, f1, $
		fname, ymd, hms, det, vdir, vspd, sigma_vspd, $
		FORMAT='(A,A,A,F,F,F,F)'
	;; vdir is angle North of East (regular x/y axes)
	;; Want the angle East of North
	vdir = (270 - vdir) mod 360	;; If E is along plus x
	;; Do we have a GL file?
	IF ( N_ELEMENTS(fA) NE 0 ) THEN BEGIN
		READCOL, fA, fnameGL, GL, FORMAT='(A,F)'
	ENDIF ELSE BEGIN
		GL = INTARR(N_ELEMENTS(fname))
	ENDELSE

	READCOL, f2, c_js, c_vs, c_vd, FORMAT='(D,F,F)'	
	js2ymds, c_js, yr, mn, dy, sc

	READCOL, f3, g_js, va,vb,vc,vd,ve,vf,vg,vh,vi,vj, da,db,dc,dd,de,df,dg,dh,di,dj, FORMAT='(D, F,F,F,F,F,F,F,F,F,F, F,F,F,F,F,F,F,F,F,F)'	
	; a=700mb, b=650, c=600, d=550, e=500, f=450, g=400, h=350, i=300, j=250mb
	; 3000masl, 3590,  4200,  4860,  5580,  6300,  7200,  8100,  9170,  10400masl

	js2ymds, g_js, gyr, gmn, gdy, gsc
	g_vs = vj
	g_vd = dj

STOP

	;WINDOW, 0, XSI=1200, YSI=450
	;!P.MULTI = [0,3,1]
	WINDOW, 0, XSI=1100*0.8, YSI=1600*0.8
	LOADCT, 39, /SILENT
	!P.MULTI = [0,2,3]
	!P.CHARSIZE = 2.5

	idx = WHERE( det GE 3 )			; All

	faidx = WHERE( det GE 3 AND GL EQ 0 )	; FA
	faidx = WHERE( det GE 3 AND GL EQ 0  AND (ABS(vdir-g_vd) LE 45) )	; FA
	faidx = WHERE( det GE 3 AND GL EQ 0  AND (ABS(vspd-g_vs)/vspd LE 0.5) )	; FA

	glidx = WHERE( det GE 3 AND GL EQ 1 )	; GL
	;glidx = WHERE( det GE 3 AND GL EQ 1 AND (ABS(vdir-40-c_vd) LE 20) )	; GL with dirs that match CFHT (e.g. within 45deg
	glidx = WHERE( det GE 3 AND GL EQ 1 AND (ABS(vdir-40-c_vd) GT 20) )	; GL with dirs that do not match CFHT
	;glidx = WHERE( det GE 3 AND GL EQ 1 AND (ABS(vspd-c_vs)/vspd LE 0.25) )	; GL with speeds that match CFHT
	help, idx, faidx, glidx

	bsz = 10.
	bins = (FINDGEN(100)*bsz)+bsz/2.
	PLOT, bins, histogram(vdir[idx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(idx)), PSYM=10, $
		TITLE = 'Distribution of Wind Directions', $
		XTITLE = 'Direction (Degrees E of N)', XRA=[0,360], $
		YTITLE = 'Fraction of detections', YRA=[0,0.20]
	OPLOT, bins, HISTOGRAM(vdir[faidx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(faidx)), PSYM=10, COLOR=250
	OPLOT, bins, HISTOGRAM(vdir[glidx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(glidx)), PSYM=10, COLOR=128
	OPLOT, bins, HISTOGRAM(c_vd[idx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(idx)), PSYM=10, LINE=2
	OPLOT, bins, HISTOGRAM(g_vd[idx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(idx)), PSYM=10, LINE=1
	AL_LEGEND, ['imaka AOCB', 'imaka(FA)', 'imaka(GL)', 'CFHT', 'GFS'], LINEs=[0,0,0,2,1], COLOR=[255,250,128,255,255], /RIGHT, CHARSIZE=1

	bsz = 1.
	bins = (FINDGEN(100)*bsz)+bsz/2.
	PLOT, bins, histogram(vspd[idx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(idx)), PSYM=10, $
		TITLE = 'Distribution of Wind Speeds', $
		XTITLE = 'Speeds (m/s)', XRA=[0,25], $
		YTITLE = 'Fraction of detections', YRA=[0,0.20]
	OPLOT, bins, HISTOGRAM(vspd[faidx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(faidx)), PSYM=10, COLOR=250
	OPLOT, bins, HISTOGRAM(vspd[glidx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(glidx)), PSYM=10, COLOR=128
	OPLOT, bins, HISTOGRAM(c_vs[idx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(idx)), PSYM=10, LINE=2
	OPLOT, bins, HISTOGRAM(g_vs[idx],BIN=bsz,MIN=0)/FLOAT(N_ELEMENTS(idx)), PSYM=10, LINE=1
	AL_LEGEND, ['imaka AOCB', 'imaka(FA)', 'imaka(GL)', 'CFHT', 'GFS'], LINEs=[0,0,0,2,1], COLOR=[255,250,128,255,255], /RIGHT, CHARSIZE=1

	res = HIST_2d(c_vd[glidx], vdir[glidx], bin1=10, bin2=10, min1=0, max1=360, min2=0, max2=360)
	plot_image, REBIN(res,(SIZE(res))[1]*10,(SIZE(res))[2]*10,/SAMPLE), MIN=0, MAX=0.15*N_ELEMENTS(glidx), $
		TITLE = 'Directions (deg E of N)', $
		XTITLE = 'CFHT Vdir (deg)', YTITLE = 'imaka AOCB Vdir (deg)'
	OPLOT, [0,360],[0,360], LINE=1

	res = HIST_2d(c_vs[glidx], vspd[glidx], bin1=1.0, bin2=1.0, min1=0, max1=25, min2=0, max2=25)
	plot_image, res, MIN=0, MAX=0.15*N_ELEMENTS(glidx), $
		TITLE = 'Comparison of Speeds', $
		XTITLE = 'CFHT Vspd (m/s)', YTITLE = 'imaka AOCB Vspd (m/s)'
	OPLOT, [0,360],[0,360], LINE=1

	res = HIST_2d( (g_vd[faidx]) mod 360, vdir[faidx], bin1=10, bin2=10, min1=0, max1=360, min2=0, max2=360)   
	plot_image, REBIN(res,(SIZE(res))[1]*10,(SIZE(res))[2]*10,/SAMPLE), MIN=0, MAX=0.15*N_ELEMENTS(faidx), $
		TITLE = 'Directions (deg E of N)', $
		XTITLE = 'GFS Vdir (deg)', YTITLE = 'imaka AOCB Vdir (deg)'
	OPLOT, [0,360],[0,360], LINE=1

	res = HIST_2d(g_vs[faidx], vspd[faidx], bin1=1, bin2=1, min1=0, max1=25, min2=0, max2=25)
	plot_image, res, MIN=0, MAX=0.15*N_ELEMENTS(faidx), $
		TITLE = 'Comparison of Speeds', $
		XTITLE = 'GFS Vspd (m/s)', YTITLE = 'imaka AOCB Vspd (m/s)'
	OPLOT, [0,360],[0,360], LINE=1
	!P.MULTI = 0
	!P.CHARSIZE = 1

	;SAVE_PNG, 'vwind.png'

STOP
END



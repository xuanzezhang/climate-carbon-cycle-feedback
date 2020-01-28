;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    year=ispan(1880,2017,1)
    nyear=dimsizes(year)

    fpath_CO2="./"
    file_CO2="CO2_global_mean_1850-2017.txt"
    data_CO2 = asciiread(fpath_CO2+file_CO2,(/nyear+30/),"float")
    data_CO2!0="time"
    data_CO2&time=ispan(1850,2017,1)

    fpath_CO2E="./"
    file_CO2E="CO2E_global_annual_cumulated_FF_plus_LUC_1850-2017.txt"
    data_CO2E = asciiread(fpath_CO2E+file_CO2E,(/nyear+30/),"float")
    data_CO2E!0="time"
    data_CO2E&time=ispan(1850,2017,1)

    file_CO2E_unc="CO2E_uncertainty_global_annual_cumulated_FF_plus_LUC_1850-2017.txt"
    data_CO2E_unc = asciiread(fpath_CO2E+file_CO2E_unc,(/nyear+30/),"float")
    data_CO2E_unc!0="time"
    data_CO2E_unc&time=ispan(1850,2017,1)

    data_CO2E_up=data_CO2E + data_CO2E_unc
    data_CO2E_dw=data_CO2E - data_CO2E_unc
    data_CO2E_up!0="time"
    data_CO2E_up&time=ispan(1850,2017,1)
    data_CO2E_dw!0="time"
    data_CO2E_dw&time=ispan(1850,2017,1)

    fpath_T="./"
    file_T="T_global_mean_1850-2017_from_HadCRUT4.nc"
    fin=addfile(fpath_T+file_T,"r")
    data_T=fin->T({time|1880:2017})

    data_T!0="time"
    data_T&time=year

    CO2  = data_CO2({time|1880:2017})
    CO2E = data_CO2E({time|1880:2017})
    T    = data_T({time|1880:2017})

    CO2E_up = data_CO2E_up({time|1880:2017})
    CO2E_dw = data_CO2E_dw({time|1880:2017})

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;****************************************
; Calculate harmonics           
;   finfo(0,:,:,nhar) - amplitude of each harmonic
;   finfo(1,:,:,nhar) - phase of each harmonic
;   finfo(2,:,:,nhar) - percent variance explained by each harmonic 
;****************************************

   nfrq=nyear/2
   K=ispan(1,nfrq,1)*1.0 
   frq=K/nyear
   print(frq)
   timescale=1.0/frq
   print(timescale)

   finfo_CO2  = fourier_info(CO2, nfrq, 1.)
   finfo_CO2E = fourier_info(CO2E, nfrq, 1.)
   finfo_T    = fourier_info(T, nfrq, 1.)        ; fourier_info works on *right* dim

   finfo_CO2E_up = fourier_info(CO2E_up, nfrq, 1.)
   finfo_CO2E_dw = fourier_info(CO2E_dw, nfrq, 1.)

   amp_CO2_in  = finfo_CO2(0,:)
   amp_CO2E_in = finfo_CO2E(0,:)
   amp_T_in    = finfo_T(0,:)

   amp_CO2E_up_in = finfo_CO2E_up(0,:)
   amp_CO2E_dw_in = finfo_CO2E_dw(0,:)

   amp_T_in    = where(amp_T_in.le.0.01,0.01,amp_T_in)

   amp_CO2=runave(amp_CO2_in,1,1)
   amp_CO2E=runave(amp_CO2E_in,1,1)
   amp_T=runave(amp_T_in,1,1)

   amp_CO2E_up=runave(amp_CO2E_up_in,1,1)
   amp_CO2E_dw=runave(amp_CO2E_dw_in,1,1)

   alpha = amp_CO2/amp_CO2E
   print(alpha)

   alpha_up = amp_CO2/amp_CO2E_up
   alpha_dw = amp_CO2/amp_CO2E_dw
   print(alpha_up)
   print(alpha_dw)

   alpha_avg= avg(alpha)
   alpha_up_avg= avg(alpha_up)
   alpha_dw_avg= avg(alpha_dw)
   print(alpha_avg)
   print(alpha_up_avg)
   print(alpha_dw_avg)
;   print(finfo_CO2(1,:))

    reg_dCO2dT=new((/nfrq/),"float")
    reg_dCO2dT!0="timescale"
    reg_dCO2dT&timescale=timescale

    reg_dCO2EdT=new((/nfrq/),"float")
    reg_dCO2EdT!0="timescale"
    reg_dCO2EdT&timescale=timescale

    reg_gamma=new((/nfrq/),"float")
    reg_gamma!0="timescale"
    reg_gamma&timescale=timescale

    reg_gamma_up=new((/nfrq/),"float")
    reg_gamma_up!0="timescale"
    reg_gamma_up&timescale=timescale

    reg_gamma_dw=new((/nfrq/),"float")
    reg_gamma_dw!0="timescale"
    reg_gamma_dw&timescale=timescale

    reg_dCO2EdT_up=new((/nfrq/),"float")
    reg_dCO2EdT_up!0="timescale"
    reg_dCO2EdT_up&timescale=timescale

    reg_dCO2EdT_dw=new((/nfrq/),"float")
    reg_dCO2EdT_dw!0="timescale"
    reg_dCO2EdT_dw&timescale=timescale

    reg_dCO2dT = amp_CO2/amp_T
    reg_dCO2EdT = amp_CO2E/amp_T
    reg_dCO2EdT_up = amp_CO2E_up/amp_T
    reg_dCO2EdT_dw = amp_CO2E_dw/amp_T

   ;;;;;  reg_dCO2EdT = reg_dCO2dT * (1+beta) + gamma

    rc= regline(reg_dCO2dT,reg_dCO2EdT)

    print(rc)

    rc_zeta  = rc@yave
    rc_eta   = rc@xave
    rc_beta  = rc -1
    rc_gamma = rc@yintercept

    data_scatter  = new ( (/2,dimsizes(reg_dCO2dT)/), typeof(reg_dCO2dT))
    data_scatter(0,:) = reg_dCO2EdT
    data_scatter(1,:) = rc*(reg_dCO2dT-rc@xave) + rc@yave
   ;;;;;  reg_gamma = reg_dCO2EdT - reg_dCO2dT * (1+beta)

   reg_gamma = reg_dCO2EdT - reg_dCO2dT * (1 + rc_beta)
   print(reg_gamma)
;;;;;;;;;;;;;;;;;;;;
   reg_dCO2EdT_70yr = reg_dCO2EdT({timescale|70})
   print(reg_dCO2EdT_70yr)
   reg_dCO2dT_70yr = reg_dCO2dT({timescale|70})
   print(reg_dCO2dT_70yr)
   reg_gamma_70yr = reg_gamma({timescale|70})*2.12
   print(reg_gamma_70yr)

   reg_gamma_avg=avg(reg_gamma({timescale|20:70}))*2.12

   print(reg_gamma_avg)
;;;;;;;;;;;;;;;;
   p = rc_beta+1
   reg_beta=rc_beta*2.12
   print(p)
   print(reg_beta)
;;;;;;;;;;;;;;;

   rc_up= regline(reg_dCO2dT,reg_dCO2EdT_up)
   rc_dw= regline(reg_dCO2dT,reg_dCO2EdT_dw)
   p_unc_up=rc_up-p
   p_unc_dw=rc_dw-p
   print(p_unc_up)
   print(p_unc_dw)
   rc_beta_up = rc_up-1
   rc_beta_dw = rc_dw-1
   reg_beta_up=rc_beta_up*2.12
   reg_beta_dw=rc_beta_dw*2.12
   print(reg_beta_up)
   print(reg_beta_dw)

   reg_gamma_up = reg_dCO2EdT_up - reg_dCO2dT * (1 + rc_beta_up)
   reg_gamma_dw = reg_dCO2EdT_dw - reg_dCO2dT * (1 + rc_beta_dw)

   reg_dCO2EdT_up_70yr = reg_dCO2EdT_up({timescale|70})
   print(reg_dCO2EdT_up_70yr)
   reg_gamma_up_70yr = reg_gamma_up({timescale|70}) *2.12
   print(reg_gamma_up_70yr)

   reg_dCO2EdT_dw_70yr = reg_dCO2EdT_dw({timescale|70})
   print(reg_dCO2EdT_dw_70yr)
   reg_gamma_dw_70yr = reg_gamma_dw({timescale|70}) *2.12
   print(reg_gamma_dw_70yr)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    system("rm -f dCO2dT_dCO2EdT_beta_gamma_using_Fourier_Analysis_1880-2017_T_from_HadCRUT4_V2.nc")
    fout=addfile("dCO2dT_dCO2EdT_beta_gamma_using_Fourier_Analysis_1880-2017_T_from_HadCRUT4_V2.nc","c")
    fout->T=T
    fout->CO2=CO2
    fout->CO2E=CO2E
    fout->CO2E_up=CO2E_up
    fout->CO2E_dw=CO2E_dw
    fout->amp_T=amp_T
    fout->amp_CO2=amp_CO2
    fout->amp_CO2E=amp_CO2E
    fout->amp_CO2E_up=amp_CO2E_up
    fout->amp_CO2E_dw=amp_CO2E_dw
    fout->dCO2dT=reg_dCO2dT
    fout->dCO2EdT=reg_dCO2EdT
    fout->dCO2EdT_up=reg_dCO2EdT_up
    fout->dCO2EdT_dw=reg_dCO2EdT_dw
    fout->alpha=alpha
    fout->alpha_up=alpha_up
    fout->alpha_dw=alpha_dw
    fout->gamma=reg_gamma
    fout->gamma_up=reg_gamma_up
    fout->gamma_dw=reg_gamma_dw
    fout->beta=rc_beta
    fout->rc_eta=rc_eta
    fout->rc_zeta=rc_zeta
    fout->rc_beta=rc_beta
    fout->rc_beta_up=rc_beta_up
    fout->rc_beta_dw=rc_beta_dw
    fout->rc_gamma=rc_gamma

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    wks = gsn_open_wks("pdf","plot_spectrum_and_scatter_of_dCO2dT_dCO2EdT_beta_gamma_using_Fourier_Analysis_1880-2017_T_from_HadCRUT4")

    res                               = True 

    res@gsnMaximize          = False
    res@gsnDraw              = False             ; don't draw yet
    res@gsnFrame             = False             ; don't advance frame yet

    res@vpKeepAspect       = True
    res@vpHeightF            = 0.5               ; change aspect ratio of plot
    res@vpWidthF             = 0.5

    res@tiXAxisFontHeightF = 0.022
    res@tiYAxisFontHeightF = 0.022

    res@trXLog               = True
    res@trXMinF              = 1.
    res@trXMaxF              = 150.

    res@tiXAxisString          = "timescale"
    res@tiYAxisString          = "Amplitude of human-emitted CO~B1~2  (ppm)"
    plot_CO2E = gsn_csm_xy (wks, 1/frq(::-1),amp_CO2E(::-1), res)
    res@tiYAxisString          = "Amplitude of atmospheric CO~B1~2  (ppm)"
    plot_CO2 = gsn_csm_xy (wks, 1/frq(::-1),amp_CO2(::-1), res)
    res@tiYAxisString          = "Amplitude of temperature anomaly (~S1~o C)"
    plot_T = gsn_csm_xy (wks, 1/frq(::-1),amp_T(::-1), res)
    res@tiYAxisString          = "Atmospheric ~F33~D~F21~CO~B1~2  /~F33~D~F21~T (ppm ~S1~o C~S1~- ~S1~1 )"
    plot_eta = gsn_csm_xy (wks, 1/frq(::-1),reg_dCO2dT(::-1), res)
    res@tiYAxisString          = "Anthropogenic ~F33~D~F21~CO~B1~2  /~F33~D~F21~T (ppm ~S1~o C~S1~- ~S1~1 )"
    plot_zeta = gsn_csm_xy (wks, 1/frq(::-1),reg_dCO2EdT(::-1), res)

;;;;;;;;;;;;;;
    res_sc =True

    res_sc@gsnDraw            = False
    res_sc@gsnFrame           = False

    res_sc@gsnMaximize        = False
    res_sc@vpKeepAspect       = True
    res_sc@vpHeightF          = 0.5               ; change aspect ratio of plot
    res_sc@vpWidthF           = 0.5

    res_sc@trXMinF            = 0            ; set x-axis minimum
    res_sc@trXMaxF            = 300
    res_sc@trYMinF            = 0
    res_sc@trYMaxF            = 600
; res_sc@tmXBLabelFontHeightF = 0.03
; res_sc@tmYLLabelFontHeightF = 0.03
    res_sc@pmLegendDisplayMode  = "Never"             ; turn on legend
    res_sc@lgPerimOn = False
    res_sc@tiMainFontHeightF  = 0.024                    ; change font heights
    res_sc@tiXAxisFontHeightF = 0.022
    res_sc@tiYAxisFontHeightF = 0.022
    res_sc@xyMarkLineModes     = (/"Markers","Lines"/)  ; choose which have markers
    res_sc@xyMarkers           = 16                     ; choose type of marker
    res_sc@xyMarkerColor       = "red"                  ; Marker color
    res_sc@xyMarkerSizeF       = 0.008                  ; Marker size (default 0.01)
    res_sc@xyDashPatterns      = 1                      ; solid line
    res_sc@xyLineThicknesses   = (/2,3/)                ; set second line to 2
; res_sc@tiYAxisString          = ""
    res_sc@tiYAxisString          = "Anthropogenic ~F33~D~F21~CO~B1~2  /~F33~D~F21~T (ppm ~S1~o C~S1~- ~S1~1 )"
    res_sc@tiXAxisString          = "Atmospheric ~F33~D~F21~CO~B1~2  /~F33~D~F21~T (ppm ~S1~o C~S1~- ~S1~1 )"
    res_sc@tiMainString           = " "

  plot_scatter  = gsn_csm_xy (wks, reg_dCO2dT,data_scatter,res_sc)

    ;************************************************
    resP  = True
    resP@gsnFrame          = False           ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod

    resP@amJust            = "TopLeft"
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.025
    resP@gsnPanelFigureStrings = (/" ~F22~a"," ~F22~b"," ~F22~c"," ~F22~d"," ~F22~e"," ~F22~f"/)

;    resP@gsnPanelCenter        = False
;    resP@txPosXF          = 0.52
;    resP@txPosYF          = 0.71
;    resP@txFontHeightF    = 0.017

    gsn_panel(wks,(/plot_CO2E,plot_CO2,plot_T,plot_zeta,plot_eta,plot_scatter/),(/2,3/),resP)
    frame(wks)
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end

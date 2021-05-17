;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    year=ispan(1850,2017,1)
    nyear=dimsizes(year)

    fpath_CO2="./"
    file_CO2="CO2_global_mean_1850-2017.txt"
    data_CO2 = asciiread(fpath_CO2+file_CO2,(/nyear/),"float")
    data_CO2!0="time"
    data_CO2&time=ispan(1850,2017,1)

    fpath_CO2E="./"
    file_CO2E="CO2E_global_annual_cumulated_FF_plus_LUC_1850-2017.txt"
    data_CO2E = asciiread(fpath_CO2E+file_CO2E,(/nyear/),"float")
    data_CO2E!0="time"
    data_CO2E&time=ispan(1850,2017,1)

    fpath_T="./"
    file_T="T_global_mean_1850-2017_from_HadCRUT4.nc"
    fin=addfile(fpath_T+file_T,"r")
    data_T=fin->T({time|1850:2017})

    T_dtrend = dtrend(data_T,False)

    asciiwrite("CO2_emissions_global_annual_cumulated_FF_plus_LUC_1850-2017.txt",data_CO2E)
    asciiwrite("T_detrended_global_mean_1850-2017_from_HadCRUT4.txt",T_dtrend)

    C0=data_CO2(0)
    T0=data_T(0)
    beta= 3.22/2.12   ; ppm/ppm
    gam= -10.9/2.12     ; ppm/K

;;;;;;;

    C_nat=new(nyear,float)
    T_nat=new(nyear,float)

    C_hum=new(nyear,float)
    T_hum=new(nyear,float)

    C_nat_up=new(nyear,float)
    T_nat_up=new(nyear,float)

    C_hum_up=new(nyear,float)
    T_hum_up=new(nyear,float)

    C_nat_dw=new(nyear,float)
    T_nat_dw=new(nyear,float)

    C_hum_dw=new(nyear,float)
    T_hum_dw=new(nyear,float)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;for s= 3. K
 
    C=new(nyear,float)
    T=new(nyear,float)
    E=new(nyear,float)

    C(0) = C0
    T(0) = T0

    C=285.
    T=0.
    s=3. ;;K
  ;;; for E=0.
    E=0.
  
   do iyr=1,nyear-1

    T(iyr) = T0 + s/log(2)*log(C(iyr-1)/C0) + T_dtrend(iyr) 

    C(iyr) = C0 + (E(iyr) - gam * (T(iyr)-T0))/(1 + beta)

    T(iyr) = T0 + s/log(2)*log(C(iyr)/C0) + T_dtrend(iyr)

   end do

  ;  print(C)
  ;  print(T)

    C_nat=C
    T_nat=T

  ;;; for E=data_CO2E
    E=data_CO2E

   do iyr=1,nyear-1

    T(iyr) = T0 + s/log(2)*log(C(iyr-1)/C0) + T_dtrend(iyr)

    C(iyr) = C0 + (E(iyr) - gam * (T(iyr)-T0))/(1 + beta)

    T(iyr) = T0 + s/log(2)*log(C(iyr)/C0) + T_dtrend(iyr)

   end do

;    print(C)
;    print(T)

    C_hum=C
    T_hum=T

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;for s= 4.5

    C=285.
    T=0.
    s=4.5 ;;K
  ;;; for E=0.
    E=0.

   do iyr=1,nyear-1

    T(iyr) = T0 + s/log(2)*log(C(iyr-1)/C0) + T_dtrend(iyr)

    C(iyr) = C0 + (E(iyr) - gam * (T(iyr)-T0))/(1 + beta)

    T(iyr) = T0 + s/log(2)*log(C(iyr)/C0) + T_dtrend(iyr)

   end do

;    print(C)
;    print(T)

    C_nat_up=C
    T_nat_up=T

  ;;; for E=data_CO2E
    E=data_CO2E

   do iyr=1,nyear-1

    T(iyr) = T0 + s/log(2)*log(C(iyr-1)/C0) + T_dtrend(iyr)

    C(iyr) = C0 + (E(iyr) - gam * (T(iyr)-T0))/(1 + beta)

    T(iyr) = T0 + s/log(2)*log(C(iyr)/C0) + T_dtrend(iyr)

   end do

;    print(C)
;    print(T)

    C_hum_up=C
    T_hum_up=T

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;for s= 1.5

    C=285.
    T=0.
    s=1.5 ;;K
  ;;; for E=0.
    E=0.

   do iyr=1,nyear-1

    T(iyr) = T0 + s/log(2)*log(C(iyr-1)/C0) + T_dtrend(iyr)

    C(iyr) = C0 + (E(iyr) - gam * (T(iyr)-T0))/(1 + beta)

    T(iyr) = T0 + s/log(2)*log(C(iyr)/C0) + T_dtrend(iyr)

   end do

;    print(C)
;    print(T)

    C_nat_dw=C
    T_nat_dw=T

  ;;; for E=data_CO2E
    E=data_CO2E

   do iyr=1,nyear-1

    T(iyr) = T0 + s/log(2)*log(C(iyr-1)/C0) + T_dtrend(iyr)

    C(iyr) = C0 + (E(iyr) - gam * (T(iyr)-T0))/(1 + beta)

    T(iyr) = T0 + s/log(2)*log(C(iyr)/C0) + T_dtrend(iyr)

   end do

;    print(C)
;    print(T)

    C_hum_dw=C
    T_hum_dw=T


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    asciiwrite("Predicted_annual_CO2_over_1850-2017_using_CO2_emissions_by_box_model.txt",C_hum)
    asciiwrite("Predicted_annual_T_over_1850-2017_using_CO2_emissions_by_box_model.txt",T_hum)
    asciiwrite("Predicted_annual_CO2_up_over_1850-2017_using_CO2_emissions_by_box_model.txt",C_hum_up)
    asciiwrite("Predicted_annual_CO2_dw_over_1850-2017_using_CO2_emissions_by_box_model.txt",C_hum_dw)
    asciiwrite("Predicted_annual_T_up_over_1850-2017_using_CO2_emissions_by_box_model.txt",T_hum_up)
    asciiwrite("Predicted_annual_T_dw_over_1850-2017_using_CO2_emissions_by_box_model.txt",T_hum_dw)

    system("rm -rf Annual_CO2_and_T_with_s_beta_and_gamma_CO2E_using_box_model_from_1850-2017.nc")
    fout=addfile("Annual_CO2_and_T_with_s_beta_and_gamma_CO2E_using_box_model_from_1850-2017.nc","c")
    fout->T_residul=T_dtrend
    fout->CO2E=data_CO2E
    fout->CO2_nat=C_nat
    fout->T_nat=T_nat
    fout->CO2_hum=C_hum
    fout->T_hum=T_hum
    fout->CO2_nat=C_nat_up
    fout->T_nat=T_nat_up
    fout->CO2_hum=C_hum_up
    fout->T_hum=T_hum_up
    fout->CO2_nat_dw=C_nat_dw
    fout->T_nat_dw=T_nat_dw
    fout->CO2_hum_dw=C_hum_dw
    fout->T_hum_dw=T_hum_dw

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    NSE_C_hum_obs=1-sum((C_hum-data_CO2)^2)/sum((data_CO2-avg(data_CO2))^2)
    NSE_T_hum_obs=1-sum((T_hum-data_T)^2)/sum((data_T-avg(data_T))^2)
    print(NSE_C_hum_obs)
    print(NSE_T_hum_obs) 

    r_C_hum_obs=escorc(C_hum,data_CO2)
    r_T_hum_obs=escorc(T_hum,data_T)

    R2_C_hum_obs=r_C_hum_obs^2
    R2_T_hum_obs=r_T_hum_obs^2

     print(R2_C_hum_obs)
     print(R2_T_hum_obs)

     nobs=nyear

     RMSE_C_hum_obs= sqrt(sum((C_hum-data_CO2)^2)/nobs)
     RMSE_T_hum_obs= sqrt(sum((T_hum-data_T)^2)/nobs)

    print(RMSE_C_hum_obs)
    print(RMSE_T_hum_obs)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;;plot
    wks = gsn_open_wks("pdf","plot_timeseries_of_CO2E_CO2_T_using_box_model_over_1850-2017_V1")

    ;************************************************
    ;;;plot CO2E,CO2 over 1850-2017
    res                      = True
    res@gsnMaximize          = False
    res@gsnDraw              = False             ; don't draw yet
    res@gsnFrame             = False             ; don't advance frame yet

    res@vpHeightF            = 0.4               ; change aspect ratio of plot
    res@vpWidthF             = 0.6

  ;  res@vpXF                 = 0.5                   ; width of legend (NDC)
  ;  res@vpYF                 = 0.8

    res@tiMainFontHeightF    = 0.015
    res@tiXAxisFontHeightF   = 0.015
    res@tiYAxisFontHeightF   = 0.015
    res@tmXBLabelFontHeightF = 0.015
    res@tmYLLabelFontHeightF = 0.015

    res@pmLegendDisplayMode  = "Always"             ; turn on legend
    res@pmLegendSide           = "Top"                ; Change location of
    res@lgPerimOn = False
    res@pmLegendParallelPosF   = 0.4                  ; move units right
    res@pmLegendOrthogonalPosF = -0.3                ; move units down
    res@pmLegendWidthF         = 0.08                  ; Change width and
    res@pmLegendHeightF        = 0.075                   ; height of legend.
    res@lgLabelFontHeightF     = .016                   ; change font height

    res@tmXBMode             = "Manual"
    res@tmXBOn               = True
    res@tmXBMinorOn          = True
    res@tmXTOn               = False
    res@tmXTMinorOn          = False
    res@tmXBTickStartF       = 1850
    res@tmXBTickEndF         = 2020
    res@tmXBTickSpacingF     = 50
    res@tmXBMajorOutwardLengthF = 0.0
    res@tmXBMinorOutwardLengthF = 0.0
    res@tmXBMinorPerMajor    = 4

    res@tiXAxisString        = "Year"  ; xaxis
    res@tiYAxisString        = "Cumulative CO~B1~2  emissions (GtC)"
    res@trYMinF              = -100.
    res@trYMaxF              = 800.
    res@trXMinF              = 1850
    res@trXMaxF              = 2020
    res@tmYLMajorOutwardLengthF = 0.00
    res@tmYLMinorOutwardLengthF = 0.00
    res@tmYLMinorPerMajor    = 1
    res@tmYRMinorPerMajor    = 1

    res@xyDashPatterns       = (/0,1/)
    res@xyLineThicknessF     = (/2.5,2./)
    res@xyLineColors         = (/"black","black"/)

    res@tiYAxisString        = "Cumulative CO~B1~2  emissions (GtC)"

    res@xyExplicitLegendLabels  = (/"no CO~B1~2  emissions",""/)

    E=0.
    plot_E = gsn_csm_xy (wks, year,(/E*2.12/),res)

    res@xyLineColors         = (/"red","black"/)
    res@xyExplicitLegendLabels  = (/"CO~B1~2  emissions",""/)
    plot_CO2E = gsn_csm_xy (wks, year,(/data_CO2E*2.12/),res)

    res@tiYAxisString        = "CO~B1~2  concentration (ppm)"
    res@trYMinF              = 280.
    res@trYMaxF              = 420.
    res@xyLineColors         = (/"black","black"/)
    res@xyExplicitLegendLabels  = (/"no CO~B1~2  emissions",""/)
    plot_Cnat_top = gsn_csm_xy (wks, year,(/C_nat/),res)

    res@xyLineColors         = (/"red","black"/)
    res@xyExplicitLegendLabels  = (/"Prediction using CO~B1~2  emissions","Observation"/)
    plot_Chum_top = gsn_csm_xy (wks, year,(/C_hum,data_CO2/),res)

    res@tiYAxisString        = "Temperature (K)"
    res@trYMinF              = -1.
    res@trYMaxF              = 2.
    res@xyLineColors         = (/"black","black"/)
    res@xyExplicitLegendLabels  = (/"no CO~B1~2  emissions",""/)
    plot_Tnat_top = gsn_csm_xy (wks, year,(/T_nat/),res)

    delete(res@xyLineColors)
    delete(res@xyExplicitLegendLabels)
    res@xyLineColors         = (/"red","black"/)
    res@xyExplicitLegendLabels  = (/"Prediction using CO~B1~2  emissions","Observation"/)
    plot_Thum_top = gsn_csm_xy (wks, year,(/T_hum,data_T/),res)

   ; Create a plot with the area between both curves filled in pink.
    res@tiYAxisString        = "CO~B1~2  concentration (ppm)"
    res@trYMinF              = 280.
    res@trYMaxF              = 420.
    res@pmLegendDisplayMode  = "Never"
    delete(res@xyLineColors)
    res@gsnXYFillColors = "LightGray"
    res@xyLineColor     = -1      ; We don't want the line, so make it transparent.
    plot_Cnat  = gsn_csm_xy (wks,year,(/C_nat_dw,C_nat_up/),res)  ; Create another filled XY plot.
    overlay(plot_Cnat,plot_Cnat_top)

    plot_Chum  = gsn_csm_xy (wks,year,(/C_hum_dw,C_hum_up/),res)  ; Create another filled XY plot.
    overlay(plot_Chum,plot_Chum_top)

   ; Create a plot with the area between both curves filled in pink.
    res@tiYAxisString        = "Temperature (K)"
    res@trYMinF              = -1.
    res@trYMaxF              = 2.
    res@pmLegendDisplayMode  = "Never"
;    delete(res@xyLineColors)
    res@gsnXYFillColors = "LightPink"
    res@xyLineColor     = -1      ; We don't want the line, so make it transparent.
    plot_Tnat  = gsn_csm_xy (wks,year,(/T_nat_dw,T_nat_up/),res)  ; Create another filled XY plot.
    overlay(plot_Tnat,plot_Tnat_top)

    plot_Thum  = gsn_csm_xy (wks,year,(/T_hum_dw,T_hum_up/),res)  ; Create another filled XY plot.
    overlay(plot_Thum,plot_Thum_top)
    delete(res@gsnXYFillColors)
    ;************************************************
    resP  = True
    resP@gsnFrame          = False           ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod

    resP@amJust            = "TopLeft"
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.02
    resP@gsnPanelFigureStrings = (/" ~F22~a"," ~F22~b"," ~F22~c"," ~F22~d"," ~F22~e"," ~F22~f"/)
    resP@gsnPanelYWhiteSpacePercent = 1
    resP@gsnPanelXWhiteSpacePercent = 3

    gsn_panel(wks,(/plot_E,plot_CO2E,plot_Cnat,plot_Chum,plot_Tnat,plot_Thum/),(/3,2/),resP)


end

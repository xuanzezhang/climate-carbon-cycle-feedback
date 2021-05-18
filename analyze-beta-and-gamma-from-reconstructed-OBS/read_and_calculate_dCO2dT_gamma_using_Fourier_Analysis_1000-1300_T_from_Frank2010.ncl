;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   styr=1000
   edyr=1300
 ;;;1000-1850
    lmyear=ispan(styr,edyr,1)
    nlmyear = dimsizes(lmyear)

    ;; CO2 from Low Dome
    lmCO2_LowDome=new(nlmyear,"float")
    lmCO2_LowDome!0="time"
    lmCO2_LowDome&time=lmyear

    fpath_lmCO2_LowDome="./"
    file_lmCO2_LowDome="CO2_global_mean_LowDome_1-2000.txt"
    data_lmCO2_LowDome= asciiread(fpath_lmCO2_LowDome+file_lmCO2_LowDome,(/2000/),"float")
    data_lmCO2_LowDome!0="time"
    data_lmCO2_LowDome&time=ispan(1,2000,1)
    lmCO2_LowDome({time|styr:edyr}) = (/data_lmCO2_LowDome({time|styr:edyr})/)

    ;; CO2 from WAIS
    lmCO2_WAIS=new(nlmyear,"float")
    lmCO2_WAIS!0="time"
    lmCO2_WAIS&time=lmyear

    fpath_lmCO2_WAIS="./"
    file_lmCO2_WAIS="CO2_global_mean_WAIS_Bauska2015_740-1940.txt"
    data_lmCO2_WAIS= asciiread(fpath_lmCO2_WAIS+file_lmCO2_WAIS,(/1201/),"float")
    data_lmCO2_WAIS!0="time"
    data_lmCO2_WAIS&time=ispan(740,1940,1)
    lmCO2_WAIS({time|styr:edyr}) = (/data_lmCO2_WAIS({time|styr:edyr})/)

    ;; CO2 from DML
    lmCO2_DML=new(nlmyear,"float")
    lmCO2_DML!0="time"
    lmCO2_DML&time=lmyear

    fpath_lmCO2_DML="./"
    file_lmCO2_DML="CO2_global_mean_DML_Siegenthaler2005_950-1880.txt"
    data_lmCO2_DML= asciiread(fpath_lmCO2_DML+file_lmCO2_DML,(/931/),"float")
    data_lmCO2_DML!0="time"
    data_lmCO2_DML&time=ispan(950,1880,1)
    lmCO2_DML({time|styr:edyr}) = (/data_lmCO2_DML({time|styr:edyr})/)

    ;;T from Frank2010
    lmCO2_Frank=new(nlmyear,"float")
    lmCO2_Frank!0="time"
    lmCO2_Frank&time=lmyear

    lmT_Frank=new(nlmyear,"float")
    lmT_Frank!0="time"
    lmT_Frank&time=lmyear

    fpath_lmCO2_Frank="./"
    file_lmCO2_Frank="CO2_global_mean_Frank2010_1000-2000.txt"
    data_lmCO2_Frank = asciiread(fpath_lmCO2_Frank+file_lmCO2_Frank,(/1001/),"float")
    lmCO2_Frank({time|styr:edyr}) = (/data_lmCO2_Frank(styr-1000:edyr-1000)/)

    fpath_lmT_Frank="./"
    file_lmT_Frank="T_10yr_smth_Frank2010_data.nc"
    fin=addfile(fpath_lmT_Frank+file_lmT_Frank,"r")
    data_lmT_Frank=fin->T ;; (year,nrun)

    nrun=dimsizes(data_lmT_Frank(0,:))
    data_lmT_Frank!0="time"
    data_lmT_Frank!1="MCrun"
    data_lmT_Frank&MCrun=ispan(1,nrun,1)

    data_lmT_Frank_avg = dim_avg_n(data_lmT_Frank,1)
    nyear_Frank = dimsizes(data_lmT_Frank_avg)
    data_lmT_Frank_avg!0="time"
    data_lmT_Frank_avg&time=ispan(1000,(nyear_Frank+1000-1),1)
    ;print(data_lmT_Frank_avg&time)

    lmT_Frank({time|styr:edyr}) = (/data_lmT_Frank_avg({time|styr:edyr})/)
    lmT_Frank = lmT_Frank - dim_avg_n(data_lmT_Frank_avg({time|1961:1990}),0)

    lmT_Frank_run = data_lmT_Frank({time|styr:edyr},{MCrun|:})

    lmT_Frank_run_30yr_avg = dim_avg_n(data_lmT_Frank({time|1961:1990},{MCrun|:}),0)
    do irun=0,nrun-1
       lmT_Frank_run(:,irun) = lmT_Frank_run(:,irun) - lmT_Frank_run_30yr_avg(irun)
    end do

;************************************************
; Calculate harmonics           
;   finfo(0,:,:,nhar) - amplitude of each harmonic
;   finfo(1,:,:,nhar) - phase of each harmonic
;   finfo(2,:,:,nhar) - percent variance explained by each harmonic 
;************************************************
    nfrq=nlmyear/2
    K=ispan(1,nfrq,1)*1.0
    frq=K/nlmyear
;   print(frq)
    timescale = 1.0/frq

    nlmCO2=3

    lmCO2 = new((/nlmCO2,nlmyear/),"float")
    lmCO2!0="lmCO2"
    lmCO2&lmCO2=ispan(1,3,1)
    lmCO2@lmCO2="1.LowDome,2.WAIS,3.DML"
    lmCO2(0,:) = lmCO2_LowDome(:)
    lmCO2(1,:) = lmCO2_WAIS(:) 
    lmCO2(2,:) = lmCO2_DML(:)

    amp_lmCO2 = new((/nlmCO2,nfrq/),"float")

    do iCO2=0,nlmCO2-1
       finfo_lmCO2  = fourier_info(lmCO2(iCO2,:), nfrq, 1.)
       amp_lmCO2(iCO2,:) = finfo_lmCO2(0,:)
    end do

    amp_lmCO2!0="lmCO2"
    amp_lmCO2&lmCO2=ispan(1,3,1)
    amp_lmCO2!1="timescale"
    amp_lmCO2&timescale=timescale
    amp_lmCO2@lmCO2="1.LowDome,2.WAIS,3.DML"


    amp_lmT = new((/nfrq/),"float")
    finfo_lmT    = fourier_info(lmT_Frank, nfrq, 1.) 
    amp_lmT = finfo_lmT(0,:)
    amp_lmT!0="timescale"
    amp_lmT&timescale=timescale


    reg_dCO2dT=new((/nlmCO2,nfrq/),"float")
    reg_dCO2dT!0="lmCO2"
    reg_dCO2dT&lmCO2=ispan(1,3,1)
    reg_dCO2dT!1="timescale"
    reg_dCO2dT&timescale=timescale
    reg_dCO2dT@lmCO2="1.LowDome,2.WAIS,3.DML"

    do iCO2=0,nlmCO2-1
       reg_dCO2dT(iCO2,:) = (amp_lmCO2(iCO2,:)/amp_lmT(:))
    end do

    rc = 2.5177
    rc_sd = 0.0131
    rc_beta = rc - 1

    rc_gamma =new((/nlmCO2/),"float")
    rc_gamma!0="lmCO2"
    rc_gamma&lmCO2=ispan(1,3,1)
    rc_gamma@lmCO2="1.LowDome,2.WAIS,3.DML"
    rc_gamma_sd =new((/nlmCO2/),"float")
    rc_gamma_sd!0="lmCO2"
    rc_gamma_sd&lmCO2=ispan(1,3,1)
    rc_gamma_sd@lmCO2="1.LowDome,2.WAIS,3.DML"

    reg_dCO2dT_100yr = reg_dCO2dT({lmCO2|:},{timescale|100})
    rc_gamma(:)= -1.*reg_dCO2dT_100yr(:)*rc*2.12   ;;; convert to GtC/K
    rc_gamma_sd(:)= -1.*reg_dCO2dT_100yr(:)*rc_sd*2.12
    print(rc_gamma)
    print(rc_gamma_sd)
    print(reg_dCO2dT_100yr)

    lmT_Frank_avg=avg(lmT_Frank)
    print(lmT_Frank_avg)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;estimate gamma from ensemble runs
    amp_lmT_run = new((/nrun,nfrq/),"float")

    num_missing=1  
    do irun=0,nrun-1
       ; print(irun)
       ;print(lmT_Frank_run(:,irun))
       if (all(ismissing(lmT_Frank_run(:,irun)))) then
          num_missing=num_missing+1
       else
          finfo_lmT_run = fourier_info(lmT_Frank_run(:,irun), nfrq, 1.)
         ; print(spec_lmT_run)
          amp_lmT_run_in      = finfo_lmT_run(0,:)
          amp_lmT_run(irun,:) = runave(amp_lmT_run_in,3,1)
       ;frq      = spec_lmT_run@frq
       end if
    end do
    print(num_missing)

    reg_dCO2dT_run=new((/nlmCO2,nrun,nfrq/),"float")

    reg_dCO2dT_run!0="lmCO2"
    reg_dCO2dT_run&lmCO2=ispan(1,3,1)
    reg_dCO2dT_run!1="MCrun"
    reg_dCO2dT_run&MCrun=ispan(1,nrun,1)
    reg_dCO2dT_run!2="timescale"
    reg_dCO2dT_run&timescale=timescale
    reg_dCO2dT_run@lmCO2="1.LowDome,2.WAIS,3.DML"

    do iCO2=0,nlmCO2-1
       do irun=0,nrun-1
          reg_dCO2dT_run(iCO2,irun,:) = (amp_lmCO2(iCO2,:)/amp_lmT_run(irun,:))
       end do
    end do

    rc_gamma_run =new((/nlmCO2,nrun/),"float")
    rc_gamma_run!0="lmCO2"
    rc_gamma_run&lmCO2=ispan(1,3,1)
    rc_gamma_run!1="MCrun"
    rc_gamma_run&MCrun=ispan(1,nrun,1)
    rc_gamma_run@lmCO2="1.LowDome,2.WAIS,3.DML"
    rc_gamma_sd_run =new((/nlmCO2,nrun/),"float")
    rc_gamma_sd_run!0="lmCO2"
    rc_gamma_sd_run&lmCO2=ispan(1,3,1)
    rc_gamma_sd_run!1="MCrun"
    rc_gamma_sd_run&MCrun=ispan(1,nrun,1)
    rc_gamma_sd_run@lmCO2="1.LowDome,2.WAIS,3.DML"


    reg_dCO2dT_run_100yr = reg_dCO2dT_run({lmCO2|:},{MCrun|:},{timescale|100})
    rc_gamma_run(:,:)= -1.*reg_dCO2dT_run_100yr(:,:)*rc*2.12   ;;; convert to GtC/K
    rc_gamma_sd_run(:,:)= -1.*reg_dCO2dT_run_100yr(:,:)*rc_sd*2.12

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    system("rm -f dCO2dT_gamma_using_Fourier_Analysis_"+styr+"-"+edyr+"_T_from_Frank2010.nc")
    fout=addfile("dCO2dT_gamma_using_Fourier_Analysis_"+styr+"-"+edyr+"_T_from_Frank2010.nc","c")
    fout->T=lmT_Frank
    fout->CO2=lmCO2
    fout->amp_T=amp_lmT
    fout->amp_CO2=amp_lmCO2
    fout->dCO2dT=reg_dCO2dT
    fout->beta=rc_beta
    fout->beta_sd=rc_sd
    fout->gamma=rc_gamma
    fout->gamma_sd=rc_gamma_sd
    fout->T_run=lmT_Frank_run
    fout->amp_T_run=amp_lmT_run
    fout->dCO2dT_run=reg_dCO2dT_run
    fout->gamma_run=rc_gamma_run
    fout->gamma_sd_run=rc_gamma_sd_run


 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end



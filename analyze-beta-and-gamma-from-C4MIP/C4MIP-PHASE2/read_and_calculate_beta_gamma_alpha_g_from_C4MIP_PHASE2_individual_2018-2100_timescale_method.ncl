;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=2018
Endyr=2100
nyear=Endyr-Startyr+1
year=ispan(Startyr,Endyr,1)

Models=(/"BERN","CCSM1","CLIMBER","FRCGC","HADLEY","IPSL","LLNL","LOOP","MPI","UMD","UVIC"/)
nmodel=dimsizes(Models)

Exp="global"
varname=(/"UNC","COU"/)
nvar=dimsizes(varname)

nfrq=nyear
nmax=nyear

BETA=new((/nfrq,nmax/),"float")
BETA!0="timescale"
BETA!1="number"
BETA&timescale=ispan(1,nyear,1)
BETA&number=ispan(1,nmax,1)
BETA@units="GtC/ppm"

GAMMA=new((/nfrq,nmax/),"float")
GAMMA!0="timescale"
GAMMA!1="number"
GAMMA&timescale=ispan(1,nyear,1)
GAMMA&number=ispan(1,nmax,1)
GAMMA@units="GtC/K"

ALPHA=new((/nfrq,nmax/),"float")
ALPHA!0="timescale"
ALPHA!1="number"
ALPHA&timescale=ispan(1,nyear,1)
ALPHA&number=ispan(1,nmax,1)

G=new((/nfrq,nmax/),"float")
G!0="timescale"
G!1="number"
G&timescale=ispan(1,nyear,1)
G&number=ispan(1,nmax,1)

one=new((/nfrq,nmax/),"float")

CO2E_FLUX  =new((/nyear/),"float")
CO2E       =new((/nyear/),"float")
CO2        =new((/nyear/),"float")
ST         =new((/nyear/),"float")
LA_CO2_FLUX=new((/nyear/),"float")
OA_CO2_FLUX=new((/nyear/),"float")

dTA=new(nyear,"float")
dTA!0="year"
dTA&year=year

dCA=new(nyear,"float")
dCA!0="year"
dCA&year=year

dCL=new(nyear,"float")
dCL!0="year"
dCL&year=year

dCO=new(nyear,"float")
dCO!0="year"
dCO&year=year

CA=new(nyear,"float")
CA!0="year"
CA&year=year

CS=new(nyear,"float")
CS!0="year"
CS&year=year

TA=new(nyear,"float")
TA!0="year"
TA&year=year

BETA=BETA@_FillValue
GAMMA=GAMMA@_FillValue

dTA=0.
dCA=0.
dCL=0.
dCO=0.

dTA_in=dTA
dCA_in=dCA
dCL_in=dCL
dCO_in=dCO



do imodel=0, nmodel-1
   Model=Models(imodel)

 do i= 0, nvar-1
   print(varname(i))

   OutFileNameAnn="C4MIP_"+Models(imodel)+"_BETA_GAMMA_ALPHA_and_G_over_timescale_2018-2100.nc"

   FilePath="./"
   FileNameAnn=systemfunc("ls C4MIP_"+Models(imodel)+"_"+Exp+"_from_*_"+varname(i)+".nc")
   print(FilePath+FileNameAnn)
   fin=addfile(FilePath+FileNameAnn,"r")
   CO2E_in=fin->CO2E
   CO2_in=fin->CO2
   ST_in=fin->ST
   LA_CO2_FLUX_in=fin->LA_CO2_FLUX
   OA_CO2_FLUX_in=fin->OA_CO2_FLUX

   ntime=dimsizes(CO2E_in)
   styr=ntime-(2100-Endyr)-nyear
   edyr=ntime-(2100-Endyr)-1

   CO2(:)        =(/CO2_in(styr:edyr)/)  ;; ppm
   ST(:)         =(/ST_in(styr:edyr)/)  
   LA_CO2_FLUX(:)=(/LA_CO2_FLUX_in(styr:edyr)/)
   OA_CO2_FLUX(:)=(/OA_CO2_FLUX_in(styr:edyr)/)

   do iyr=1,nyear-1
      dTA_in(iyr)=ST(iyr)-ST(iyr-1)
      dCA_in(iyr)=CO2(iyr)-CO2(iyr-1)
      dCL_in(iyr)=LA_CO2_FLUX(iyr)
      dCO_in(iyr)=OA_CO2_FLUX(iyr)
   end do

   dTA= runave(dTA_in,3,0)
   dCA= runave(dCA_in,3,0)
   dCL= runave(dCL_in,3,0)
   dCO= runave(dCO_in,3,0)

     ;; calculate ALPHA from fully-coupled simulations
     if (varname(i).eq."COU") then
      do iyr=1,nyear-1
        CA(iyr)=sum(dCA(0:iyr))*2.12 ;; ppm to GtC
        CS(iyr)=sum(dCL(0:iyr))+sum(dCO(0:iyr)) ;; GtC
      end do

      do ifrq=1,nfrq-1
        n=0
        do iyr=1,nyear-1-ifrq
          a = CS(iyr+ifrq)-CS(iyr)
          b = CA(iyr+ifrq)-CA(iyr)
          c = a + b
          if (b.ne.0) then
          ALPHA(ifrq,iyr) = b/c
          n = n +1
          end if
        end do
;        print("ifrq = "+ifrq)
;        print(n)
      end do
     end if

     ;; calculate BETA from biogeochemically-coupled simulations
     if (varname(i).eq."UNC") then
      do iyr=1,nyear-1
        CA(iyr)=sum(dCA(0:iyr)) ;; ppm
        CS(iyr)=sum(dCL(0:iyr))+sum(dCO(0:iyr)) ;; GtC
      end do
      print("CA for 2100 in COU is")
      print(CA(nyear-1))

      do ifrq=1,nfrq-1
        n=0
        do iyr=1,nyear-1-ifrq
          a = CS(iyr+ifrq)-CS(iyr)
          b = CA(iyr+ifrq)-CA(iyr)
          if (b.ne.0) then
          BETA(ifrq,iyr) = a/b
          n = n +1
          end if
        end do
 ;       print("ifrq = "+ifrq)
 ;       print(n)
      end do
     end if
 
     ;; calculate GAMMA from fully-coupled simulations
     if (varname(i).eq."COU") then
      do iyr=1,nyear-1
        TA(iyr)=sum(dTA(0:iyr))
        CA(iyr)=sum(dCA(0:iyr));; GtC
        CS(iyr)=sum(dCL(0:iyr))+sum(dCO(0:iyr))
      end do
      print("TA for 2100 in COU is")
      print(TA(nyear-1))

      do ifrq=1,nfrq-1
        n=0
        do iyr=1,nyear-1-ifrq
          t = TA(iyr+ifrq)-TA(iyr)
          a = CS(iyr+ifrq)-CS(iyr)
          b = CA(iyr+ifrq)-CA(iyr)
          if (t.ne.0) then
          GAMMA(ifrq,iyr)=(a-BETA(ifrq,iyr)*b)/t
          n = n +1
;          print(a)
;          print(b)
;          print(t)
;          print(GAMMA(ifrq,iyr))
          end if
        end do
;        print("ifrq = "+ifrq)
;        print(n)
      end do
     end if
 
   delete(CO2E_in)
   delete(CO2_in)
   delete(ST_in)
   delete(LA_CO2_FLUX_in)
   delete(OA_CO2_FLUX_in)

 end do ; :

   ALPHA=where(abs(ALPHA).ge.1.,ALPHA@_FillValue,ALPHA)
   BETA=where(abs(BETA).ge.10,BETA@_FillValue,BETA)
   GAMMA=where(abs(GAMMA).ge.200,GAMMA@_FillValue,GAMMA)

   one = 1.
   printVarSummary(one)
   C = ALPHA * (one + BETA/2.12)
   G = (C - one)/C

   BETA_avg = dim_avg_n(BETA,1)
   BETA_std = dim_stddev_n(BETA,1)
   GAMMA_avg = dim_avg_n(GAMMA,1)
   GAMMA_std = dim_stddev_n(GAMMA,1)

   ALPHA_avg = dim_avg_n(ALPHA,1)
   ALPHA_std = dim_stddev_n(ALPHA,1)
   one_avg = dim_avg_n(one,1)
   C_avg = ALPHA_avg * (one_avg + BETA_avg/2.12)
   G_avg = (C_avg - one_avg)/C_avg
   ;G_avg = dim_avg_n(G,1)
   G_std = dim_stddev_n(G,1)
   BETA_avg!0="timescale"
   BETA_avg&timescale=ispan(1,nyear,1)
   copy_VarCoords(BETA_avg,BETA_std)
   copy_VarCoords(BETA_avg,GAMMA_avg)
   copy_VarCoords(BETA_avg,GAMMA_std)
   copy_VarCoords(BETA_avg,ALPHA_avg)
   copy_VarCoords(BETA_avg,ALPHA_std)
   copy_VarCoords(BETA_avg,G_avg)
   copy_VarCoords(BETA_avg,G_std)
;   print(G_avg)

   system("rm -rf "+OutFileNameAnn)
   fout=addfile(OutFileNameAnn,"c")
   fout->BETA=BETA
   fout->GAMMA=GAMMA
   fout->ALPHA=ALPHA
   fout->G=G
   fout->BETA_avg=BETA_avg
   fout->BETA_std=BETA_std
   fout->GAMMA_avg=GAMMA_avg
   fout->GAMMA_std=GAMMA_std
   fout->ALPHA_avg=ALPHA_avg
   fout->ALPHA_std=ALPHA_std
   fout->G_avg=G_avg
   fout->G_std=G_std

end do ; i



end

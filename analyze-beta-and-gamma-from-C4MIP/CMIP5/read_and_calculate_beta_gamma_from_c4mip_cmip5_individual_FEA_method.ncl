;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1850
Endyr=1989
nyear=Endyr-Startyr+1
year=ispan(Startyr,Endyr,1)

Models=(/"BCC","CCCma","CCSM","Hadley","IPSL","JAMSTEC","MPI","NCC","UVic"/)
nmodel=dimsizes(Models)

beta_model=new(nmodel,float)
gamma_model=new(nmodel,float)
fbg_model=new(nmodel,float)
gamma_COU_BIO_model=new(nmodel,float)
TA_COU=new(nmodel,float)
CA_COU=new(nmodel,float)
CS_COU=new(nmodel,float)
TA_BIO=new(nmodel,float)
CA_BIO=new(nmodel,float)
CS_BIO=new(nmodel,float)
TA_RAD=new(nmodel,float)
CA_RAD=new(nmodel,float)
CS_RAD=new(nmodel,float)

expname=(/"CONTROL","1pc_BIO","1pc_RAD","1pc_COU"/)
nexp=dimsizes(expname)

do imodel=0, nmodel-1
   Model=Models(imodel)
   OutFileNameAnn=Models(imodel)+"_beta_and_gamma_using_FEA_method.txt"
 
 do i= 0, nexp-1
   ;print(expname(i))

   FilePath="./"
   FileNameAnn=Models(imodel)+"_"+expname(i)+".nc"
   fin=addfile(FilePath+FileNameAnn,"r")
   ST=fin->ST
   CO2=fin->CO2
   LA_CO2_FLUX=fin->LA_CO2_FLUX
   OA_CO2_FLUX=fin->OA_CO2_FLUX

   dTA_in=ST(nyear-1)-ST(0)
   dCA_in=CO2(nyear-1)-CO2(0)   
   dCL_in=sum(LA_CO2_FLUX)
   dCO_in=sum(OA_CO2_FLUX)

;   print(dTA)
;   print(dCA)  

   if (expname(i).eq."1pc_BIO") then
        CA_BIO(imodel)=dCA_in
        CS_BIO(imodel)=dCL_in + dCO_in
        beta_model(imodel)=CS_BIO(imodel)/CA_BIO(imodel)
   end if
   if (expname(i).eq."1pc_RAD") then
        TA_RAD(imodel)=dTA_in
        CS_RAD(imodel)=dCL_in + dCO_in
        gamma_model(imodel)=CS_RAD(imodel)/TA_RAD(imodel)
   end if

   if (expname(i).eq."1pc_COU") then
        TA_COU(imodel)=dTA_in
        CA_COU(imodel)=dCA_in 
        CS_COU(imodel)=dCL_in + dCO_in
   end if

   delete(dTA_in)
   delete(dCA_in)
   delete(dCL_in)
   delete(dCO_in)

 end do ; i

end do ; imodel

gamma_COU_BIO_model=(CS_COU-CS_BIO)/TA_COU

fbg_model=(CS_COU-CS_BIO-CS_RAD)/CA_COU/TA_COU
fbg_dCA_model=fbg_model*CA_COU

beta_model_avg=avg(beta_model)
beta_model_std=stddev(beta_model)

fbg_model=fbg_model*1000.
gamma_model_avg=avg(gamma_model)
gamma_model_std=stddev(gamma_model)

gamma_COU_BIO_model_avg=avg(gamma_COU_BIO_model)
gamma_COU_BIO_model_std=stddev(gamma_COU_BIO_model)

fbg_model_avg=avg(fbg_model)
fbg_model_std=stddev(fbg_model)

fbg_dCA_model_avg=avg(fbg_dCA_model)
fbg_dCA_model_std=stddev(fbg_dCA_model)

print(beta_model)
print(beta_model_avg)
print(beta_model_std)
print(gamma_model)
print(gamma_model_avg)
print(gamma_model_std)
print(gamma_COU_BIO_model)
print(gamma_COU_BIO_model_avg)
print(gamma_COU_BIO_model_std)
print(fbg_model)
print(fbg_model_avg)
print(fbg_model_std)
print(fbg_dCA_model)
print(fbg_dCA_model_avg)
print(fbg_dCA_model_std)



end

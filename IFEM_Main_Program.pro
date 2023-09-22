
;===============================IFEM======================================
;An Independent Framework-based Evapotranspiration Model (IFEM)
;Based on the Fractional Vegetation Cover--Land Surface Temperature space
;
;Function:
;    Program can be used for estimating surface energy components
;    A dual-source scheme for soil Evaporation and canopy Transpiration
;
;Inputs:
;    Program requires for basic scalar ground meteorological data
;    OR spatial raster fields of air temperature and Vapor Pressure Deficit
;
;Update history:
;    2023-Sept.-13 Code simplification
;    2023-July-03  Debug and annotate
;    2023-Apr.-09  Initialization code
;
;Please Cite:  Wang, Shuai, Wang, Chaozi, et al, 2023. DOI: 10.1016/j.rse.2023.113792
;
;E-mail: 214544015@qq.com (Dr. Shuai,Wang)
;        huozl@cau.edu.cn (Prof. Zailin,Huo)
;        
;This program has been tested in ENVI/IDL 5.3/5.6 environment
;=========================================================================



;===========================================================================
;==============================Main Program Part============================
;===========================================================================

;Main Program
PRO IFEM_Main_Program

  ;==============Input Parameters=============
  ;[x] represent essential parameters, others optional parameters (if ignored, set to !NULL)

  ;UTC time of overpass [Required]
  Year = 2017   ;[x]Year
  Month = 7     ;[x]Month
  Day = 4       ;[x]Day
  Time = 3.5    ;[x]Time[HH]

  ;Instantaneous meteorological observations [Required]
  u_ref = 0.5    ;[x]average wind velocity [m s-1]
  z_ref = 2.0    ;[x]height of wind velocity observation [m]
  Ta = 28.7      ;[x]average air temperature [℃] ;Scalar or raster file path
  zT = 2.0       ;[x]height of air temperature observation [m]
  VPD = 2.01     ;[x]average Vapor Pressure Deficit [kPa] ;Scalar or raster file path

  ;Daytime meteorological parameters for Daily ET caculation [Optional]
  Ta_max = 33.76   ;Air temperature maximum [℃]
  Ta_min = 14.51   ;Air temperature minimum [℃]
  RH_max = 93.0    ;Relative humidity maximum [%]
  RH_min = 37.3    ;Relative humidity minimum [%]
  Daylight = 13.0   ;Sunshine duration [h]

  ;Basic geographic information [Required]
  DEM = 1040.0   ;[x]Mean altitude [m] ;Scalar or raster file path

  ;Paths for input raster file [Required]
  Albedo_File = 'H:\Data\2017007004_IFEM_Albedo.tif'
  LST_File = 'H:\Data\2017007004_IFEM_Land_Surface_Temperature.tif'
  NDVI_File = 'H:\Data\2017007004_IFEM_NDVI.tif'

  ;Path for output raster file [Required]
  Out_path = 'H:\Data\'    ;[x]End with '\'

  ;===========================================


  ;==============Constants List===============
  
  DEFSYSV,'!Cp',1013.0      ;specific heat of air at constant pressure  [J kg-1 K-1]
  DEFSYSV,'!Sigma',5.67E-8  ;Stefan-Boltzmann constant  [W m-2 K-4]
  DEFSYSV,'!k',0.41         ;von Karman constant
  DEFSYSV,'!g',9.8          ;gravitational acceleration [m s-2]
  DEFSYSV,'!Z_Blend',200.0  ;Atmospheric Surface Layer(ASL) height [m]
  NDVI_max = 0.9            ;NDVI maximum for fully vegetated canopy
  NDVI_min = 0.08           ;NDVI minmum for bare soil
  
  ;===========================================


  ;=============Environment Settings===============

  ;Open Envi
  IF e EQ !NULL THEN BEGIN
    e = ENVI(/headless)
  ENDIF ELSE BEGIN
    e = ENVI(/current)
  ENDELSE

  ;Code specification
  COMPILE_OPT IDL2
  !EXCEPT=0

  ;Number of processes
  CPU,TPOOL_NTHREADS=!CPU.HW_NCPU/2

  ;==============================================


  Print,'======Program start======'

  ;==============Open Raster Files===============
  
  Print,'Start reading raster files...'
  
  ;Read NDVI file
  if ~FILE_TEST(NDVI_File) then begin
    print,"Cannot found NDVI raster file!"
    goto,END_Program
  endif
  NDVI_Raster = e.openraster(NDVI_file)
  NDVI = NDVI_Raster.getdata()
  NDVI[where(NDVI LT -1 OR NDVI GT 1)] = !VALUES.F_NAN
  
  ;Getting infromations of spatial reafrance
  Spatialref = NDVI_Raster.SPATIALREF
  NumRows = NDVI_Raster.NROWS
  Numcols = NDVI_Raster.NCOLUMNS

  ;Getting attributes of Raster
  Pixelsize = Spatialref.PIXEL_SIZE
  Coordsysstring = Spatialref.COORD_SYS_STR
  TiePointMap = Spatialref.TIE_POINT_MAP
  TiePointPixel = Spatialref.TIE_POINT_PIXEL

  ;Definition of spatial grid
  Coord_sys = ENVICOORDSYS(COORD_SYS_STR=Coordsysstring)
  Grid = ENVIGRIDDEFINITION(Coord_sys,$
    Pixel_size = Pixelsize,$
    Nrows = Numrows,$
    Ncolumns = Numcols,$
    TIE_POINT_MAP = TiePointMap,$
    TIE_POINT_PIXEL = TiePointPixel)
  
  NDVI_Raster.close
    
  ;Read albedo raster file
  IF ~FILE_TEST(Albedo_File) THEN BEGIN
    PRINT,"Cannot found Albedo raster file!"
    GOTO,END_Program
  ENDIF
  Albedo_raster_temp = e.openraster(Albedo_File)
  Albedo_raster = ENVISPATIALGRIDRASTER(Albedo_raster_temp,GRID_DEFINITION=Grid,RESAMPLING="bilinear")
  Albedo = Albedo_raster.getdata()
  Albedo[where(Albedo LT 0 OR Albedo GT 1)] = !VALUES.F_NAN
  Albedo_raster.close
  Albedo_raster_temp.close
  
  ;Read raster file for Land Surface Temperature
  IF ~FILE_TEST(LST_File) THEN BEGIN
    PRINT,"Cannot found LST raster file!"
    GOTO,END_Program
  ENDIF
  LST_raster_temp = e.OPENRASTER(LST_File)
  LST_raster = ENVISPATIALGRIDRASTER(LST_raster_temp,GRID_DEFINITION=Grid,RESAMPLING="bilinear")
  Trad = LST_raster.GETDATA()
  IF mean(Trad,/NAN) LT 200 THEN Trad = Trad + 273.15
  Trad[where(Trad LT 200 OR Trad GT 350)] = !VALUES.F_NAN
  LST_raster.CLOSE
  LST_raster_temp.CLOSE
  
  ;Read raster file for Air temperature
  IF Ta.TYPENAME EQ 'STRING' THEN BEGIN
    ;Check the Air temperature raster file
    IF ~FILE_TEST(Ta) THEN BEGIN
      PRINT,"Cannot found the Air temperature raster file! A scalar value can be allowed!"
      GOTO,END_Program
    ENDIF
    ;Open raster file
    Ta_raster_temp = e.OPENRASTER(Ta)
    Ta_raster = ENVISPATIALGRIDRASTER(Ta_raster_temp,GRID_DEFINITION=Grid,RESAMPLING="bilinear")
    Ta = Ta_raster.GETDATA()
    Ta_raster.CLOSE
    Ta_raster_temp.CLOSE
  ENDIF
  IF Ta LT 200 THEN Ta = Ta + 273.15

  
  ;Read VPD raster file
  IF VPD.TYPENAME EQ 'STRING' THEN BEGIN
    ;Check the VPD raster file
    IF ~FILE_TEST(VPD) THEN BEGIN
      PRINT,"Cannot found the VPD raster file! A scalar value can be allowed!"
      GOTO,END_Program
    ENDIF
    ;Open raster file
    VPD_raster_temp = e.OPENRASTER(VPD)
    VPD_raster = ENVISPATIALGRIDRASTER(VPD_raster_temp,GRID_DEFINITION=Grid,RESAMPLING="bilinear")
    VPD = VPD_raster.GETDATA()
    VPD_raster.CLOSE
    VPD_raster_temp.CLOSE
  ENDIF
  IF MEAN(VPD,/NAN) GT 10 THEN BEGIN
    Print,"The unit entered for VPD may be incorrect, please convert to kPa!"
    GOTO,END_Program
  ENDIF
    
  ;Mask NaN values
  NaN_Mask,NDVI,Albedo,Trad,Ta,VPD,Count
  
  ;Error message for input files
  IF Count GE N_ELEMENTS(NDVI)-1 THEN BEGIN
    PRINT,'Error: There is no intersection of input raster files!'
    GOTO,END_Program
  ENDIF
  
  ;================================================



  ;===================Calculation for basic parameters=================
  
  Print,'Start calculating...'
  
  ;Geographic longitude and latitude [rad]
  Ncolumns = Grid.NCOLUMNS
  Nrows = Grid.NROWS   
  SPATIALREF = Grid.SPATIALREF    ;Spatial reference
  Spatialref.CONvertfiletomap,Ncolumns/2,Nrows/2,Mapx,Mapy  ;Get the Map coordinates of the center pixel
  Spatialref.CONvertmaptolonlat,Mapx,Mapy,Lon,Lat           ;Get the Longitude and Latitude of the center pixel
  Lat = Lat[0]*!DTOR   ;Latitude [rad]
  Lon = Lon[0]*!DTOR   ;Longitude [rad]

  ;Solar radiation parameters calculation
  DOY = day-32+FIX(275*month/9)+2*FIX(3/(month+1))+FIX(month/100-(year MOD 4)/4+0.975)  ;Day of year [1,366]
  dr = 1+0.033*COS(DOY*2*!PI/365.0)   ;relative distance between sun and earth
  SolDec = 0.409*SIN(2*!PI*DOY/365.0-1.39)    ;Solar declination [rad]
  Sunrise_ang = ACOS(-TAN(Lat)*TAN(SolDec))   ;Sunset angle [rad]
  N = 24*Sunrise_ang/(!PI)         ;Maximum daylight hours [h]
  day_ang = ((2*!PI/365)*(DOY-1))  ;Day angle [rad]
  EquTime =(0.000075+0.001868*COS(day_ang)-0.032077*SIN(day_ang)-0.014615*COS(2*day_ang)-0.04089*SIN(2*day_ang))*229.18
  LoT= Time+(4*Lon*!RADEG/60)+EquTime/60   ;Local time [HH]
  Hour_ang = 15*(LoT-12)*(!PI/180)         ;Hour angle [rad]
  SolAlt = SIN(Lat)*SIN(SolDec)+COS(Lat)*COS(SolDec)*COS(Hour_ang)   ;Solar altitude cosine [rad]

  ;Surface characteristic parameters calculation
  es = 0.96  ;Surface emissivity of bare soil
  ec = 0.98  ;Surface emissivity of fully vegetated canopy
  FRACTIONAL_VAGETATION_COVER,NDVI,NDVI_max,NDVI_min,Fc,Fs
  ACTUAL_VAPOR_PRESSURE,Ta,VPD,Trad,eact   ;Vapor Pressure Deficit [kPa]
  ALBEDO_DECOMPOSITION,Albedo,SolAlt,DEM,Fc,Fs,Albedo_s,Albedo_c  ;Component Albedo calculation

  ;Atmospheric parameters calculation
  ;Refer to Allen et al.(2007),  http:/doi.org/10.1061/(asce)0733-9437(2007)133:4(380)
  Pa = 101.3*((293-0.0065*DEM)/293)^5.26                    ;Atmospheric pressure [kPa]
  Transmvty = TRANSMVTY_CALCULATION(Trad,Pa,eact,SolAlt)    ;Atmospheric transmvty
  Sd = 1367.0*SolAlt*Transmvty/dr^2     ;Downgoing shortwave radiation [W m-2]
  L_e = (2.501-2.361e-3*(Ta-273.15))    ;latent heat of vaporization [MJ kg-1]
  gamma_hy = (!CP*Pa)/(0.662*L_e*1E6)   ;Hygrometer constant [kPa ¡æ-1]
  rho = 1000*Pa/(1.01*Ta*287)           ;Air density [kg m-3]
  ;Slope of the saturation vapor pressure curve at air temperature [kPa ¡æ-1]
  Delta = (4098.0*0.6108*EXP((17.27*(Ta-273.5))/(Ta-273.15+237.3)))/(Ta-273.15+237.3)^2
  U_Blend = WIND_VELOCITY_ASL(NDVI,NDVI_max,u_ref,z_ref)   ;Wind velocity at ASL [m s-1]
  ;ea = 1.24*(10*eact/Ta)^(1.0/7.0)    ;Atmospheric emissivity refer to Brutsaert et al., 1975
  ea = 0.85*(-ALOG(Transmvty))^0.09    ;Atmospheric emissivity refer to Bastiaanssen et al.,1995

  ;====================================================================


  ;=====Iterative solving for Trapezoid framework and Soil Moisture Availability=====
  
  Print,'Start iterating...'

  ;Initial value setting
  Gamma_s = 0.315  ;Initial Gamma for bare soil
  Gamma_c = 0.05   ;Gamma for canopy
  Ts = Trad        ;Initial Surface temperature for bare soil
  Tc = Trad        ;Initial Surface temperature for canopy

  ;Net radiation independent on LST [W m-2]
  Rna_s = (1-Albedo_s)*Sd + es*ea*!SIGMA*Ta^4-es*!SIGMA*Ta^4
  Rna_c = (1-Albedo_c)*Sd + ec*ea*!SIGMA*Ta^4-ec*!SIGMA*Ta^4

  ;Radiation resistance for longwave [s m-1]
  rR_s = rho*!CP/(4*es*!SIGMA*Ta^3)
  rR_c = rho*!CP/(4*ec*!SIGMA*Ta^3)

  ;zT correction based on mean canopy height
  zT = ZT_CORRECTION(NDVI,NDVI_max,zT)

  FOR i=1,20 DO BEGIN

    ;Aerodynamic Resistance calculation
    Ras = AERODYNAMICS_RESISTANCE_SOIL(Ts,Ta,rho,u_Blend,zT)
    Rac = AERODYNAMICS_RESISTANCE_CANOPY(Tc,Ta,rho,u_Blend,zT)
    ;PRINT,'Ras mean: ' + STRTRIM(STRING(MEAN(Ras[where(NDVI GT 0)],/nan)),2)
    ;PRINT,'Rac mean: ' + STRTRIM(STRING(MEAN(Rac[where(NDVI GT 0)],/nan)),2)

    SUB = WHERE(NDVI LT 0)
    ;Independent trapezoid framework for each pixel
    ;Maximum surface temperature for dry bare soil
    Tsmax = (1-Gamma_s)*Rna_s/(rho*!CP*(1/Ras+(1-Gamma_s)/rR_s))+Ta
    ;Maximum surface temperature for dry canopy
    Tcmax = (1-Gamma_c)*Rna_c/(rho*!CP*(1/Rac+(1-Gamma_c)/rR_c))+Ta
    Tcmax[sub] = (1-0.5)*Rna_c[sub]/(rho[sub]*!CP*(1/Rac[sub]+(1-0.5)/rR_c[sub]))+Ta[sub]
    ;Minmum surface temperature for wet bare soil
    Tsmin = ((1-Gamma_s)*Rna_s-rho*!CP*VPD/(gamma_hy*Ras))/(rho*!CP*((gamma_hy+Delta)/(gamma_hy*Ras)+(1-Gamma_s)/rR_s))+Ta
    ;Minmum surface temperature for wet canopy
    Tcmin = ((1-Gamma_c)*Rna_c-rho*!CP*VPD/(gamma_hy*Rac))/(rho*!CP*((gamma_hy+Delta)/(gamma_hy*Rac)+(1-Gamma_c)/rR_c))+Ta
    Tcmin[sub] = ((1-0.5)*Rna_c[sub]-rho[sub]*!CP*VPD[sub]/(gamma_hy[sub]*Rac[sub]))/(rho[sub]*!CP*((gamma_hy[sub]+Delta[sub])/(gamma_hy[sub]*Rac[sub])+(1-0.5)/rR_c[sub]))+Ta[sub]

    ;Soil Moisture availability calculation
    Lambda_SM = SOIL_MOISTURE_AVAILABILITY(Trad,Tsmax,Tcmax,Tsmin,Tcmin,Fc,Fs)

    ;Gamma (the ratio of G/Rn)
    Gamma_sd = 0.2   ;Gamma for dry bare soil
    Gamma_sw = 0.5   ;Gamma for wet bare soil
    Gamma_s = Gamma_sw*Lambda_SM + Gamma_sd*(1-Lambda_SM)  ;Gamma for bare soil

    ;Interpolation of the Slope 
    K_cold = Tcmin-Tsmin     ;Slope of the cold edge
    K_warm = Tcmax-Tsmax     ;Slope of the warm edge
    K_slope = K_warm+Lambda_SM*(K_cold-K_warm)  ;Slope of the Lambda_SM isoline

    ;Let the slope K converge faster
    IF K_slope_old EQ !NULL THEN K_slope_old = K_slope*!VALUES.F_NAN
    K_slope = MEAN([[[K_slope]],[[K_slope_old]]],DIMENSION=3,/NAN)
    
    ;Prompt information
    K_Slope_MEAN = MEAN(K_slope[where(NDVI GT 0 AND Fc LT 0.8)],/NAN)
    Print,Strtrim(string(i),2)+'-th iterating, mean slope = ' + STRTRIM(STRING(K_Slope_MEAN),2)
    
    ;Break out of the loop
    IF ABS(MEAN(K_slope_old-K_slope,/NAN)) LT 0.1 THEN BREAK
    
    ;Trad decomposition
    Ts = Trad - K_slope*Fc   ;Surface tamperature for bare soil component [K]
    Tc = Ts + K_slope*1.0    ;Surface tamperature for fully canopy component [K]
    
    K_slope_old = K_slope

  ENDFOR

  ;==================================================================================


  Print,'Results calculation...'

  ;=====================Instantaneous Energy Fluxes====================

  ;Instantaneous latent heat flux [W m-2]
  LEs_max =(1-Gamma_s)*Rna_s - rho*!CP*(1/Ras+(1-Gamma_s)/rR_s)*(Tsmin-Ta)
  LEc_max =(1-Gamma_c)*Rna_c - rho*!CP*(1/Rac+(1-Gamma_c)/rR_c)*(Tcmin-Ta)
  LE_s = Lambda_SM*LEs_max     ;For bare soil component
  LE_c = Lambda_SM*LEc_max     ;For canopy component
  LE_inst = Fc*LE_c + Fs*LE_s  ;Aggregate surface

  ;Instantaneous net radiation [W m-2]
  Rn_s = (1-Albedo_s)*Sd+es*ea*!SIGMA*Ta^4-es*!SIGMA*Ts^4  ;For bare soil component
  Rn_c = (1-Albedo_c)*Sd+ec*ea*!SIGMA*Ta^4-ec*!SIGMA*Tc^4  ;For canopy component
  Rn = Fc*Rn_c + Fs*Rn_s                                   ;Aggregate surface
  
  ;Instantaneous soil heat flux [W m-2]
  G = Fc*Gamma_c*Rn_c + Fs*Gamma_s*Rn_s
  ;Aggregate ratio of G to Rn  [0,1]
  Gamma = Fc*Gamma_c+Fs*Gamma_s

  ;Instantaneous sensible heat flux [W m-2]
  H_s = (1-Gamma_s)*Rn_s-LE_s  ;For bare soil component
  H_c = (1-Gamma_c)*Rn_c-LE_c  ;For canopy component
  H = Fc*H_c + Fs*H_s          ;Aggregate
  
  ;==============================================================


  ;========Extrapolated to daily Evaporation and daily Transpiration=======

  ;IF Daylight equal to zero, don't extension to daily scale
  IF Daylight NE !NULL AND Daylight NE 0 THEN BEGIN

    ;Daily net radiation [MJ m-2 day-1]
    Rn_day = DAILY_NET_RADIATION(DOY,LAT,DEM,Daylight,Ta_max,Ta_min,RH_max,RH_min,Albedo,Pa,SolAlt)
    
    ;Evaporation Fraction [0,1]
    EF_s = LE_s/((1-Gamma_s)*Rn_s)  ;For bare soil component
    EF_c = LE_c/((1-Gamma_c)*Rn_c)  ;For canopy component
    
    ;Ratio of component net radiation to aggregate net radiation [0,1]
    q_s = Rn_s/Rn   ;For bare soil component
    q_c = Rn_c/Rn   ;For canopy component
    
    ;Daily Evaporation and Transpiration  [mm day-1]
    E_daily = (1-Fc)*EF_s*q_s*Rn_day/L_e   ;Evaporation
    T_daily = Fc*EF_c*q_c*Rn_day/L_e       ;Transpiration
    E_daily[WHERE(E_daily LT 0)] = 0
    T_daily[WHERE(T_daily LT 0)] = 0
    
    ;Daily Evapotranspiration [mm day-1]
    ET_daily = E_daily + T_daily 
    
  ENDIF

  ;========================================================================


  Print,'Saving outputs...'

  ;=====================Output Files========================

  ;Convert date to String
  YEAR = STRTRIM(STRING(YEAR),2)
  IF MONTH LT 10 THEN BEGIN
    MONTH = '0'+STRTRIM(STRING(MONTH),2)
  ENDIF ELSE BEGIN
    MONTH = STRTRIM(STRING(MONTH),2)
  ENDELSE
  IF DAY LT 10 THEN BEGIN
    DAY = '0'+STRTRIM(STRING(DAY),2)
  ENDIF ELSE BEGIN
    DAY = STRTRIM(STRING(DAY),2)
  ENDELSE

;  ;Output Gamma
;  Gamma_URI = Out_path+Year+Month+Day+'_IFEM_Gamma.dat'
;  Gamma_Raster = ENVIRASTER(Gamma,SPATIALREF=Grid.SPATIALREF,URI=Gamma_URI)
;  Gamma_Raster.SAVE
;  Gamma_Raster.CLOSE

;  ;Output Ta
;  Ta_URI = Out_path+Year+Month+Day+'_IFEM_Air_Temperature.dat'
;  Ta_Raster = ENVIRASTER(Ta-273.15,SPATIALREF=Grid.SPATIALREF,URI=Ta_URI)
;  Ta_Raster.SAVE
;  Ta_Raster.CLOSE

  ;Output Soil Moistuare Availability
  Lambda_SM_URI = Out_path+Year+Month+Day+'_IFEM_Soil_Moisture_Availability.dat'
  Lambda_SM_Raster = ENVIRASTER(Lambda_SM,SPATIALREF=Grid.SPATIALREF,URI=Lambda_SM_URI)
  Lambda_SM_Raster.SAVE
  Lambda_SM_Raster.CLOSE

  ;Output Rn
  Rn_URI = Out_path+Year+Month+Day+'_IFEM_net_radiation.dat'
  Rn_Raster = ENVIRASTER(Rn,SPATIALREF=Grid.SPATIALREF,URI=Rn_URI)
  Rn_Raster.SAVE
  Rn_Raster.CLOSE

;  ;Output LE_s and LE_c
;  LEs_URI = Out_path+Year+Month+Day+'_IFEM_Soil_Latent.dat'
;  LEs_Raster = ENVIRASTER(LE_s,SPATIALREF=Grid.SPATIALREF,URI=LEs_URI)
;  LES_Raster.SAVE
;  LES_Raster.CLOSE
;  LEc_URI = Out_path+Year+Month+Day+'_IFEM_Vegetation_Latent.dat'
;  LEc_Raster = ENVIRASTER(LE_c,SPATIALREF=Grid.SPATIALREF,URI=LEc_URI)
;  LEc_Raster.SAVE
;  LEc_Raster.CLOSE
  LE_URI = Out_path+Year+Month+Day+'_IFEM_Latent.dat'
  LE_Raster = ENVIRASTER(LE_inst,SPATIALREF=Grid.SPATIALREF,URI=LE_URI)
  LE_Raster.SAVE
  LE_Raster.CLOSE

  ;Output G
  G_URI = Out_path+Year+Month+Day+'_IFEM_Soil_Heat_Flux.dat'
  G_Raster = ENVIRASTER(G,SPATIALREF=Grid.SPATIALREF,URI=G_URI)
  G_Raster.SAVE
  G_Raster.CLOSE

  ;Output H
  H_URI = Out_path+Year+Month+Day+'_IFEM_Sensible_Heat_Flux.dat'
  H_Raster = ENVIRASTER(H,SPATIALREF=Grid.SPATIALREF,URI=H_URI)
  H_Raster.SAVE
  H_Raster.CLOSE

  IF Daylight NE !NULL AND Daylight NE 0 THEN BEGIN
    ;Output Component Evaporationn
    E_URI = Out_path+Year+Month+Day+'_IFEM_Daily_Evaporation.dat'
    E_Raster = ENVIRASTER(E_Daily,SPATIALREF=Grid.SPATIALREF,URI=E_URI)
    E_Raster.SAVE
    E_Raster.CLOSE

    ;Output Component Transpiration
    T_URI = Out_path+Year+Month+Day+'_IFEM_Daily_Transpiration.dat'
    T_Raster = ENVIRASTER(T_Daily,SPATIALREF=Grid.SPATIALREF,URI=T_URI)
    T_Raster.SAVE
    T_Raster.CLOSE
    
    ;Output Component Transpiration
    ET_URI = Out_path+Year+Month+Day+'_IFEM_Daily_ET.dat'
    ET_Raster = ENVIRASTER(ET_Daily,SPATIALREF=Grid.SPATIALREF,URI=ET_URI)
    ET_Raster.SAVE
    ET_Raster.CLOSE
  ENDIF


;  ;Output Trad
;  Trad_URI = Out_path+Year+Month+Day+'_IFEM_Trad.dat'
;  Trad_Raster = ENVIRASTER(Trad,SPATIALREF=Grid.SPATIALREF,URI=Trad_URI)
;  Trad_Raster.SAVE
;  Trad_Raster.CLOSE
;  ;Output Fc
;  Fc_URI = Out_path+Year+Month+Day+'_IFEM_Fc.dat'
;  Fc_Raster = ENVIRASTER(Fc,SPATIALREF=Grid.SPATIALREF,URI=Fc_URI)
;  Fc_Raster.SAVE
;  Fc_Raster.CLOSE
;  ;Output Tsmax
;  Tsmax_URI = Out_path+Year+Month+Day+'_IFEM_Tsmax.dat'
;  Tsmax_Raster = ENVIRASTER(Tsmax,SPATIALREF=Grid.SPATIALREF,URI=Tsmax_URI)
;  Tsmax_Raster.SAVE
;  Tsmax_Raster.CLOSE
;  ;Output Tsmin
;  Tsmin_URI = Out_path+Year+Month+Day+'_IFEM_Tsmin.dat'
;  Tsmin_Raster = ENVIRASTER(Tsmin,SPATIALREF=Grid.SPATIALREF,URI=Tsmin_URI)
;  Tsmin_Raster.SAVE
;  Tsmin_Raster.CLOSE
;  ;Output Tcmax
;  Tcmax_URI = Out_path+Year+Month+Day+'_IFEM_Tcmax.dat'
;  Tcmax_Raster = ENVIRASTER(Tcmax,SPATIALREF=Grid.SPATIALREF,URI=Tcmax_URI)
;  Tcmax_Raster.SAVE
;  Tcmax_Raster.CLOSE
;  ;Output Tcmin
;  Tcmin_URI = Out_path+Year+Month+Day+'_IFEM_Tcmin.dat'
;  Tcmin_Raster = ENVIRASTER(Tcmin,SPATIALREF=Grid.SPATIALREF,URI=Tcmin_URI)
;  Tcmin_Raster.SAVE
;  Tcmin_Raster.CLOSE
  
  ;=======================================================
  
  Print,'======End of program======'

  END_Program:
END

;========================================================================
;=========================End of Main Program============================
;========================================================================





;========================================================================
;=============================Functions Part=============================
;========================================================================


PRO NAN_MASK,NDVI,Albedo,Trad,Ta,VPD,Count
  ;Mask NaN values in rasters
  
  ;Create an index array
  Size = size(NDVI,/dime)
  Index = make_array(Size,VALUE=1)
  
  ;Mask NaN values
  Index[WHERE(FINITE(NDVI,/NAN),/NULL)] = 0
  Index[WHERE(FINITE(Albedo,/NAN),/NULL)] = 0
  Index[WHERE(FINITE(Trad,/NAN),/NULL)] = 0
  IF N_ELEMENTS(Ta) GT 1 THEN Index[WHERE(FINITE(Ta,/NAN),/NULL)] = 0
  IF N_ELEMENTS(VPD) GT 1 THEN Index[WHERE(FINITE(VPD,/NAN),/NULL)] = 0
  
  ;Align available dataset
  Sub = WHERE(Index EQ 0,Count,/NULL)
  NDVI[sub] = !VALUES.F_NAN
  Albedo[sub] = !VALUES.F_NAN
  Trad[sub] = !VALUES.F_NAN
  IF N_ELEMENTS(Ta) GT 1 THEN Ta[Sub] = !VALUES.F_NAN
  IF N_ELEMENTS(VPD) GT 1 THEN VPD[Sub] = !VALUES.F_NAN

END


FUNCTION DAILY_NET_RADIATION,DOY,Lat,DEM,Daylight,Ta_max,Ta_min,RH_max,RH_min,Albedo,Pa,SolAlt
;Calculation for daily net radiation according to FAO-56

  ;Convert temperature unit to Kelvin [K]
  IF Ta_max LT 100 THEN Ta_max = Ta_max+273.15
  IF Ta_min LT 100 THEN Ta_min = Ta_min+273.15

  ;Calculation for the solar position
  dr = 1+0.033*COS(DOY*2*!PI/365.0)         ;Relative distance between sun and earth
  SolDec = 0.409*SIN(2*!PI*DOY/365.0-1.39)  ;Solar declination [rad]
  Sunrise_ang = ACOS(-TAN(Lat)*TAN(SolDec)) ;Sunset angle [rad]

  ;Vapor pressure near surface
  e0_Ta_max = 0.611*EXP(17.27*(Ta_max-273.15)/(Ta_max-273.15+240.97))  ;[kPa]
  e0_Ta_min = 0.611*EXP(17.27*(Ta_min-273.15)/(Ta_min-273.15+240.97))  ;[kPa]
  esat = (e0_Ta_max+e0_Ta_min)/2  ;Saturation vapor pressure [kPa]
  eact = (e0_Ta_max*RH_min/100 + e0_Ta_min*RH_max/100)/2   ;Actual vapor pressure [kPa]

  ;Broad-band atmospheric transmissivity according to ASCE-EWRI(2005)
  Transmvty = TRANSMVTY_CALCULATION(Albedo,Pa,eact,SolAlt)

  ;Extraterrestrial radiation [MJ m-2 day-1]
  Ra = 24*60*0.082*dr*(Sunrise_ang*SIN(LAT)*SIN(SolDec)+COS(LAT)*COS(SolDec)*SIN(Sunrise_ang))/(!PI)
  N = 24*Sunrise_ang/(!PI)     ;Maximum daylight hours [h]

  ;Solar radition [MJ m-2 day-1]
  as = 0.15
  bs = 0.6
  Rs = (as+bs*Daylight/N)*Ra

  ;Daily net radiation [MJ m-2 day-1]
  RnS = (1-1.1*Albedo)*Rs    ;Net shortwave radiation [MJ m-2 day-1]
  Rso = Transmvty*Ra         ;Clear-sky solar radiation [MJ m-2 day-1]
  RnL = 4.903e-9*(Ta_max^4+Ta_min^4)/2*(0.34-0.14*SQRT(eact))*(1.35*Rs/Rso-0.35)  ;Net longwave radiation  [MJ m-2 day-1]
  Rn_day = RnS-RnL    ;Daily net radiation  [MJ m-2 day-1]
  Rn_day[where(Rn_day LT 0)] = 0
  
  RETURN,RN_DAY
  
END


PRO FRACTIONAL_VAGETATION_COVER,NDVI,NDVI_max,NDVI_min,Fc,Fs
;Calculation for Fractional Vagetation Cover (FVC) [0,1]
;Refer to Li et al.(2005), https://doi.org/10.1175/JHM464.1

  k = 0.6   ;Experience coefficient,[0.6,1.25]

  ;Fractional Vagetation Cover [0,1]
  Fc = 1-((NDVI_max - NDVI)/(NDVI_max - NDVI_min))^k
  Fc[WHERE(NDVI LT NDVI_min)] = 0
  Fc[WHERE(NDVI GT NDVI_max)] = 1

  ;Complement of FVC
  Fs = 1 - Fc

END


PRO ACTUAL_VAPOR_PRESSURE,Ta,VPD,Trad,eact
;Actual vapor pressure calculation [kPa]

  ;Saturated vapor pressure [kPa]
  esat = 0.6108*EXP(17.27*(Ta-273.15)/(Ta-273.15+237.3))
  
  ;Actual vapor pressure [kPa]
  eact = esat - VPD
  eact = eact > 0.1
  eact = eact < esat

END


PRO ALBEDO_DECOMPOSITION,Albedo,SolAlt,DEM,Fc,Fs,Albedo_s,Albedo_c
;Decomposition Albedo into bare soil and canopy components
;The Albedo of canopy simulated based on uniform leaf optical properties
;Refer to Jacobs and VanPul.(1990)
;http://doi.org/10.1016/0168-1923(90)90006-R

  ;Uniform leaf optical properties
  Al_PAR = 0.1  ;leaf reflectivity in PAR wavelength
  Tl_PAR = 0.1  ;leaf transmittance in PAR wavelength
  Al_NIR = 0.4  ;leaf reflectivity in NIR wavelength
  Tl_NIR = 0.4  ;leaf transmittance in NIR wavelength

  ;Fraction of the total incident shortwave irradiation that is direct
  Transmvty = 0.75+2.1e-5*MEAN(DEM,/nan)   ;Clear sky condition
  R = 0.847-1.61*SolAlt+1.04*SolAlt^2
  K = (1.47-R)/1.66
  IF Transmvty LE K THEN f_a = 1.66*Transmvty - 0.47
  IF Transmvty GT K THEN f_a = 1 - R

  ;Compontent albedo for fully vegetated canopy
  As_PAR = Al_PAR/(1-Tl_PAR+SQRT((1-Tl_PAR)^2-Al_PAR^2))  ;Reflection coefficient in PAR wavelength
  As_NIR = Al_NIR/(1-Tl_NIR+SQRT((1-Tl_NIR)^2-Al_NIR^2))  ;Reflection coefficient in NIR wavelength
  As = (As_PAR+As_NIR)/2  ;Reflection coefficient in shortwave
  As_k = 2*f_a/(1+1.6*SolAlt)+(1-f_a)*1.11
  Albedo_c = As_k*As  ;Albedo for fully vegetated surface

  ;Compontent albedo for bare soil
  Albedo_s = (Albedo-Albedo_c*Fc)/Fs  ;Albedo for bare soil component

  ;Avoid outliers from Decomposition
  As_dry_soil = 0.13
  As_wet_dry = 0.05
  Albedo_dry_soil = As_k*As_dry_soil
  Albedo_wet_soil = As_k*As_wet_dry
  Albedo_s[WHERE(Albedo_s GT 1.0)] = 1.0
  Albedo_s[WHERE(Albedo_s LT As_wet_dry)] = Albedo_wet_soil
  Albedo_s[WHERE(Fs LT 0.2 AND Albedo_s GT 0.2)] = 0.2

END


FUNCTION SOIL_MOISTURE_AVAILABILITY,Trad,Tsmax,Tcmax,Tsmin,Tcmin,Fc,Fs
;Calculation for Soil Moisture availability 
;Equal to evapotranspiration efficiency

  Tmin = (Tsmin*Fs+Tcmin*Fc)   ;Trad higher limit [K]
  Tmax = (Tsmax*Fs+Tcmax*Fc)   ;Trad lower limit [K]
  a = Tmax - Trad         ;the distance from pixel to warm edge
  b = Trad - Tmin         ;the distance from pixel to cold edge
  Lambda_SM = a/(a+b)     ;Soil Moisture availability

  ;Range constraints
  Lambda_SM[WHERE(Lambda_SM GT 1)] = 1
  Lambda_SM[WHERE(Lambda_SM LT 0)] = 0

  RETURN,Lambda_SM

END


FUNCTION WIND_VELOCITY_ASL,NDVI,NDVI_max,u_ref,z_ref
;Calculating the wind speed constant at the height of the Atmospheric Surface Layer(ASL)

  z0m = 0.005+0.5*(NDVI/NDVI_max)^2.5              ;Rough length [m]
  z0m_Mean = MEAN(z0m,/NAN)                        ;Assuming a scalar rough length at meteorological station location [m]
  Ustar_ref = !K*u_ref/ALOG(z_ref/z0m_Mean)        ;Assumed wind velocity profile
  U_Blend = Ustar_ref*ALOG(!z_Blend/z0m_Mean)/(!K) ;Wind velocity at ASL [m s-1]

  RETURN,U_blend
END


FUNCTION ZT_CORRECTION,NDVI,NDVI_max,zT
;IFEM assumpted a fully vegetated canopy height of 1 m according Long and Singh(2012)
;The Ta observation height would be corrected base on the height difference----zT_Diff

  hc_canopy = 1.0                      ;Assumping height of fully vegetated canopy [m]
  z0m = 0.005+0.5*(NDVI/NDVI_max)^2.5  ;Rough length [m]
  z0m_Mean = MEAN(z0m,/NAN)            ;Assuming a scalar rough length [m]
  hc_mean = z0m_Mean/0.136             ;Assuming canopy height [m]
  zT_Diff = zT-hc_mean                 ;Air tamperature observation height diff above canopy [m]
  IF zT_Diff LT 1.0 THEN zT = 2.0      ;Lowest height of air tamperature observation [m]
  IF zT_Diff GT 1.0 THEN zT = hc_canopy+zT_Diff  ;Correcting zT above canopy [m]
  RETURN,zT

END


FUNCTION TRANSMVTY_CALCULATION,Trad,Pa,eact,SolAlt
;Calculating broad-band atmospheric transmissivity according to ASCE-EWRI(2005)

  sub = WHERE(FINITE(Trad),Count)   ;Number of non-cloud cover pixels
  Kt = 1.0*Count/N_ELEMENTS(Trad)   ;Unitless turbidity coefficient estimation [0-1],Kt=1 for clean air
  W = 0.14*eact*Pa+2.1              ;Water in the atmosphere [mm] according to Garrison and Adler(1990)
  term1 = -0.00146*Pa/(Kt*SolAlt)
  term2 = 0.075*(W/SolAlt)^0.4
  Transmvty = 0.35+0.627*EXP(term1-term2)  ;Atmospheric transmissivity
  RETURN,Transmvty

END


FUNCTION AERODYNAMICS_RESISTANCE_SOIL,Ts,Ta,rho,U_Blend,zT
;Iterative solution for aerodynamic resistance of bare soil component
;Computition with Kondo's (1994) formula in conjunction with Paulson (1970) and Webb's (1970) formula

  ;Define bare soil properties
  d0 = 0.0     ;Zero plane displacement [m]
  z0m = 0.01   ;Roughness length for momentum transfer [m], according to Garratt et al.(1973)
  z0h = z0m/3  ;Roughness length for Heat transfer [m], according to Verseghy et al.(1993)

  ;Initial calculation
  Ustar = !K*U_Blend/ALOG((!z_Blend-d0)/z0m)   ;Friction velocity [m s-1]
  Ras = (ALOG((zT-d0)/z0h))/(Ustar*!K)         ;Aerodynamic resistance for non-vegetated surface [s m-1]

  ;Stability correction factors
  sz = SIZE(Ts,/DIMENSIONS)
  Psim_zm = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)
  Psim_z0m = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)
  Psih_zT = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)
  Psih_z0h = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)

  FOR i=1,10 DO BEGIN

    ;Sensible heat flux [W m-2]
    H = !CP*rho*(Ts-Ta)/Ras
    
    ;Monin-Obukhov length [m]
    L = -(rho*!CP*Ustar^3*Ta)/(!K*!G*H)

    ;Index for each pixel
    Index_GT = WHERE(L GT 0,N_GT)   ;Index for L greater than zero
    Index_LT = WHERE(L LT 0,N_LT)   ;Index for L less than zero
    Index_EQ = WHERE(L EQ 0,N_EQ)   ;Index for L equal to zero

    ;Unstable conditions (L<0)
    IF N_LT GT 0 THEN BEGIN
      x_zm = (1-16*(!z_Blend/L[Index_LT]))^0.25
      x_z0m = (1-16*(z0m/L[Index_LT]))^0.25
      Psim_zm[Index_LT] = 2.0*ALOG((1.0+x_zm)/2.0)+ALOG((1.0+x_zm^2)/2.0)-2.0*ATAN(x_zm)+0.5*!PI
      Psim_z0m[Index_LT] = 2.0*ALOG((1.0+X_z0m)/2.0)+ALOG((1.0+X_z0m^2)/2.0)-2.0*ATAN(X_z0m)+0.5*!PI
      X_zT = (1-16*zT/L[Index_LT])^0.25
      X_z0h = (1-16*z0h/L[Index_LT])^0.25
      Psih_zT[Index_LT] = 2.0*ALOG((1+X_zT^2)/2)
      Psih_z0h[Index_LT] = 2.0*ALOG((1+X_z0h^2)/2)
    ENDIF

    ;Stable conditions (L>0)
    IF N_GT GT 0 THEN BEGIN
      Psim_zm[Index_GT] = -5.0*2.0/L[Index_GT]
      Psim_z0m[Index_GT] = -5.0*z0m/L[Index_GT]
      Psih_zT[Index_GT] = -5.0*zT/L[Index_GT]
      Psih_z0h[Index_GT] = -5.0*z0h/L[Index_GT]
    ENDIF

    ;Neutral condition (L=0)
    IF N_EQ GT 0 THEN BEGIN
      Psim_zm[Index_EQ] = 0
      Psim_z0m[Index_EQ] = 0
      Psih_zT[Index_EQ] = 0
      Psih_z0h[Index_EQ] = 0
    ENDIF

    ;Correction of friction velocity and aerodynamic resistance
    Ustar = U_Blend*!K/(ALOG((!z_Blend-d0)/z0m)-Psim_zm+Psim_z0m)  ;Friction velocity [m s-1]
    Ras = (ALOG((zT-d0)/z0h)-Psih_zT+Psih_z0h)/(Ustar*!K)          ;Aerodynamic resistance for non-vegetated surface [s m-1]

  ENDFOR

  ;Aerodynamic Impedance Limitation [s m-1]
  Ras[WHERE(Ras GT 300)] = 300
  Ras[WHERE(Ras LT 10)] = 10
  RETURN,Ras

END



FUNCTION AERODYNAMICS_RESISTANCE_CANOPY,Tc,Ta,rho,U_Blend,zT
;Iterative solution for aerodynamic resistance of canopy component
;Computition with Kondo's (1994) formula in conjunction with Paulson (1970) and Webb's (1970) formula

  ;Define dense canopy properties
  hc = 1.0        ;Assumed height of canopy surface [m], according to Long and Singh. (2012)
  d0 = 0.66*hc    ;Zero plane displacement [m]
  z0m = 0.136*hc  ;Roughness length for momentum transfer [m], according to Garratt et al.(1973)
  z0h = z0m/7     ;Roughness length for Heat transfer [m], according to Verseghy et al.(1993)

  ;Initial calculation
  Ustar = !K*U_Blend/ALOG((!z_Blend-d0)/z0m)   ;Friction velocity [m s-1]
  Rac = (ALOG((zT-d0)/z0h))/(Ustar*!K)         ;Aerodynamic resistance for vegetated surface [s m-1]

  ;Stability correction factors
  sz = SIZE(Tc,/DIMENSIONS)
  Psim_zm = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)
  Psim_z0m = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)
  Psih_zT = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)
  Psih_z0h = MAKE_ARRAY(DIMENSION=sz,value = 0.0,/FLOAT)

  FOR i=1,10 DO BEGIN

    ;Sensible heat flux [W m-2]
    H = !CP*rho*(Tc-Ta)/Rac
    
    ;Monin-Obukhov length [m]
    L = -(rho*!CP*Ustar^3*Ta)/(!K*!G*H)
    
    ;Index for each pixel
    Index_GT = WHERE(L GT 0,N_GT)   ;Index for L greater than zero
    Index_LT = WHERE(L LT 0,N_LT)   ;Index for L less than zero
    Index_EQ = WHERE(L EQ 0,N_EQ)   ;Index for L equal to zero

    ;Unstable conditions (L<0)
    IF N_LT GT 0 THEN BEGIN
      x_zm = (1-16*(!z_Blend/L[Index_LT]))^0.25
      x_z0m = (1-16*(z0m/L[Index_LT]))^0.25
      Psim_zm[Index_LT] = 2.0*ALOG((1.0+x_zm)/2.0)+ALOG((1.0+x_zm^2)/2.0)-2.0*ATAN(x_zm)+0.5*!PI
      Psim_z0m[Index_LT] = 2.0*ALOG((1.0+X_z0m)/2.0)+ALOG((1.0+X_z0m^2)/2.0)-2.0*ATAN(X_z0m)+0.5*!PI
      X_zT = (1-16*zT/L[Index_LT])^0.25
      X_z0h = (1-16*z0h/L[Index_LT])^0.25
      Psih_zT[Index_LT] = 2.0*ALOG((1+X_zT^2)/2)
      Psih_z0h[Index_LT] = 2.0*ALOG((1+X_z0h^2)/2)
    ENDIF

    ;Stable conditions (L>0)
    IF N_GT GT 0 THEN BEGIN
      Psim_zm[Index_GT] = -5.0*2.0/L[Index_GT]
      Psim_z0m[Index_GT] = -5.0*z0m/L[Index_GT]
      Psih_zT[Index_GT] = -5.0*zT/L[Index_GT]
      Psih_z0h[Index_GT] = -5.0*z0h/L[Index_GT]
    ENDIF

    ;Neutral condition (L=0)
    IF N_EQ GT 0 THEN BEGIN
      Psim_zm[Index_EQ] = 0
      Psim_z0m[Index_EQ] = 0
      Psih_zT[Index_EQ] = 0
      Psih_z0h[Index_EQ] = 0
    ENDIF

    ;Correction of friction velocity and aerodynamic resistance
    Ustar = U_Blend*!K/(ALOG((!z_Blend-d0)/z0m)-Psim_zm+Psim_z0m)  ;Friction velocity [m s-1]
    Rac = (ALOG((zT-d0)/z0h)-Psih_zT+Psih_z0h)/(Ustar*!K)          ;Aerodynamic resistance for non-vegetated surface [s m-1]

  ENDFOR

  ;Aerodynamic Impedance Limitation [s m-1]
  Rac[WHERE(Rac GT 300)] = 300  
  Rac[WHERE(Rac LT 10)] = 10
  RETURN,Rac

END

;========================================================================
;==========================End of Functions Part=========================
;========================================================================

#ifndef _ARRAYS_H
#define _ARRAYS_H

//! arrays for array pool
enum ARRAY_NAME
{
  //! Integer Arrays
  MPFANUM = 0,      //0 Array with numbers of mpfa region
  EQLNUM,       //1 Array with numbers of equilibrium region
  SATNUM,           //2 Array with numbers of saturation region
  PVTNUM,           //3 Array with numbers of pvt region
  ACTNUM,           //4 active cells
  FIPNUM,           //5 Array with numbers of FIPNUM region
  BNDNUM,           //6 boundary condition
  EOSNUM,						//7 Array with numbers of EOS region
  ROCKNUM,          //8 Array with numbers of rock region
  ARR_TOTAL_INT ,
  //! Double Arrays
/*  
  PERMXY = 0,       //0 Permeability in XY direction (K_xy)
  PERMYZ,           //1 Permeability in YZ direction (K_yz)
  PERMZX,           //2 Permeability in ZX direction (K_zx)
  PERMX,            //3 Permeability in X direction (K_xx)
  PERMY,            //4 Permeability in Y direction (K_yy)
  PERMZ,            //5 Permeability in Z direction (K_zz)
  PORO,             //6 Porosity
  NTG,              //7 Net to gross ratio7
*/  
  SGL = 0,              //8 minimum saturation in SCAL for water oil system
  SGU,              //9 maximum saturation in SCAL for water oil system
  SOGCR,            //10 Critical oil saturation in water oil system
  SGCR,             //11 Critical water saturation in water oil system
  SWL,              //12 minimum saturation in SCAL for water oil system
  SWU,              //13 maximum saturation in SCAL for water oil system
  SOWCR,            //14 Critical oil saturation in water oil system
  SWCR,             //15 Critical water saturation in water oil system
  PBUB,             //16
  RS,               //17 initial gas oil ration array
  PCW,              //18
  SWATINIT,         //19
  SOIL,             //20 Oil saturation initial array
  SWAT,             //21 Water saturation initial array
  PRESSURE,         //22 Pressure initial array
/*
  DX,               //23 Lenght of blocks in X direction
  DY,               //24 Lenght of blocks in Y direction
  DZ,               //25 Lenght of blocks in Z direction
  MULTX,            //26
  MULTY,            //27
  MULTZ,            //28
  TOPS,             //29
  MULTPV,             //30 poro volume multiplier
*/
  SGAS,               //31
  ARR_TOTAL_DOUBLE,   //32
};

//! arrays for compositional model pool
enum
{
  ACF = 0,          // Compositional: accentric factor
  MW, 					    // Compositional: molecular weights
  PCRIT,						// Compositional: critical pressures
  TCRIT,						// Compositional: critical temperatures
  VCRIT,						// Compositional: critical volumes
  ZCRIT,            // Compositional: critical Z-factor
  RTEMP,            // Compositional: reservoir temperature
  ACFS,             // Compositional: accentric factor at surface
  MWS,              // Compositional: molecular weights at surface
  PCRITS,						// Compositional: critical pressures at surface
  TCRITS,						// Compositional: critical temperatures at surface
  VCRITS,						// Compositional: critical volumes at surface
  ZCRITS,           // Compositional: critical Z-factor at surface
  XMF,							// Compositional: cell initial oil composition
  YMF,							// Compositional: cell initial gas composition
  ZMF,							// Compositional: cell initial total composition
  ZI,               // Compositional: cell initial total composition for each region of EOS
  COMP_ARR_TOTAL_DOUBLE,
};
#endif //_ARRAYS_H

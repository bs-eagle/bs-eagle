/**
 * \file arrays_tables.cpp
 * \brief
 * \author Sergey Miryanov
 * \date 21.08.2009
 * */
#include "bs_bos_core_data_storage_stdafx.h"
#include "arrays_tables.h"

int i_pool_default_values[ARR_TOTAL_INT]=
{
  0,                 // MPFANUM
  0,                 // EQLNUM
  0,                 // SATNUM
  0,                 // PVTNUM
  1,                 // ACTNUM
  0,                 // FIPNUM
  0,                 // BNDNUM
  0,		             // EOSNUM
  0,                 // ROCKNUM
};
double d_pool_default_values[ARR_TOTAL_DOUBLE]=
{
/*
  0,                 // PERMXY
  0,                 // PERMYZ
  0,                 // PERMZX
  0,                 // PERMX
  0,                 // PERMY
  0,                 // PERMZ
  0,                 // PORO
  1,                 // NTG
*/
  0,                 // SGL
  0,                 // SGU
  0,                 // SOGCR
  0,                 // SGCR
  0,                 // SWL
  0,                 // SWU
  0,                 // SOWCR
  0,                 // SWCR
  0,                 // PBUB
  0,                 // RS
  0,                 // PCW
  0,                 // SWATINIT
  0,                 // SOIL
  0,                 // SWAT
  0,                 // PRESSURE
/*
  0,                 // DX
  0,                 // DY
  0,                 // DZ
  1,                 // MULTX
  1,                 // MULTY
  1,                 // MULTZ
  0,                 // TOPS
  1,                 // MULTPV
*/
  0,                 // SGAS
};

double c_pool_default_values[COMP_ARR_TOTAL_DOUBLE]=
{
  0,								 // ACF
  0,								 // MW
  0,								 // PCRIT
  0,								 // TCRIT
  0,								 // VCRIT
  0,								 // ZCRIT
  0,								 // RTEMP
  0,								 // ACFS
  0,                 // MWS
  0,                 // PCRITS
  0,                 // TCRITS
  0,                 // VCRITS
  0,                 // ZCRITS
  0,                 // XMF
  0,                 // YMF
  0,                 // ZMF
  0,                 // ZI
};
int i_pool_sizes[ARR_TOTAL_INT * 6] =    //!< begin and end chek array dimmention 1,2,3
{
  1, 0, 1, 0, 1, 0,    // MPFANUM
  1, 0, 1, 0, 1, 0,    // EQLNUM
  1, 0, 1, 0, 1, 0,    // SATNUM
  1, 0, 1, 0, 1, 0,    // PVTNUM
  1, 0, 1, 0, 1, 0,    // ACTNUM
  1, 0, 1, 0, 1, 0,    // FIPNUM
  1, 0, 1, 0, 1, 0,    // BNDNUM
  1, 0, 1, 0, 1, 0,    // EOSNUM
  1, 0, 1, 0, 1, 0     // ROCKNUM
};
int d_pool_sizes[ARR_TOTAL_DOUBLE * 6] = //!< begin and end chek array dimmension 1,2,3
{
/*
  1, 0, 1, 0, 1, 0,    // PERMXY
  1, 0, 1, 0, 1, 0,    // PERMYZ
  1, 0, 1, 0, 1, 0,    // PERMZX
  1, 0, 1, 0, 1, 0,    // PERMX
  1, 0, 1, 0, 1, 0,    // PERMY
  1, 0, 1, 0, 1, 0,    // PERMZ
  1, 0, 1, 0, 1, 0,    // PORO
  1, 0, 1, 0, 1, 0,    // NTG
*/
  1, 0, 1, 0, 1, 0,    // SGL
  1, 0, 1, 0, 1, 0,    // SGU
  1, 0, 1, 0, 1, 0,    // SOGCR
  1, 0, 1, 0, 1, 0,    // SGCR
  1, 0, 1, 0, 1, 0,    // SWL
  1, 0, 1, 0, 1, 0,    // SWU
  1, 0, 1, 0, 1, 0,    // SOWCR
  1, 0, 1, 0, 1, 0,    // SWCR
  1, 0, 1, 0, 1, 0,    // PBUB
  1, 0, 1, 0, 1, 0,    // RS
  1, 0, 1, 0, 1, 0,    // PCW
  1, 0, 1, 0, 1, 0,    // SWATINIT
  1, 0, 1, 0, 1, 0,    // SOIL
  1, 0, 1, 0, 1, 0,    // SWAT
  1, 0, 1, 0, 1, 0,    // PRESSURE
/*
  1, 0, 1, 0, 1, 0,    // DX
  1, 0, 1, 0, 1, 0,    // DY
  1, 0, 1, 0, 1, 0,    // DZ
  1, 0, 1, 0, 1, 0,    // MULTX
  1, 0, 1, 0, 1, 0,    // MULTY
  1, 0, 1, 0, 1, 0,    // MULTZ
  1, 0, 1, 0, 0, 1, 	  // TOPS
  1, 0, 1, 0, 1, 0,    // MULTPV
*/
  1, 0, 1, 0, 1, 0,    // SGAS
};

/*int c_pool_sizes[COMP_ARR_TOTAL_DOUBLE * COMP_ARRAY_POOL_TOTAL] = //!< begin and end chek array dimmension 1,2,3
{
 // nx_a  nx_b  ny_a  ny_b  nz_a  nz_b  nc_a  nc_b  nr_a  nr_b  ns_a  ns_b
    0,    1,    0,    1,    0,    1,    1,    0,	  1,    0,    0,    1,      // ACF    [nc * n_eos_rsv]
    0,    1,    0,    1,    0,    1,    1,    0,	  1,    0,    0,    1,      // MW		  [nc * n_eos_rsv]
    0,    1,    0,    1,    0,    1,    1,    0,	  1,    0,    0,    1,      // PCRIT  [nc * n_eos_rsv]
    0,    1,    0,    1,    0,    1,    1,    0,	  1,    0,    0,    1,      // TCRIT  [nc * n_eos_rsv]
    0,    1,    0,    1,    0,    1,    1,    0,	  1,    0,    0,    1,      // VCRIT  [nc * n_eos_rsv]
    0,    1,    0,    1,    0,    1,    1,    0,	  1,    0,    0,    1,      // ZCRIT  [nc * n_eos_rsv]
    0,    1,    0,    1,    0,    1,    0,    1,	  1,    0,    0,    1,      // RTEMP  [n_eos_rsv]
    0,    1,    0,    1,    0,    1,    1,    0,	  0,    1,    1,    0,      // ACFS    [nc * n_eos_srf]
    0,    1,    0,    1,    0,    1,    1,    0,	  0,    1,    1,    0,      // MWS		 [nc * n_eos_srf]
    0,    1,    0,    1,    0,    1,    1,    0,	  0,    1,    1,    0,      // PCRITS  [nc * n_eos_srf]
    0,    1,    0,    1,    0,    1,    1,    0,	  0,    1,    1,    0,      // TCRITS  [nc * n_eos_srf]
    0,    1,    0,    1,    0,    1,    1,    0,	  0,    1,    1,    0,      // VCRITS  [nc * n_eos_srf]
    0,    1,    0,    1,    0,    1,    1,    0,    0,    1,    1,    0,      // ZCRITS  [nc * n_eos_srf]
    1,    0,    1,    0,    1,    0,    1,    0,	  0,    1,    0,    1,      // XMF    [nx*ny*nz*nc]
    1,    0,    1,    0,    1,    0,    1,    0,	  0,    1,    0,    1,      // YMF    [nx*ny*nz*nc]
    1,    0,    1,    0,    1,    0,    1,    0,    0,    1,    0,    1,      // ZMF    [nx*ny*nz*nc]
    0,    1,    0,    1,    0,    1,    1,    0,    1,    0,    0,    1,      // ZI     [nc * n_eos_rsv]
};*/

const char *int_names_table [ARR_TOTAL_INT] =
{
  "MPFANUM",
  "EQLNUM",
  "SATNUM",
  "PVTNUM",
  "ACTNUM",
  "FIPNUM",
  "BNDNUM",
  "EOSNUM",
  "ROCKNUM"
};//!< Names Of int arrays

const char *double_names_table [ARR_TOTAL_DOUBLE] =
{
//  "PERMXY",
//  "PERMYZ",
//  "PERMZX",
//  "PERMX",
//  "PERMY",
//  "PERMZ",
//  "PORO",
//  "NTG",
  "SGL",
  "SGU",
  "SOGCR",
  "SGCR",
  "SWL",
  "SWU",
  "SOWCR",
  "SWCR",
  "PBUB",
  "RS",
  "PCW",
  "SWATINIT",
  "SOIL",
  "SWAT",
  "PRESSURE",
//  "DX",
//  "DY",
//  "DZ",
//  "MULTX",
//  "MULTY",
//  "MULTZ",
//  "TOPS",
//  "MULTPV",
  "SGAS",
};//!< Names Of double arrays


const char *comp_names_table [COMP_ARR_TOTAL_DOUBLE] =
{
  "ACF",              // Compositional: accentric factor
  "MW", 					    // Compositional: molecular weights
  "PCRIT",						// Compositional: critical pressures
  "TCRIT",						// Compositional: critical temperatures
  "VCRIT",						// Compositional: critical volumes
  "ZCRIT",            // Compositional: critical Z-factor
  "RTEMP",            // Compositional: reservoir temperature
  "ACFS",             // Compositional: accentric factor at surface
  "MWS",              // Compositional: molecular weights at surface
  "PCRITS",						// Compositional: critical pressures at surface
  "TCRITS",						// Compositional: critical temperatures at surface
  "VCRITS",						// Compositional: critical volumes at surface
  "ZCRITS",           // Compositional: critical Z-factor at surface
  "XMF",							// Compositional: cell initial oil composition
  "YMF",							// Compositional: cell initial gas composition
  "ZMF",							// Compositional: cell initial total composition
  "ZI",               // Compositional: cell initial total composition for each region of EOS
};


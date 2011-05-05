/**
 *       \file  well_events.h
 *      \brief  Declaration of well events
 *     \author  Morozov Andrey
 *       \date  07.06.2008
 *  \copyright  This source code is released under the terms of
 *              the BSD License. See LICENSE for more details.
 * */
#ifndef WELL_EVENTS_H_
#define WELL_EVENTS_H_

#include "event_base.h"
#include "well_event_params_decl.h"

namespace blue_sky
  {

  class well;
  class reservoir;
  class calc_model;

  namespace wells
    {
    class connection;
    class well_controller;
    class well_limit_operation;
  }

  /**
   * \class well_event
   * \brief Base class for well events
   * */
  class BS_API_PLUGIN well_event : public event_base
    {
    public:

      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
    public:

      //! dtor
      virtual ~well_event () {}

      /**
       * \brief  Applies event
       * \param  top Instance of reservoir_simulator
       * \param  msh
       * \param  calc_model
       * */
      virtual void
      apply (const sp_top &top, const sp_mesh_iface_t &msh, const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &data) const;

      /**
       * \brief  Returns name of well in event
       * \return Name of well in event
       * */
      virtual std::string
      get_well_name () const;

      /**
       * \brief  Returns name of group in event
       * \return Name of group in event
       * */
      virtual std::string
      get_group_name () const;

      /**
       * \brief  Returns name of event
       * \return Name of event
       * */
      virtual std::string
      get_event_name () const;

    protected:
      /**
       * \brief  Applies event
       * \param  top Instance of reservoir_simulator
       * \param  msh
       * \param  calc_model
       * */
      virtual void
      apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
      const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &data) const;

    };

//============================================================================================
  /**
   * \class   WELSPECS_event
   * \brief   WELSPECS keyword
   * \details WELLSPEC Introduces a new well 3
   *          This keyword introduces a new well, specifying the
   *          name and position of the well head, the BHP reference
   *          depth and the separator used.
   * */
  class BS_API_PLUGIN WELSPECS_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model            calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WELSPECS_event);
    public:

      MAIN_PARAMS (
        ((WELL_NAME,        "Well name", PT_STR))
        ((WELL_GROUP_NAME,  "Name of the group to which the well belongs", PT_STR))
        ((I,                "I - location of well head or heel", PT_INT))
        ((J,                "J - location of well head or heel", PT_INT))
        ((BHP_DEPTH,        "", PT_FLOAT))
      );

      //NEXT_LINE_PARAMS (
      //  ((WELL_NAME,        "Well name", PT_STR))
      //  ((WELL_GROUP_NAME,  "Name of the group to which the well belongs", PT_STR))
      //  ((I,                "I - location of well head or heel", PT_INT))
      //  ((J,                "J - location of well head or heel", PT_INT))
      //  ((BHP_DEPTH,        "", PT_FLOAT))
      //);

      //! destructor
      virtual ~WELSPECS_event () {}

      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   WELLCON_event
   * \brief   WELLCON keyword
   * \details Specifies the position and properties of one
   *          or more well completions. This must be
   *          entered after the WELLSPEC keyword defining
   *          the appropriate well.
   * */
  class BS_API_PLUGIN WELLCON_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model            calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WELLCON_event);
    public:

      MAIN_PARAMS (
        ((WELL_NAME,        "Well name", PT_STR))
        ((WELL_GROUP_NAME,  "Name of the group to which the well belongs", PT_STR))
        ((I,                "I - location of well head or heel", PT_INT))
        ((J,                "J - location of well head or heel", PT_INT))
        ((FIRST_LAYER,      "First layer", PT_INT))
        ((LAST_LAYER,       "Last layer", PT_INT))
        ((WELL_STATE,       "Open/shut flag for the well", PT_STR))
        ((WELL_DIAMETER,    "Well bore diameter at the connection", PT_FLOAT))
        ((SKIN,             "Skin factor", PT_FLOAT))
        ((KH,               "Effective Kh (permeability x thickness) value of the connection", PT_FLOAT))
        ((DIRECTION,        "Direction in which the well penetrates the grid block", PT_STR))
        );

      //! destructor
      virtual ~WELLCON_event () {}

      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================
  /**
   * \class   COMPDAT_event
   * \brief   COMPDAT keyword
   * \details Specifies the position and properties of one or
   *          more well completions. This must be entered after
   *          the WELSPECS keyword defining the appropriate well.
   * */
  class BS_API_PLUGIN COMPDAT_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef t_long                    index_t;
      typedef t_double                     item_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model            calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>          sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(COMPDAT_event);
    public:
      MAIN_PARAMS (
        ((WELL_NAME,          "Well name, well name root or well list name", PT_STR))
        ((I,                  "I - location of well head or heel", PT_INT))
        ((J,                  "J - location of well head or heel", PT_INT))
        ((FIRST_LAYER,        "K - location of upper connecting block in this set of data", PT_INT))
        ((LAST_LAYER,         "K - location of lower connecting block in this set of data", PT_INT))
        ((PERFORATION_STATUS, "Open/shut flag for the well", PT_STR))
        ((PARAM7,             "Well bore diameter at the connection", PT_INT))
        ((PERFORATION_FACTOR, "Transmissibility factor for the connection", PT_FLOAT))
        ((WELL_DIAMETER,      "Well bore diameter at the connection", PT_FLOAT))
        ((KH,                 "Effective Kh (permeability x thickness) value of the connection", PT_FLOAT))
        ((SKIN,               "Skin factor", PT_FLOAT))
        ((PARAM12,            "", PT_STR))
        ((DIRECTION,          "Direction in which the well penetrates the grid block", PT_STR))
        ((PARAM13,            "", PT_STR))
        ((SEG_NUMBER,         "Segment which contains this perforation", PT_INT))
        );

      //! destructor
      virtual ~COMPDAT_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   WCONPROD_event
   * \brief   WCONPROD keyword
   * \details Control data for production wells
   * */
  class BS_API_PLUGIN WCONPROD_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model            calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WCONPROD_event);
    public:
      MAIN_PARAMS (
        ((WELL_NAME,          "Well name, well name root or well list name", PT_STR))
        ((WELL_STATE,         "Open/shut flag for the well", PT_STR))
        ((WELL_CONTROL_MODE,  "Control mode", PT_STR))
        ((OIL_RATE,           "Oil rate target or upper limit", PT_FLOAT))
        ((WATER_RATE,         "Water rate target or upper limit", PT_FLOAT))
        ((GAS_RATE,           "Gas rate target or upper limit", PT_FLOAT))
        ((OUTER_RATE,         "Liquid rate target or upper limit", PT_FLOAT))
        ((INNER_RATE,         "Reservoir fluid volume rate target or upper limit", PT_FLOAT))
        ((BHP,                "BHP target or lower limit", PT_FLOAT))
      );

      //! destructor
      virtual ~WCONPROD_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };


//============================================================================================

  /**
   * \class   WCONHIST_event
   * \brief   WCONHIST keyword
   * \details Observed rates for history matching production wells.
   *          This keyword is used in place of WCONPROD to declare
   *          production wells as special history matching wells,
   *          and to enter their observed flow rates (and optionally
   *          their measured BHP and THP values). The equivalent
   *          keyword for defining history matching injection wells
   *          is WCONINJH.
   * */
  class BS_API_PLUGIN WCONHIST_event: public well_event
    {
    public:
      typedef reservoir                             reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WCONHIST_event);
    public:

      //TODO: check for one more param!
      MAIN_PARAMS (
        ((WELL_NAME,          "Well name, well name root or well list name", PT_STR))
        ((WELL_STATE,         "Open/shut flag for the well", PT_STR))
        ((WELL_CONTROL_MODE,  "Control mode", PT_STR))
        ((OIL_RATE,           "Observed oil production rate", PT_FLOAT))
        ((WATER_RATE,         "Observed water production rate", PT_FLOAT))
        ((GAS_RATE,           "Observed gas production rate", PT_FLOAT))
        ((PARAM7,             "", PT_STR))
        ((PARAM8,             "", PT_STR))
        ((PARAM9,             "", PT_STR))
        ((BHP,                "Observed bottom hole pressure (BHP)", PT_FLOAT))
      );

      //! destructor
      virtual ~WCONHIST_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };


//============================================================================================

  /**
   * \class   WCONINJE_event
   * \brief   WCONINJE keyword
   * \details Injection well control data, with no group control.
   *          This keyword can be used to set individual control
   *          targets and limits for injection wells
   * */
  class BS_API_PLUGIN WCONINJE_event: public well_event
    {
    public:
      typedef reservoir                             reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WCONINJE_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((WELL_NAME,          "Well name, well name root or well list name", PT_STR))
        ((INJECTOR_TYPE,      "Injector type", PT_STR))
        ((WELL_STATE,         "Open/shut flag for the well", PT_STR))
        ((WELL_CONTROL_MODE,  "Control mode", PT_STR))
        ((OUTER_RATE,         "Surface flow rate target or upper limit", PT_FLOAT))
        ((INNER_RATE,         "Reservoir fluid volume rate target or upper limit", PT_FLOAT))
        ((BHP,                "BHP target or lower limit", PT_FLOAT))
        ((PARAM1,             "dummy for bug", PT_STR))
      );

      //! destructor
      virtual ~WCONINJE_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;
      //float def_outer_rate(){return 0;}
      //float def_inner_rate(){return 0;}

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   WECON_event
   * \brief   WECON keyword
   * \details Economic limit data for production wells
   * */
  class BS_API_PLUGIN WECON_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WECON_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((WELL_NAME,    "Well name, well name root or well list name", PT_STR))
        ((MIN_OIL_RATE, "Minimum oil production rate", PT_FLOAT))
        ((WATER_CUT,    "Maximum water cut", PT_FLOAT))
        ((OPERATION,    "Workover procedure on exceeding limit", PT_STR))
      );

      //! destructor
      virtual ~WECON_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   WECONINJ_event
   * \brief   WECONINJ keyword
   * \details Economic limit data for injection wells
   * */
  class BS_API_PLUGIN WECONINJ_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WECONINJ_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((WELL_NAME,  "Well name, well name root or well list name", PT_STR))
        ((MIN_RATE,   "Minimum economic water injection rate", PT_FLOAT))
      );

      //! destructor
      virtual ~WECONINJ_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   WEFAC_event
   * \brief   WEFAC keyword
   * \details Sets well efficiency factors (for downtime)
   * */
  class BS_API_PLUGIN WEFAC_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WEFAC_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((WELL_NAME,        "Well name, well name root or well list name", PT_STR))
        ((OPERATION_FACTOR, "Efficiency factor for the well", PT_FLOAT))
      );

      //! destructor
      virtual ~WEFAC_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   WELTARG_event
   * \brief   WELTARG keyword
   * \details Resets a well operating target or limit.
   *          This keyword can be used to reset a target or
   *          limit value for a well, without having to re-specify
   *          all the other quantities required by the control
   *          keywords WCONPROD or WCONINJE. These other
   *          quantities are left unchanged, including the open/shut
   *          status of the well. If the well has been declared a
   *          history matching well (see keywords WCONHIST and WCONINJH)
   *          the WELTARG keyword may be used to modify its BHP limit,
   *          VFP table number, and artificial lift quantity. The other
   *          quantities should not be modified with this keyword.
   * */
  class BS_API_PLUGIN WELTARG_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WELTARG_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((WELL_NAME,    "Well name, well name root or well list name", PT_STR))
        ((WELL_CONTROL, "Definition of the control or constraint quantity to be changed", PT_STR))
        ((VALUE,        "New value", PT_FLOAT))
      );

      //! destructor
      virtual ~WELTARG_event () {}
      void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
           const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   WPIMULT_event
   * \brief   WPIMULT keyword
   * \details Multiplies well connection factors by a given value.
   *          This keyword can be used to multiply the connection
   *          transmissibility factors of selected well connections
   *          by a specified value. To multiply the transmissibility
   *          factors of all the connections in a well, leave
   *          items 3 - 7 defaulted. To multiply the transmissibility
   *          factors of a subset of connections in a well, you can
   *          identify the subset by their I,J,K location (items 3 - 5).
   *          A subset of connections can also be identified by their
   *          completion.
   * */
  class BS_API_PLUGIN WPIMULT_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(WPIMULT_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((WELL_NAME,    "Well name, well name root or well list name", PT_STR))
        ((PERM_FACTOR,  "Multiplier to act on the well’s connection transmissibility factors and Kh values", PT_FLOAT))
        ((I,            "I - location of connecting grid block(s)", PT_INT))
        ((J,            "J - location of connecting grid block(s)", PT_INT))
        ((K,            "K - location of connecting grid block(s)", PT_INT))
        ((Z1,           "Number of first completion in range", PT_INT))
        ((Z2,           "Number of last completion in range", PT_INT))
      );

      //! destructor
      virtual ~WPIMULT_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_well_name () const;
      virtual std::string get_event_name () const;
    };

//============================================================================================

  /**
   * \class   COMPENSATION_event
   * \brief   COMPENSATION keyword
   * \details Defines a compensation type and rate for well group
   * */
  class BS_API_PLUGIN COMPENSATION_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(COMPENSATION_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((WELL_GROUP,         "Name of the group to which the well belongs", PT_STR))
        ((COMPENSATION_RATE,  "Compensation rate", PT_FLOAT))
        ((COMPENSATION_TYPE,  "Compensation type", PT_STR))
      );

      //! destructor
      virtual ~COMPENSATION_event () {}
      virtual void apply_internal (const sp_top &top, const sp_mesh_iface_t &msh,
                   const sp_calc_model_t &calc_model, const smart_ptr <idata, true> &/*data*/) const;

    protected:
      virtual std::string get_event_name () const;
    };

  /**
   * \class   PERMFRAC_event
   * \brief   PERMFRAC keyword
   * \details Defines the fracture in cells
   * \todo    Obsolete
   * */
  class BS_API_PLUGIN PERMFRAC_event: public well_event
    {
    public:
      typedef reservoir                           reservoir_t;

      typedef well                                well_t;
      typedef wells::connection                   connection_t;
      typedef wells::well_controller              well_controller_t;
      typedef wells::well_limit_operation                     well_limit_operation_t;
      typedef rs_mesh_iface                       mesh_iface_t;
      typedef rs_smesh_iface                      smesh_iface_t;
      typedef calc_model                          calc_model_t;

      typedef smart_ptr <well_t, true>                        sp_well_t;
      typedef smart_ptr <connection_t, true>                  sp_connection_t;
      typedef smart_ptr <mesh_iface_t, true>                  sp_mesh_iface_t;
      typedef smart_ptr <smesh_iface_t, true>                 sp_smesh_iface_t;
      typedef smart_ptr <well_controller_t, true>             sp_well_controller_t;
      typedef smart_ptr <well_limit_operation_t, true>        sp_well_limit_operation_t;
      typedef smart_ptr <calc_model_t, true>                  sp_calc_model_t;

      typedef smart_ptr <reservoir_t, true>                   sp_top;
      BLUE_SKY_TYPE_DECL(PERMFRAC_event);
    public:

      //TODO: there are a number of additional params in eclipse !
      MAIN_PARAMS (
        ((I,                  "I - location", PT_INT))
        ((J,                  "J - location", PT_INT))
        ((KW1,                "fracture head location", PT_INT))
        ((KW2,                "fracture tail location", PT_INT))
        ((NUM_RIGHT,          "Number of cell placed at right side of central fracture cell", PT_INT))
        ((NUM_LEFT,           "Number of cell placed at left side of central fracture cell", PT_INT))
        ((FRAC_PERM,          "Fracture permeability", PT_FLOAT))
        ((SKIN,               "Skin factor", PT_FLOAT))
        ((FRACTURE_DIRECTION, "Fracture direction", PT_STR))
        ((PARAM10,            "", PT_STR))
      );

      //! destructor
      virtual ~PERMFRAC_event () {}

    protected:
      virtual std::string get_event_name () const;
    };

}//namespace blue_sky

#endif//WELL_EVENTS_H_


/**
 * \file well_rate_control_impl_type.h
 * \brief types for parametrization well_rate_control_impl
 * \author Sergey Miryanov
 * \date 25.11.2008
 * */
#ifndef BS_WELLS_WELL_RATE_CONTROL_IMPL_TYPE_H_
#define BS_WELLS_WELL_RATE_CONTROL_IMPL_TYPE_H_

namespace blue_sky
{
  namespace wells
  {
    template <typename rate_type, typename bhp_deriv_type, typename rate_deriv_type>
    struct compute_impl
    {
      typedef typename rate_type::strategy_t    strategy_t;
      typedef typename rate_type::mobility_t    mobility_t;
      typedef rate_type                         rate_t;
      typedef bhp_deriv_type                    bhp_deriv_t;
      typedef rate_deriv_type                   rate_deriv_t;
    };

    template <typename strategy_type>
    struct compute_dummy_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <dummy_deriv <mobility_calc_inj_t>, dummy_deriv <mobility_calc_inj_t>, dummy_deriv  <mobility_calc_inj_t> >    inj_impl_t;
      typedef compute_impl <dummy_deriv <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t> >  prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_bhp_3p_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_3p <mobility_calc_inj_t>, compute_bhp_deriv_3p <mobility_calc_inj_t>, dummy_deriv  <mobility_calc_inj_t> >   inj_impl_t;
      typedef compute_impl <compute_rate_3p <mobility_calc_prod_t>, compute_bhp_deriv_3p <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_rate_3p_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_3p <mobility_calc_inj_t>, compute_bhp_deriv_3p <mobility_calc_inj_t>, compute_rate_deriv_3p <mobility_calc_inj_t> >    inj_impl_t;
      typedef compute_impl <compute_rate_3p <mobility_calc_prod_t>, compute_bhp_deriv_3p <mobility_calc_prod_t>, compute_rate_deriv_3p <mobility_calc_prod_t> > prod_impl_t;
    };
 
    template <typename strategy_type>
    struct compute_bhp_2p_ow_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_2p_ow <mobility_calc_inj_t>, compute_bhp_deriv_2p_ow <mobility_calc_inj_t>, dummy_deriv  <mobility_calc_inj_t> >   inj_impl_t;
      typedef compute_impl <compute_rate_2p_ow <mobility_calc_prod_t>, compute_bhp_deriv_2p_ow <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_rate_2p_ow_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_2p_ow <mobility_calc_inj_t>, compute_bhp_deriv_2p_ow <mobility_calc_inj_t>, compute_rate_deriv_2p_ow <mobility_calc_inj_t> >    inj_impl_t;
      typedef compute_impl <compute_rate_2p_ow <mobility_calc_prod_t>, compute_bhp_deriv_2p_ow <mobility_calc_prod_t>, compute_rate_deriv_2p_ow <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_bhp_2p_og_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_2p_og <mobility_calc_inj_t>, compute_bhp_deriv_2p_og <mobility_calc_inj_t>, dummy_deriv  <mobility_calc_inj_t> >   inj_impl_t;
      typedef compute_impl <compute_rate_2p_og <mobility_calc_prod_t>, compute_bhp_deriv_2p_og <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_rate_2p_og_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_2p_og <mobility_calc_inj_t>, compute_bhp_deriv_2p_og <mobility_calc_inj_t>, compute_rate_deriv_2p_og <mobility_calc_inj_t> >    inj_impl_t;
      typedef compute_impl <compute_rate_2p_og <mobility_calc_prod_t>, compute_bhp_deriv_2p_og <mobility_calc_prod_t>, compute_rate_deriv_2p_og <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_bhp_1p_o_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_1p_o <mobility_calc_inj_t>, compute_bhp_deriv_1p_o <mobility_calc_inj_t>, dummy_deriv  <mobility_calc_inj_t> >   inj_impl_t;
      typedef compute_impl <compute_rate_1p_o <mobility_calc_prod_t>, compute_bhp_deriv_1p_o <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_rate_1p_o_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_1p_o <mobility_calc_inj_t>, compute_bhp_deriv_1p_o <mobility_calc_inj_t>, compute_rate_deriv_1p_o <mobility_calc_inj_t> >    inj_impl_t;
      typedef compute_impl <compute_rate_1p_o <mobility_calc_prod_t>, compute_bhp_deriv_1p_o <mobility_calc_prod_t>, compute_rate_deriv_1p_o <mobility_calc_prod_t> > prod_impl_t;
    };
    template <typename strategy_type>
    struct compute_bhp_1p_w_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_1p_w <mobility_calc_inj_t>, compute_bhp_deriv_1p_w <mobility_calc_inj_t>, dummy_deriv  <mobility_calc_inj_t> >   inj_impl_t;
      typedef compute_impl <compute_rate_1p_w <mobility_calc_prod_t>, compute_bhp_deriv_1p_w <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_rate_1p_w_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_1p_w <mobility_calc_inj_t>, compute_bhp_deriv_1p_w <mobility_calc_inj_t>, compute_rate_deriv_1p_w <mobility_calc_inj_t> >    inj_impl_t;
      typedef compute_impl <compute_rate_1p_w <mobility_calc_prod_t>, compute_bhp_deriv_1p_w <mobility_calc_prod_t>, compute_rate_deriv_1p_w <mobility_calc_prod_t> > prod_impl_t;
    };
    template <typename strategy_type>
    struct compute_bhp_1p_g_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_1p_g <mobility_calc_inj_t>, compute_bhp_deriv_1p_g <mobility_calc_inj_t>, dummy_deriv  <mobility_calc_inj_t> >   inj_impl_t;
      typedef compute_impl <compute_rate_1p_g <mobility_calc_prod_t>, compute_bhp_deriv_1p_g <mobility_calc_prod_t>, dummy_deriv <mobility_calc_prod_t> > prod_impl_t;
    };

    template <typename strategy_type>
    struct compute_rate_1p_g_type
    {
      typedef strategy_type                   strategy_t;
      typedef mobility_calc_inj <strategy_t>  mobility_calc_inj_t;
      typedef mobility_calc_prod <strategy_t> mobility_calc_prod_t;
      typedef compute_impl <compute_rate_1p_g <mobility_calc_inj_t>, compute_bhp_deriv_1p_g <mobility_calc_inj_t>, compute_rate_deriv_1p_g <mobility_calc_inj_t> >    inj_impl_t;
      typedef compute_impl <compute_rate_1p_g <mobility_calc_prod_t>, compute_bhp_deriv_1p_g <mobility_calc_prod_t>, compute_rate_deriv_1p_g <mobility_calc_prod_t> > prod_impl_t;
    };
  } // namespace wells
} // namespace blue_sky


#endif  // #ifndef BS_WELLS_WELL_RATE_CONTROL_TYPE_H_


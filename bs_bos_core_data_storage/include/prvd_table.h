#ifndef PRVD_TABLE_H
#define PRVD_TABLE_H

/*!
  \file prvd_table.h
  \brief declaration of table class for storing prvd
				 and rsvd keyword content
*/

#ifdef EPS
#undef EPS
#endif // EPS
#define EPS 1.e-10

namespace blue_sky
  {
  /*!
  	\class val_versus_depth
  	\brief class for storing values versus depth keyword content
  */
  class BS_API_PLUGIN val_vs_depth
    {
    public:
      // typedefs
      typedef std::vector<double> d_vec;

      //! default ctor
      val_vs_depth();
      //! dtor
      ~val_vs_depth();

      //! length of the table
      size_t get_table_len () const;
      //! set table len
      void set_table_len (unsigned);

      //! return reference to the depth array
      d_vec &tdepth();
      //! -//- values array
      d_vec &tvalues();

      //! calculate pressure using linear interpolation
      /*inline*/ double interpolate_linear (double) const;

      //! check monotonic of pressure and depth
      int check_monotonic (int region_number = 1) const;

      // Data
    protected:
      d_vec depth,values;
    };


}

#endif // PRVD_TABLE_H

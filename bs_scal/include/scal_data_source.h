/**
 * \file scal_data_source.h
 * \brief data sources for scal tables
 * \author Sergey Miryanov
 * \date 21.05.2008
 * */
#ifndef BS_SCAL_DATA_SOURCE
#define BS_SCAL_DATA_SOURCE

namespace blue_sky
  {

  class spof_data_source
    {
    public:

      virtual ~spof_data_source () {}

      virtual double Sp () const                = 0;
      virtual double So () const                = 0;
      virtual double Krp () const               = 0;
      virtual double Krop () const              = 0;
      virtual double Pcp () const               = 0;

      virtual void first ()                     = 0;
      virtual void next ()                      = 0;
      virtual bool end ()                       = 0;
    };

  class spfn_data_source
    {
    public:

      virtual ~spfn_data_source () {}

      virtual double Sp () const                = 0;
      virtual double Krp () const               = 0;
      virtual double Pcp () const               = 0;

      virtual void first ()                     = 0;
      virtual void next ()                      = 0;
      virtual bool end ()                       = 0;
    };

  class sof3_data_source
    {
    public:

      virtual ~sof3_data_source () {}

      virtual double So () const                = 0;
      virtual double Krow () const              = 0;
      virtual double Krog () const              = 0;

      virtual void first ()                     = 0;
      virtual void next ()                      = 0;
      virtual bool end ()                       = 0;
    };

  struct input_scal_region
    {
      typedef shared_vector <double> vector_t;

      vector_t    So;                       ///< oil saturation
      vector_t    Sp;                       ///< phase (water or gas) saturation
      vector_t    Krp;
      vector_t    Krop;
      vector_t    Pcp;                      ///< phase (water or gas) capillarity, sizeof same to sizeof Sp

      void clear ()
      {
        Sp.clear ();
        So.clear ();
        Krp.clear ();
        Krop.clear ();
        Pcp.clear ();
      }

      size_t size () const
        {
          return Sp.size () + So.size () + Krp.size () + Krop.size () + Pcp.size ();
        }
    };

} // namespace blue_sky
#endif // #ifndef BS_SCAL_DATA_SOURCE

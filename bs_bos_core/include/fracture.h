#ifndef FRACTURE_H
#define FRACTURE_H

namespace blue_sky
  {
  class fracture
    {
    public:
      //! default ctor
      fracture();
      //! copy ctor
      fracture(const fracture&);
      //! dtor
      ~fracture();

      //! Assignment operator
      fracture &operator=(const fracture&);

      //! initialization
      void init();

      //! for get and set well's name
      std::string &tname();

      //! check data
      int check_data() const;

      //data
      int dim[4];         //!< coordinates accessed as array
      double frac_len;    //!< fracture lenght
      double frac_angle;  //!< fracture angle
      double frac_skin;   //!< fracture skin factor
      double frac_perm;   //!< fracture permeability
      double frac_width;  //!< fracture width

      int frac_flag;

      int connection_type;

      //! FRACTURE creation time as amount of days from
      //! starting date, negative value means not
      //! initialized
      double start_time;

      int is_activated;    //!< show whether fracture already activated
      int ignore_fracture; //!< ignore this fracture (for example, if
      //!  fracture placed in inactive block)

    private:
      //! name of well with fracture
      std::string name;
    };
};

#endif // FRACTURE_H

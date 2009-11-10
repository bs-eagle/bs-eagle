/** 
 * \file make_me_happy.h
 * \brief 
 * \author Sergey Miryanov
 * \date 08.05.2009
 * */
#ifndef BS_BOS_CORE_BASE_MAKE_ME_HAPPY_H_
#define BS_BOS_CORE_BASE_MAKE_ME_HAPPY_H_

#define MAKE_ME_HAPPY(type_t, base_t, type_name)      \
  BLUE_SKY_TYPE_STD_CREATE_T_MEM (type_t)             \
  BLUE_SKY_TYPE_STD_COPY_T_MEM (type_t)               \
  BS_LOCK_THIS_DECL(type_t);                          \
                                                      \
  friend class type_descriptor;                       \
                                                      \
  static const type_descriptor &                      \
  td_maker (const std::string &stype_postfix)         \
  {                                                   \
    static blue_sky::type_descriptor td(Loki::Type2Type<type_t> ()  \
      , Loki::Type2Type <base_t> ()                   \
      , Loki::Int2Type <false> ()                     \
      , stype_postfix                                 \
      , ""                                            \
      , "");                                          \
                                                      \
    return td;                                        \
  }                                                   \
                                                      \
  static blue_sky::type_descriptor bs_type()          \
  {                                                   \
    return td_maker (std::string (type_name) + "_" + BOOST_CURRENT_FUNCTION); \
  }                                                   \
  virtual blue_sky::type_descriptor bs_resolve_type() const \
  {                                                   \
    return bs_type ();                                \
  }                                                   \
                                                      \
  type_t (bs_type_ctor_param param)                   \
  : bs_refcounter (), base_t (param)                                    \
  {                                                   \
  }                                                   \
  type_t (const type_t &solver)                       \
  : bs_refcounter (), base_t (solver)                                   \
  {                                                   \
  }


#endif  // #ifndef BS_BOS_CORE_BASE_MAKE_ME_HAPPY_H_

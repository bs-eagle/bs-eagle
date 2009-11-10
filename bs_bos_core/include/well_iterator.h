///**
// * \file well_iterator.h
// * \brief class for iterate well in facility_manager
// * \author Sergey Miryanov
// * \date 08.08.2008
// * */
//#ifndef BS_WELL_ITERATOR_H_
//#define BS_WELL_ITERATOR_H_
//
//#include <iterator>
//
//namespace blue_sky
//  {
//
//  template <typename strategy_t>
//  class facility_manager;
//
//  template <typename strategy_t>
//  class well;
//
//  template <typename list_t, typename iterator_t, typename value_type_t>
//  class BS_API_PLUGIN well_iterator : public std::iterator <std::forward_iterator_tag, value_type_t>
//    {
//    public:
//
//      typedef well_iterator <list_t, iterator_t, value_type_t>  this_t;
//      typedef value_type_t																			value_type;
//
//    public:
//      well_iterator (const list_t &list, const iterator_t &it);
//
//      bool operator == (const this_t &rhs) const;
//      bool operator != (const this_t &rhs) const;
//
//      const this_t &operator = (const this_t &src);
//
//      value_type operator * () const;
//
//      this_t& operator ++ ();
//      this_t operator ++ (int);
//
//    private:
//      list_t			list_;
//      iterator_t	it_;
//    };
//
//
//
//} // namespace blue_sky
//
//#endif  // #ifndef BS_WELL_ITERATOR_H_
//

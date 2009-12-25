/**
 *       \file  exp_temp_mx.h
 *      \brief  Matrix arithmetics on expression templates
 *     \author  Sergey Miryanov (sergey-miryanov), sergey.miryanov@gmail.com
 *       \date  20.10.2008
 *  \copyright  This source code is released under the terms of 
 *              the BSD License. See LICENSE for more details.
 *       \todo  Obsolete, should be removed or redisigned
 * */
#ifndef BS_EXPRESSION_TEMPLATES_MATRIX_ARITHMETIC_H_
#define BS_EXPRESSION_TEMPLATES_MATRIX_ARITHMETIC_H_

namespace blue_sky
  {

  struct ref_trait {};
  struct ref_unsafe_trait {};
  struct const_trait {};
  struct const_unsafe_trait {};

  template <typename array_t_, typename trait_t = const_trait>
  struct mx_holder
    {
      typedef array_t_  array_t;
      typedef typename array_t::value_type value_type;
      typedef mx_holder <array_t, trait_t> this_t;

      mx_holder (const array_t &array_, size_t size)
          : array_ (array_)
          , block_size_ (size * size)
      {
        BS_ASSERT (block_size_ <= 9) (block_size_);
      }

      value_type
      operator [] (size_t i) const
        {
          return array_[i];
        }

      size_t
      block_size () const
        {
          return block_size_;
        }

      const array_t &array_;
      size_t          block_size_;
    };
  template <typename array_t_>
  struct mx_holder <array_t_, ref_trait>
    {

      typedef array_t_ array_t;
      typedef typename array_t::value_type value_type;
      typedef mx_holder <array_t, ref_trait> this_t;

      mx_holder (array_t &array_, size_t idx, size_t size)
          : array_ (array_)
          , block_size_ (size * size)
          , idx_ (idx * block_size_)
      {
        BS_ASSERT (block_size_ <= 9) (block_size_);
      }

      template <typename op_t>
      this_t &
      operator= (const op_t &op)
      {
        for (size_t i = 0, cnt = block_size_; i < cnt; ++i)
          {
            array_[i + idx_] = op[i];
          }

        return *this;
      }

      template <typename op_t>
      this_t &
      operator+= (const op_t &op)
      {
        for (size_t i = 0, cnt = block_size_; i < cnt; ++i)
          {
            array_[i + idx_] += op[i];
          }

        return *this;
      }

      value_type &
      operator [] (size_t index) 
      {
        if (idx_ + index >= array_.size ())
          throw bs_exception ("mx_holder[]", "index out of range");

        return array_[idx_ + index];
      }

      array_t   &array_;
      size_t      block_size_;
      size_t      idx_;
    };
  template <typename array_t_, typename trait_t = const_trait>
  struct clmn_holder
    {

      typedef array_t_  array_t;
      typedef typename array_t::value_type value_type;
      typedef clmn_holder <array_t, trait_t> this_t;

      clmn_holder (const array_t &array_, size_t size)
          : array_ (array_)
          , block_size_ (size)
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      value_type
      operator [] (size_t i) const
        {
          return array_[i % block_size_];
        }

      size_t
      block_size () const
        {
          return block_size_;
        }

      const array_t  &array_;
      size_t           block_size_;
    };
  template <typename array_t_>
  struct clmn_holder <array_t_, const_unsafe_trait>
    {
      typedef array_t_                                    array_t;
      typedef typename array_t::value_type                value_type;
      typedef clmn_holder <array_t, const_unsafe_trait>   this_t;

      clmn_holder (const array_t &array_, size_t index, size_t size)
          : array_ (array_)
          , index_ (index * size)
          , block_size_ (size)
          , unsafe_array_ (&array_[index_])
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      value_type
      operator [] (size_t i) const
        {
          return unsafe_array_[i % block_size_];
        }

      size_t
      block_size () const
        {
          return block_size_;
        }

      const array_t       &array_;
      size_t              index_;
      size_t              block_size_;
      const value_type    *unsafe_array_;
    };
  template <typename array_t_>
  struct clmn_holder <array_t_, ref_trait>
    {
      typedef array_t_                          array_t;
      typedef typename array_t::value_type      value_type;
      typedef clmn_holder <array_t, ref_trait>  this_t;

      clmn_holder (array_t &array, size_t idx, size_t size)
          : array_ (array)
          , idx_ (idx * size)
          , block_size_ (size)
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      template <typename op_t>
      this_t &
      operator= (const op_t &op)
      {
        for (size_t i = 0, cnt = block_size_; i < cnt; ++i)
          {
            array_[idx_ + i] = op[i];
          }

        return *this;
      }

      array_t &array_;
      size_t    idx_;
      size_t    block_size_;
    };

  template <typename array_type, typename trait_t = const_trait>
  struct line_holder
    {
      typedef array_type                      array_t;
      typedef typename array_t::value_type    value_type;
      typedef line_holder <array_t, trait_t>  this_t;

      line_holder (const array_t &array_, size_t size)
          : array_ (array_)
          , block_size_ (size)
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      value_type
      operator [] (size_t i) const
        {
          return array_[i % block_size_];
        }
      size_t
      block_size () const
        {
          return block_size_;
        }

      const array_t &array_;
      size_t          block_size_;
    };

  template <typename array_type>
  struct line_holder <array_type, const_unsafe_trait>
    {
      typedef array_type                      array_t;
      typedef typename array_t::value_type    value_type;
      typedef line_holder <array_t, const_unsafe_trait>  this_t;

      line_holder (const array_t &array_, size_t index, size_t size)
          : array_ (array_)
          , index_ (index)
          , block_size_ (size)
          , unsafe_array_ (&array_[index_ * block_size_])
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      value_type
      operator [] (size_t i) const
        {
          return unsafe_array_[i % block_size_];
        }
      size_t
      block_size () const
        {
          return block_size_;
        }

      const array_t       &array_;
      size_t              block_size_;
      size_t              index_;
      const value_type    *unsafe_array_;
    };

  template <typename array_type>
  struct line_holder <array_type, ref_unsafe_trait>
    {
      typedef array_type                      array_t;
      typedef typename array_t::value_type    value_type;
      typedef line_holder <array_t, ref_unsafe_trait>  this_t;

      line_holder (const array_t &array_, size_t index, size_t size)
          : array_ (array_)
          , index_ (index)
          , block_size_ (size)
          , unsafe_array_ (&array_[index_ * block_size_])
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      template <typename op_t>
      this_t &
      operator= (const op_t &op)
      {
        for (size_t i = 0, cnt = block_size_; i < cnt; ++i)
          {
            unsafe_array_[i] = op [i];
          }

        return *this;
      }

      const array_t       &array_;
      size_t              block_size_;
      size_t              index_;
      const value_type    *unsafe_array_;
    };

  template <typename array_type>
  struct line_holder <array_type, ref_trait>
    {
      typedef array_type                        array_t;
      typedef typename array_t::value_type      value_type;
      typedef line_holder <array_t, ref_trait>  this_t;

      line_holder (array_t &array_, size_t size)
          : array_ (array_)
          , block_size_ (size)
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      template <typename op_t>
      this_t &
      operator= (const op_t &op)
      {
        for (size_t i = 0, cnt = block_size_; i < cnt; ++i)
          {
            array_[i] = op[i];
          }

        return *this;
      }

      array_t             &array_;
      size_t              block_size_;
    };

  template <typename mx_op_t>
  struct mx_op
    {
      typedef typename mx_op_t::value_type  value_type;
      typedef mx_op_t                       op_t;

      mx_op (const mx_op_t &mx_op_)
          : mx_op_ (mx_op_)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return mx_op_[i];
        }

      size_t
      block_size () const
        {
          return mx_op_.block_size ();
        }

      mx_op_t mx_op_;
    };

  template <typename clmn_op_t>
  struct clmn_op
    {
      typedef typename clmn_op_t::value_type  value_type;
      typedef clmn_op_t                       op_t;

      clmn_op (const clmn_op_t &clmn_op_)
          : clmn_op_ (clmn_op_)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return clmn_op_ [i];
        }

      size_t
      block_size () const
        {
          return clmn_op_.block_size ();
        }

      clmn_op_t clmn_op_;
    };

  template <typename line_op_t>
  struct line_op
    {
      typedef typename line_op_t::value_type value_type;
      typedef line_op_t                       op_t;

      line_op (const line_op_t &line_op)
          : line_op_ (line_op)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return line_op_[i];
        }

      size_t
      block_size () const
        {
          return line_op_.block_size ();
        }

      line_op_t line_op_;
    };

  template <typename scalar_op_t>
  struct scalar_op
    {
      typedef scalar_op_t                   op_t;
      typedef typename op_t::value_type     value_type;

      scalar_op (const op_t &op)
          : op_ (op)
      {
      }

      operator value_type () const
        {
          return op_;
        }

      op_t op_;
    };

  ///
  template <typename lhs_t, typename rhs_t>
  struct mx_mul
    {
      typedef typename lhs_t::value_type value_type;

      mx_mul (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      template <typename array_t>
      mx_mul (const lhs_t &lhs, const array_t &array)
          : lhs_ (lhs)
          , rhs_ (rhs_t (array))
      {
      }

      template <typename array_t>
      mx_mul (const array_t &array, const rhs_t &rhs)
          : lhs_ (lhs_t (array))
          , rhs_ (rhs)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          size_t bs = lhs_.block_size ();
          size_t ri = i / bs;
          size_t ci = i % bs;

          value_type res = 0;
          for (size_t j = 0; j < bs; ++j)
            {
              res += lhs_[ri * bs + j] * rhs_[j * bs + ci];
            }

          return res;
        }

      size_t
      block_size () const
        {
          return lhs_.block_size ();
        }

      lhs_t lhs_;
      rhs_t rhs_;
    };

  template <typename lhs_t, typename rhs_t>
  struct mx_sub
    {
      typedef typename lhs_t::value_type value_type;

      mx_sub (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      template <typename array_t>
      mx_sub (const lhs_t &lhs, const array_t &array)
          : lhs_ (lhs)
          , rhs_ (rhs_t (array))
      {
      }

      template <typename array_t>
      mx_sub (const array_t &array, const rhs_t &rhs)
          : lhs_ (lhs_t (array))
          , rhs_ (rhs)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return lhs_[i] - rhs_[i];
        }

      size_t
      block_size () const
        {
          return lhs_.block_size ();
        }

      lhs_t lhs_;
      rhs_t rhs_;
    };

  template <typename lhs_t, typename rhs_t>
  struct mx_clmn_mul
    {
      typedef typename lhs_t::value_type value_type;

      mx_clmn_mul (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return lhs_[i] * rhs_[i / lhs_.block_size ()];
        }

      size_t
      block_size () const
        {
          return lhs_.block_size ();
        }

      lhs_t lhs_;
      rhs_t rhs_;
    };

  template <typename lhs_t, typename rhs_t>
  struct clmn_line_mul
    {
      typedef typename lhs_t::value_type value_type;

      clmn_line_mul (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          size_t bs = lhs_.block_size ();
          size_t ri = i / bs;
          size_t ci = i % bs;

          return lhs_[ri] * rhs_[ci];
        }

      size_t
      block_size () const
        {
          return lhs_.block_size ();
        }

      const lhs_t &lhs_;
      const rhs_t &rhs_;
    };
  template <typename lhs_t, typename rhs_t>
  struct line_clmn_mul
    {
      typedef typename lhs_t::value_type value_type;

      line_clmn_mul (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      operator value_type () const
        {
          value_type res = 0;
          for (size_t i = 0, cnt = lhs_.block_size (); i < cnt; ++i)
            {
              res += lhs_[i] * rhs_[i];
            }

          return res;
        }

      lhs_t lhs_;
      rhs_t rhs_;
    };


  template <typename lhs_t, typename rhs_t>
  struct clmn_sub
    {
      typedef typename lhs_t::value_type value_type;

      clmn_sub (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return lhs_[i] - rhs_[i];
        }

      size_t
      block_size () const
        {
          return lhs_.block_size ();
        }

      lhs_t lhs_;
      rhs_t rhs_;
    };

  template <typename array_t_>
  struct mx_scalar_mul
    {
      typedef array_t_ array_t;
      typedef typename array_t::value_type  value_type;
      typedef mx_scalar_mul <array_t>  this_t;

      mx_scalar_mul (const array_t &lhs, value_type rhs, size_t size)
          : lhs_ (lhs)
          , rhs_ (rhs)
          , block_size_ (size * size)
      {
        BS_ASSERT (block_size_ <= 9) (block_size_);
      }

      value_type
      operator [] (size_t i) const
        {
          return lhs_[i] * rhs_;
        }

      size_t
      block_size () const
        {
          return block_size_;
        }

      array_t       lhs_;
      value_type    rhs_;
      size_t          block_size_;
    };
  template <typename array_type>
  struct clmn_scalar_mul
    {
      typedef array_type                    array_t;
      typedef typename array_t::value_type  value_type;
      typedef clmn_scalar_mul <array_t>     this_t;

      clmn_scalar_mul (const array_t &lhs, value_type rhs, size_t size)
          : lhs_ (lhs)
          , rhs_ (rhs)
          , block_size_ (size)
      {
        BS_ASSERT (block_size_ <= 3) (block_size_);
      }

      value_type
      operator [] (size_t i) const
        {
          BS_ASSERT (i < block_size_);
          return lhs_[i] * rhs_;
        }

      size_t
      block_size () const
        {
          return block_size_;
        }

      const array_t &lhs_;
      value_type    rhs_;
      size_t          block_size_;
    };

  template <typename lhs_t, typename rhs_t>
  struct line_mul
    {
      typedef typename lhs_t::value_type value_type;

      line_mul (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return lhs_[i] * rhs_[i];
        }

      size_t
      block_size () const
        {
          return lhs_.block_size ();
        }

      lhs_t lhs_;
      rhs_t rhs_;
    };
  template <typename lhs_t, typename rhs_t>
  struct line_sub
    {
      typedef typename lhs_t::value_type value_type;

      line_sub (const lhs_t &lhs, const rhs_t &rhs)
          : lhs_ (lhs)
          , rhs_ (rhs)
      {
      }

      value_type
      operator [] (size_t i) const
        {
          return lhs_[i] - rhs_[i];
        }

      size_t
      block_size () const
        {
          return lhs_.block_size ();
        }

      lhs_t lhs_;
      rhs_t rhs_;
    };


  template <typename array_t_>
  struct mx_inverse
    {
      typedef array_t_                        array_t;
      typedef typename array_t::value_type    value_type;
      typedef mx_inverse <array_t>            this_t;

      mx_inverse (const array_t &mx, size_t size)
          : mx_ (mx)
          , block_size_ (size * size)
      {
        BS_ASSERT (block_size_ <= 9) (block_size_);
        calc_inverse_matrix (block_size_, mx_);
      }

      value_type
      operator [] (size_t i) const
        {
          return mx_[i];
        }

      size_t
      block_size () const
        {
          return block_size_;
        }

      array_t       mx_;
      size_t        block_size_;
    };

  ///
  template <typename lhs_t, typename rhs_t>
  mx_op <mx_mul <lhs_t, rhs_t> >
  operator* (const mx_op <lhs_t> &lhs, const mx_op <rhs_t> &rhs)
  {
    return mx_op <mx_mul <lhs_t, rhs_t> > (mx_mul <lhs_t, rhs_t> (lhs.mx_op_, rhs.mx_op_));
  }
  template <typename lhs_t, typename array_t>
  mx_op <mx_mul <lhs_t, mx_holder <array_t> > >
  operator* (const mx_op <lhs_t> &lhs, const mx_holder <array_t> &rhs)
  {
    return mx_op <mx_mul <lhs_t, mx_holder <array_t> > > (mx_mul <lhs_t, mx_holder <array_t> > (lhs.mx_op_, rhs));
  }
  template <typename array_t, typename rhs_t>
  mx_op <mx_mul <mx_holder <array_t>, rhs_t> >
  operator* (const mx_holder <array_t> &lhs, const mx_op <rhs_t> &rhs)
  {
    return mx_op <mx_mul <mx_holder <array_t>, rhs_t> > (mx_mul <mx_holder <array_t>, rhs_t> (lhs, rhs.mx_op_));
  }
  template <typename array_t, typename rhs_t>
  mx_op <mx_sub <mx_holder <array_t>, rhs_t> >
  operator- (const mx_holder <array_t> &lhs, const mx_op <rhs_t> &rhs)
  {
    return mx_op <mx_sub <mx_holder <array_t>, rhs_t> > (mx_sub <mx_holder <array_t>, rhs_t> (lhs, rhs.mx_op_));
  }
  template <typename lhs_t, typename array_t>
  mx_op <mx_sub <lhs_t, mx_holder <array_t> > >
  operator- (const mx_op <lhs_t> &lhs, const mx_holder <array_t> &rhs)
  {
    return mx_op <mx_sub <lhs_t, mx_holder <array_t> > > (mx_sub <lhs_t, mx_holder <array_t> > (lhs.mx_op_, rhs));
  }
  template <typename lhs_t, typename rhs_t>
  mx_op <mx_sub <lhs_t, rhs_t> >
  operator- (const mx_op <lhs_t> &lhs, const mx_op <rhs_t> &rhs)
  {
    return mx_op <mx_sub <lhs_t, rhs_t> > (mx_sub <lhs_t, rhs_t> (lhs.mx_op_, rhs.mx_op_));
  }

  template <typename lhs_t, typename rhs_t>
  mx_op <clmn_line_mul <clmn_holder <lhs_t>, line_holder <rhs_t> > >
  operator* (const clmn_holder <lhs_t> &lhs, const line_holder <rhs_t> &rhs)
  {
    return mx_op <clmn_line_mul <clmn_holder <lhs_t>, line_holder <rhs_t> > > (clmn_line_mul <clmn_holder <lhs_t>, line_holder <rhs_t> > (lhs, rhs));
  }
  template <typename lhs_t, typename rhs_t>
  mx_op <clmn_line_mul <clmn_holder <lhs_t>, line_holder <rhs_t, const_unsafe_trait> > >
  operator* (const clmn_holder <lhs_t> &lhs, const line_holder <rhs_t, const_unsafe_trait> &rhs)
  {
    return mx_op <clmn_line_mul <clmn_holder <lhs_t>, line_holder <rhs_t, const_unsafe_trait> > > (clmn_line_mul <clmn_holder <lhs_t>, line_holder <rhs_t, const_unsafe_trait> > (lhs, rhs));
  }
  template <typename lhs_t, typename rhs_t>
  mx_op <clmn_line_mul <lhs_t, line_holder <rhs_t> > >
  operator* (const clmn_op <lhs_t> &lhs, line_holder <rhs_t> &rhs)
  {
    return mx_op <clmn_line_mul <lhs_t, line_holder <rhs_t> > > (clmn_line_mul <lhs_t, line_holder <rhs_t> > (lhs.clmn_op_, rhs));
  }
  template <typename lhs_t, typename rhs_t>
  line_op <line_mul <line_holder <lhs_t>, line_holder <rhs_t, const_unsafe_trait> > >
  operator* (const line_holder <lhs_t> &lhs, const line_holder <rhs_t, const_unsafe_trait> &rhs)
  {
    return line_op <line_mul <line_holder <lhs_t>, line_holder <rhs_t, const_unsafe_trait> > > (line_mul <line_holder <lhs_t>, line_holder <rhs_t, const_unsafe_trait> > (lhs, rhs));
  }
  template <typename lhs_t, typename lhs_trait_t, typename rhs_t, typename rhs_trait_t>
  scalar_op <line_clmn_mul <line_holder <lhs_t, lhs_trait_t>, clmn_holder <rhs_t, rhs_trait_t> > >
  operator* (const line_holder <lhs_t, lhs_trait_t> &lhs, const clmn_holder <rhs_t, rhs_trait_t> &rhs)
  {
    return scalar_op <line_clmn_mul <line_holder <lhs_t, lhs_trait_t>, clmn_holder <rhs_t, rhs_trait_t> > > (line_clmn_mul <line_holder <lhs_t, lhs_trait_t>, clmn_holder <rhs_t, rhs_trait_t> > (lhs, rhs));
  }
  template <typename lhs_t, typename rhs_t, typename rhs_trait_t>
  scalar_op <line_clmn_mul <lhs_t, clmn_holder <rhs_t, rhs_trait_t> > >
  operator* (const line_op <lhs_t> &lhs, const clmn_holder <rhs_t, rhs_trait_t> &rhs)
  {
    return scalar_op <line_clmn_mul <lhs_t, clmn_holder <rhs_t, rhs_trait_t> > > (line_clmn_mul <lhs_t, clmn_holder <rhs_t, rhs_trait_t> > (lhs.line_op_, rhs));
  }
  template <typename lhs_t, typename rhs_t>
  line_op <line_mul <lhs_t, rhs_t> >
  operator* (const line_op <lhs_t> &lhs, const line_op <rhs_t> &rhs)
  {
    return line_op <line_mul <lhs_t, rhs_t> > (line_mul <lhs_t, rhs_t> (lhs.line_op_, rhs.line_op_));
  }
  template <typename lhs_t, typename rhs_t, typename rhs_trait_t>
  line_op <line_mul <lhs_t, line_holder <rhs_t, rhs_trait_t> > >
  operator* (const line_op <lhs_t> &lhs, const line_holder <rhs_t, rhs_trait_t> &rhs)
  {
    return line_op <line_mul <lhs_t, line_holder <rhs_t, rhs_trait_t> > > (line_mul <lhs_t, line_holder <rhs_t, rhs_trait_t> > (lhs.line_op_, rhs));
  }
  template <typename lhs_t, typename lhs_trait_t, typename rhs_t>
  line_op <line_sub <line_holder <lhs_t, lhs_trait_t>, rhs_t> >
  operator- (const line_holder <lhs_t, lhs_trait_t> &lhs, const line_op <rhs_t> &rhs)
  {
    return line_op <line_sub <line_holder <lhs_t, lhs_trait_t>, rhs_t> > (line_sub <line_holder <lhs_t, lhs_trait_t>, rhs_t> (lhs, rhs.line_op_));
  }
  template <typename lhs_t, typename rhs_t>
  line_op <line_sub <lhs_t, rhs_t> >
  operator- (const line_op <lhs_t> &lhs, const line_op <rhs_t> &rhs)
  {
    return line_op <line_sub <lhs_t, rhs_t> > (line_sub <lhs_t, rhs_t> (lhs.line_op_, rhs.line_op_));
  }
  template <typename array_t>
  clmn_op <mx_clmn_mul <mx_holder <array_t>, clmn_holder <array_t> > >
  operator* (const mx_holder <array_t> &mx, const clmn_holder<array_t> &column)
  {
    return clmn_op <mx_clmn_mul <mx_holder <array_t>, clmn_holder <array_t> > > (mx_clmn_mul <mx_holder <array_t>, clmn_holder <array_t> > (mx, column));
  }
  template <typename lhs_t, typename array_t>
  clmn_op <mx_clmn_mul <lhs_t, clmn_holder <array_t> > >
  operator* (const mx_op <lhs_t> &lhs, const clmn_holder <array_t> &column)
  {
    return clmn_op <mx_clmn_mul <lhs_t, clmn_holder <array_t> > > (mx_clmn_mul <lhs_t, clmn_holder <array_t> > (lhs.mx_op_, column));
  }
  //
  template <typename lhs_t, typename rhs_t>
  clmn_op <clmn_sub <lhs_t, rhs_t> >
  operator- (const clmn_op <lhs_t> &lhs, const clmn_op <rhs_t> &rhs)
  {
    return clmn_op <clmn_sub <lhs_t, rhs_t> > (clmn_sub <lhs_t, rhs_t> (lhs.clmn_op_, rhs.clmn_op_));
  }
  template <typename lhs_t, typename lhs_trait_t, typename rhs_t, typename rhs_trait_t>
  clmn_op <clmn_sub <clmn_holder <lhs_t, lhs_trait_t>, clmn_holder <rhs_t, rhs_trait_t> > >
  operator- (const clmn_holder <lhs_t, lhs_trait_t> &lhs, const clmn_holder <rhs_t, rhs_trait_t> &rhs)
  {
    return clmn_op <clmn_sub <clmn_holder <lhs_t, lhs_trait_t>, clmn_holder <rhs_t, rhs_trait_t> > > (clmn_sub <clmn_holder <lhs_t, lhs_trait_t>, clmn_holder <rhs_t, rhs_trait_t> > (lhs, rhs));
  }
  template <typename lhs_t, typename lhs_trait_t, typename rhs_t>
  clmn_op <clmn_sub <clmn_holder <lhs_t, lhs_trait_t>, rhs_t> >
  operator- (const clmn_holder <lhs_t, lhs_trait_t> &lhs, const clmn_op <rhs_t> &rhs)
  {
    return clmn_op <clmn_sub <clmn_holder <lhs_t, lhs_trait_t>, rhs_t> > (clmn_sub <clmn_holder <lhs_t, lhs_trait_t>, rhs_t> (lhs, rhs.clmn_op_));
  }
  //
  template <typename mx_t>
  mx_op <mx_scalar_mul <typename mx_t::array_t> >
  operator* (const mx_t &mx, typename mx_t::value_type scalar)
  {
    return mx_op <mx_scalar_mul <typename mx_t::array_t> > (mx_scalar_mul <typename mx_t::array_t> (mx.array_, scalar, mx.block_size_));
  }
  template <typename lhs_t>
  mx_op <mx_scalar_mul <lhs_t> >
  operator* (const mx_op <lhs_t> &mx, typename lhs_t::value_type scalar)
  {
    return mx_op <mx_scalar_mul <lhs_t> > (mx_scalar_mul <lhs_t> (mx.mx_op_, scalar, (size_t)mx.block_size ()));
  }
  template <typename array_t>
  clmn_op <clmn_scalar_mul <array_t> >
  operator* (const clmn_holder <array_t> &clmn, typename array_t::value_type scalar)
  {
    return clmn_op <clmn_scalar_mul <array_t> > (clmn_scalar_mul <array_t> (clmn.array_, scalar, clmn.block_size_));
  }
  template <typename lhs_t>
  clmn_op <clmn_scalar_mul <lhs_t> >
  operator* (const clmn_op <lhs_t> &clmn, typename lhs_t::value_type scalar)
  {
    return clmn_op <clmn_scalar_mul <lhs_t> > (clmn_scalar_mul <lhs_t> (clmn.clmn_op_, scalar, (size_t)clmn.block_size ()));
  }
  template <typename array_t>
  clmn_op <clmn_scalar_mul <array_t> >
  operator* (typename array_t::value_type scalar, const clmn_holder <array_t> &clmn)
  {
    return clmn_op <clmn_scalar_mul <array_t> > (clmn_scalar_mul <array_t> (clmn.array_, scalar, clmn.block_size_));
  }
  template <typename lhs_t, typename rhs_t, typename rhs_trait_t>
  line_op <clmn_scalar_mul <rhs_t> >
  operator* (lhs_t scalar, const line_holder <rhs_t, rhs_trait_t> &rhs)
  {
    return line_op <clmn_scalar_mul <rhs_t> > (clmn_scalar_mul <rhs_t> (rhs.array_, scalar, rhs.block_size ()));
  }
  template <typename lhs_t, typename lhs_trait_t>
  line_op <clmn_scalar_mul <lhs_t> >
  operator* (const line_holder <lhs_t, lhs_trait_t> &lhs, typename lhs_t::value_type scalar)
  {
    return line_op <clmn_scalar_mul <lhs_t> > (clmn_scalar_mul <lhs_t> (lhs.array_, scalar, lhs.block_size ()));
  }
  template <typename item_t, typename rhs_t>
  item_t
  operator- (item_t lhs, const scalar_op <rhs_t> &rhs)
  {
    return lhs - static_cast <item_t> (rhs);
  }
  //
  template <typename mx_t>
  mx_op <mx_inverse <typename mx_t::array_t> >
  inverse (const mx_t &mx)
  {
    return mx_op <mx_inverse <typename mx_t::array_t> > (mx_inverse <typename mx_t::array_t> (mx.array_, mx.block_size_));
  }

} // namespace blue_sky


#endif  // #ifndef BS_EXPRESSION_TEMPLATES_MATRIX_ARITHEMTIC_H_


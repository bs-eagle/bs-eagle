/**
 * */
#include "bs_hdf5_storage_v2.h"
#include "bos_report.h"
#include "bs_kernel_tools.h"

#include "hdf5_hid_holder.hpp"
#include "hdf5_group_impl.hpp"
#include "hdf5_type_to_hid.hpp"
#include "hdf5_functions.h"

#include <boost/tokenizer.hpp>

namespace blue_sky {

  namespace detail 
  {
    hid_holder <hdf5_property_type>
    create_dataset_property ()
    {
      hid_t plist = H5Pcreate (H5P_DATASET_CREATE);
      if (plist < 0)
        {
          bs_throw_exception ("Can't create property for create dataset");
        }

      return hid_holder <hdf5_property_type> (plist);
    }

    hid_holder <hdf5_property_type>
    raw_transfer_property ()
    {
      hid_t plist = H5Pcreate (H5P_DATASET_XFER);
      if (plist < 0)
        {
          bs_throw_exception ("Can't create property for raw data transfer");
        }

      return plist;
    }
  }


  struct hdf5_storage_v2
  {
  public:

    static hdf5_storage_v2 *
    instance ();

    hdf5_storage_v2 ();
    ~hdf5_storage_v2 ();

    struct impl;
    impl *impl_;

    friend struct hdf5_file;
    friend struct hdf5_group_v2;
    friend struct impl;
  };


  typedef hdf5::private_::hdf5_buffer__ hdf5_buffer_t;

  struct hdf5_storage_v2::impl 
  {
    typedef std::map <std::string, hid_t>   file_id_map_t;

    impl ()
    : grub_call_stack_ (false)
    , error_context_ ("")
    {
      // we use in hdf5_file int instead of hid_t so we should check both types on eq
      BOOST_STATIC_ASSERT (sizeof (hid_t) == sizeof (int));
    }

    ~impl ()
    {
      file_id_map_t::iterator it = file_map_.begin (), e = file_map_.end ();
      for (; it != e; ++it)
        {
          if (it->second >= 0)
            {
              H5Fclose (it->second);
            }
        }
    }

    void
    register_file (const std::string &name, hid_t id)
    {
      file_id_map_t::iterator it = file_map_.find (name);
      if (it == file_map_.end ())
        {
          file_map_.insert (std::make_pair (name, id));
        }
      else
        {
          if (it->second != -1)
            {
              BOSWARN (section::app_info, level::warning) 
                << boost::format ("Register already opened file %s, (old id: %d, new id: %d)")
                  % name % it->second % id
                << bs_end;
            }
          else
            {
              it->second = id;
            }
        }
    }

    hid_t
    get_file_id (const std::string &name) const
    {
      BS_ASSERT (name.size ());

      file_id_map_t::const_iterator it = file_map_.find (name);
      if (it == file_map_.end ())
        {
          return -1;
        }

      return it->second;
    }

    hid_t
    get_file_id (hdf5_file &file)
    {
      BS_ASSERT (file.file_name_.size ());

      if (file.file_id_ < 0)
        {
          file.file_id_ = get_file_id (file.file_name_);
        }

      if (file.file_id_ < 0)
        {
          file.file_id_ = H5Fcreate (file.file_name_.c_str (), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          if (file.file_id_ < 0)
            {
              bs_throw_exception (boost::format ("Can't open file %s") % file.file_name_);
            }

          hid_group_t group = H5Gcreate (file.file_id_, "/results", 0);
          if (group < 0)
            {
              bs_throw_exception (boost::format ("Can't create group '/results' in file %s") % file.file_name_);
            }

          register_file (file.file_name_, file.file_id_);
        }

      return file.file_id_;
    }

    void
    create_group_hierarchy (hid_t file_id, const std::string &group_name)
    {
      typedef boost::tokenizer <boost::char_separator <char> > tokenizer_t;
      boost::char_separator <char> s ("/");
      tokenizer_t t (group_name, s);
      std::string gn = "";
      for (tokenizer_t::iterator it = t.begin (); it != t.end (); ++it)
        {
          gn += "/" + *it;
          set_error_context (gn);
          if (!detail::is_object_exists (file_id, gn.c_str ()))
            {
              detail::hid_holder <detail::hdf5_group_type> group = H5Gcreate (file_id, gn.c_str (), H5P_DEFAULT);
              if (group < 0)
                {
                  bs_throw_exception (boost::format ("Can't create group %s") % gn);
                }
            }
        }
    }

    template <typename data_t>
    void
    write_to_hdf5 (hdf5_file &file, 
                   const std::string &group_name, 
                   const std::string &dataset_name, 
                   data_t &buffer)
    {
      BS_ASSERT (group_name.size ()) (dataset_name);

      hid_t file_id               = get_file_id (file);
      const int rank              = 2;
      hsize_t chunk_dimens[rank]  = {buffer.size, 1};
      hsize_t fspace_dimens[rank] = {buffer.size, 0};
      hsize_t mspace_dimens[rank] = {buffer.size, 1};
      hsize_t max_dimens[rank]    = {buffer.size, H5S_UNLIMITED};

      if (set_error_context (group_name) && !detail::is_object_exists (file_id, group_name.c_str ()))
        {
          create_group_hierarchy (file_id, group_name);
        }
      hid_group_t group = detail::open_group (file_id, group_name.c_str ());
      BS_ASSERT (group >= 0) (group_name) (dataset_name);

      if (set_error_context (dataset_name) && !detail::is_object_exists (group, dataset_name.c_str ()))
        {
          hid_property_t plist = detail::create_dataset_property ();
          if (H5Pset_chunk (plist, rank, chunk_dimens) < 0)
            {
              bs_throw_exception (boost::format ("Can't set chunk for dataset %s in group %s") % dataset_name % group_name);
            }

          hid_dspace_t fspace = H5Screate_simple (rank, fspace_dimens, max_dimens);
          if (fspace < 0)
            {
              bs_throw_exception (boost::format ("Can't create simple dataspace for dataset %s in group %s") % dataset_name % group_name);
            }

          hid_dset_t ds_value = H5Dcreate (group, dataset_name.c_str (), get_hdf5_type (buffer.type), fspace, plist);
          if (ds_value < 0)
            {
              bs_throw_exception (boost::format ("Can't create dataset %s in group %s") % dataset_name % group_name);
            }
        }

      hid_dset_t dataset = detail::open_dataset (group, dataset_name.c_str ());
      hid_dspace_t space = H5Dget_space (dataset);
      if (space < 0)
        {
          bs_throw_exception (boost::format ("Can't get space for dataset %s in group %s") % dataset_name % group_name);
        }

      hsize_t size_value[] = {0, 0};
      if (H5Sget_simple_extent_dims (space, size_value, NULL) < 0)
        {
          bs_throw_exception (boost::format ("Can't get dataspace dims for dataset %s in group %s") % dataset_name % group_name);
        }

      hsize_t size[] = {size_value[0], size_value[1] + 1};
      if (H5Dextend (dataset, size) < 0)
        {
          bs_throw_exception (boost::format ("Can't extend dataset %s in group %s") % dataset_name % group_name);
        }

      hid_dspace_t fspace = H5Dget_space (dataset);
      if (fspace < 0)
        {
          bs_throw_exception (boost::format ("Can't get dataspace for dataset %s in group %s") % dataset_name % group_name);
        }

      hsize_t offset[] = {0, size_value[1]};
      if (H5Sselect_hyperslab (fspace, H5S_SELECT_SET, offset, NULL, mspace_dimens, NULL) < 0)
        {
          bs_throw_exception (boost::format ("Can't select hyperslab for dataset %s in group %s") % dataset_name % group_name);
        }

      hid_dspace_t mspace = H5Screate_simple (rank, mspace_dimens, max_dimens);
      if (mspace < 0)
        {
          bs_throw_exception (boost::format ("Can't create simple dataspace for dataset %s in group %s") % dataset_name % group_name);
        }

      hid_t status = H5Dwrite (dataset, get_hdf5_type (buffer.type), mspace, fspace, detail::raw_transfer_property (), buffer.data ());
      if (status < 0)
        {
          bs_throw_exception (boost::format ("Can't write data (dataset %s, group %s)") % dataset_name % group_name);
        }
    }

    void
    write_string_to_hdf5 (hdf5_file &file,
                          const std::string &group_name,
                          const std::string &dataset_name,
                          const std::string &value)
    {
      BS_ASSERT (group_name.size ()) (dataset_name);

      hid_t file_id = get_file_id (file);

      if (set_error_context (group_name) && !detail::is_object_exists (file_id, group_name.c_str ()))
        {
          create_group_hierarchy (file_id, group_name);
        }
      hid_group_t group = detail::open_group (file_id, group_name.c_str ());
      BS_ASSERT (group >= 0) (group_name) (dataset_name);

      if (set_error_context (dataset_name) && !detail::is_object_exists (group, dataset_name.c_str ()))
        {
          hid_dspace_t dspace = H5Screate (H5S_SCALAR);
          if (dspace < 0)
            {
              bs_throw_exception (boost::format ("Can't create dataspace for dataset %s in group %s") % dataset_name % group_name);
            }

          hid_dtype_t dtype = H5Tcopy (H5T_C_S1);
          if (dtype < 0)
            {
              bs_throw_exception (boost::format ("Can't create datatype for dataset %s in group %s") % dataset_name % group_name);
            }

          if (H5Tset_size (dtype, value.length ()) < 0)
            {
              bs_throw_exception (boost::format ("Can't set size of datatype for dataset %s in group %s") % dataset_name % group_name);
            }

          hid_dset_t dset = H5Dcreate (group, dataset_name.c_str (), dtype, dspace, H5P_DEFAULT);
          if (dset < 0)
            {
              bs_throw_exception (boost::format ("Can't create dataset %s in group %s") % dataset_name % group_name);
            }

          hid_dspace_t fspace = H5Dget_space (dset);
          if (fspace < 0)
            {
              bs_throw_exception (boost::format ("Can't get space for dataset %s in group %s") % dataset_name % group_name);
            }

          hid_t status = H5Dwrite (dset, dtype, dspace, fspace, detail::raw_transfer_property (), value.c_str ());
          if (status < 0)
            {
              bs_throw_exception (boost::format ("Can't write data (dataset %s, group %s)") % dataset_name % group_name);
            }
        }
    }

    static herr_t
    hdf5_walk_handler (unsigned n, const H5E_error2_t *err, void *p)
    {
      //BOSOUT (section::save_data, level::critical) 
      //  << "walk: " << n
      //  << bs_end;
      

      size_t class_name_size = H5Eget_class_name (err->cls_id, 0, 0);
      size_t major_size      = H5Eget_msg (err->maj_num, 0, 0, 0);
      size_t minor_size      = H5Eget_msg (err->min_num, 0, 0, 0);

      if (class_name_size <= 0 || major_size <= 0 || minor_size <= 0)
        {
          bs_throw_exception (boost::format ("class_name_size: %d, major_size: %d, minor_size: %d") 
            % class_name_size % major_size % minor_size);
        }

      shared_vector <char> class_name (2 * class_name_size, 0);
      shared_vector <char> major_name (2 * major_size, 0);
      shared_vector <char> minor_name (2 * minor_size, 0);

      if (H5Eget_class_name (err->cls_id, &class_name[0], 2 * class_name_size - 1) <= 0)
        {
          bs_throw_exception ("Can't get class name");
        }
      if (H5Eget_msg (err->maj_num, NULL, &major_name[0], 2 * major_size - 1) <= 0)
        {
          bs_throw_exception (boost::format ("Can't get major message (class: %s)") % &class_name[0]);
        }
      if (H5Eget_msg (err->min_num, NULL, &minor_name[0], 2 * minor_size - 1) <= 0)
        {
          bs_throw_exception (boost::format ("Can't get minor message (class: %s, major: %s)") 
            % &class_name[0] % &major_name[0]);
        }

      std::string *str = static_cast <std::string *> (p);


      (*str) += (boost::format ("#%u: %s in %s(): line %u\n\t\tclass: %s\n\t\tmajor: %s\n\t\tminor: %s\n") 
            % n % err->file_name % err->func_name % err->line
            % &class_name[0]
            % &major_name[0]
            % &minor_name[0]).str ();

      return 0;
    }

    static herr_t 
    hdf5_error_handler (hid_t, void *p)
    {
      impl *impl_ = static_cast <impl *> (p);

      hid_t stack_id = H5Eget_current_stack ();
      if (stack_id < 0)
        {
          bs_throw_exception ("Can't get current stack");
        }

      std::string hdf5_stack_walk;
      hid_t status = H5Ewalk2 (stack_id, H5E_WALK_DOWNWARD, impl::hdf5_walk_handler, &hdf5_stack_walk);
      if (status < 0)
        {
          H5E_type_t type;
          size_t s = H5Eget_msg (status, &type, NULL, 0);
          if (!s)
            {
              bs_throw_exception ("Can't get error message");
            }

          shared_vector <char> buffer (s + 1, 0);
          s = H5Eget_msg (status, &type, &buffer[0], s);

          bs_throw_exception (boost::format ("Can't set walk2 handler: %s") % &buffer[0]);
        }

      H5Eclose_stack (stack_id);

      if (impl_->grub_call_stack_)
        {
          BOSOUT (section::save_data, level::critical)
            << boost::format ("HDF5 error occured (context: %s): \n%sCall stack which prevents to error: %s")
                % impl_->error_context_ % hdf5_stack_walk % kernel_tools::get_backtrace (128)
            << bs_end;
        }
      else
        {
          BOSOUT (section::save_data, level::critical)
            << boost::format ("HDF5 error occured (context: %s): \n%s")
                % impl_->error_context_ % hdf5_stack_walk 
            << bs_end;
        }

      bs_throw_exception ("HDF5 error v2");
      return 0;
    }

    bool 
    set_error_context (const std::string &s)
    {
      error_context_ = s;
      return true;
    }

    bool            grub_call_stack_;
    std::string     error_context_;
    file_id_map_t   file_map_;
  };

  hdf5_storage_v2::hdf5_storage_v2 ()
  : impl_ (new impl ())
  {
    H5Eset_auto2 (H5E_DEFAULT, impl::hdf5_error_handler, impl_);
  }

  hdf5_storage_v2::~hdf5_storage_v2 ()
  {
    delete impl_;
  }

  hdf5_storage_v2 *
  hdf5_storage_v2::instance ()
  {
    static hdf5_storage_v2 instance_;
    return &instance_;
  }

  hdf5_group_impl::hdf5_group_impl (bs_type_ctor_param)
  : file_ ("")
  , name_ ("")
  {
  }

  hdf5_group_impl::hdf5_group_impl (const hdf5_group_impl &src)
  : bs_refcounter (src)
  , file_ (src.file_)
  , name_ (src.name_)
  {
  }

  void
  hdf5_group_impl::init (hdf5_file const &file, std::string const &name)
  {
    file_ = file;
    name_ = name;
  }

  void
  hdf5_group_impl::write_buffer (const char *dataset, const hdf5_buffer_t &buffer)
  {
    if (buffer.size)
      hdf5_storage_v2::instance ()->impl_->write_to_hdf5 (file_, name_, dataset, buffer);
  }

  void
  hdf5_group_impl::write_pod (const char *dataset, const hdf5_pod_t &pod)
  {
    hdf5_storage_v2::instance ()->impl_->write_to_hdf5 (file_, name_, dataset, pod);
  }

  void
  hdf5_group_impl::write_struct (const char *dataset, const hdf5_struct_t &s)
  {
    hdf5_storage_v2::instance ()->impl_->write_to_hdf5 (file_, name_, dataset, s);
  }

  void
  hdf5_group_impl::write_string (const char *dataset, const std::string &value)
  {
    hdf5_storage_v2::instance ()->impl_->write_string_to_hdf5 (file_, name_, dataset, value);
  }

  BLUE_SKY_TYPE_STD_CREATE (hdf5_group_impl);
  BLUE_SKY_TYPE_STD_COPY (hdf5_group_impl);
  BLUE_SKY_TYPE_IMPL (hdf5_group_impl, objbase, "hdf5_group_impl", "hdf5_group_impl", "hdf5_group_impl");

} // namespace blue_sky


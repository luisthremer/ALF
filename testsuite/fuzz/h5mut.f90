! Corrupt an ALF confin_0.h5 for fuzz_inputs.py (built with ALF's own h5fc).
! Usage: h5mut FILE OP   OP = del_conf | del_seed | short_seed | conf_zero | conf_nan | conf_int
program h5mut
  use hdf5
  use h5lt
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  character(len=256) :: fn, op
  integer(hid_t) :: fid, did
  integer(hsize_t) :: dims(3), d1(1)
  integer(size_t) :: tsize
  integer :: err, rank, tclass, n
  real(kind(0.d0)), allocatable :: x(:)
  integer, allocatable :: iv(:)

  call get_command_argument(1, fn)
  call get_command_argument(2, op)
  call h5open_f(err)
  call h5fopen_f(trim(fn), H5F_ACC_RDWR_F, fid, err)
  if (err /= 0) stop "h5mut: cannot open file"
  dims = 1
  if (op(1:4) == "conf") then
     call h5ltget_dataset_ndims_f(fid, "configuration", rank, err)
     call h5ltget_dataset_info_f(fid, "configuration", dims, tclass, tsize, err)
     n = int(product(dims(1:rank)))
  endif
  select case (trim(op))
  case ("del_conf")
     call h5ldelete_f(fid, "configuration", err)
  case ("del_seed")
     call h5ldelete_f(fid, "seed", err)
  case ("short_seed")
     call h5ldelete_f(fid, "seed", err)
     d1 = 2
     call h5ltmake_dataset_int_f(fid, "seed", 1, d1, [12345, 678], err)
  case ("conf_zero", "conf_nan")
     allocate(x(n))
     x = 0.d0
     if (op == "conf_nan") x = ieee_value(1.d0, ieee_quiet_nan)
     call h5dopen_f(fid, "configuration", did, err)
     call h5dwrite_f(did, H5T_NATIVE_DOUBLE, x, dims(1:rank), err)
     call h5dclose_f(did, err)
  case ("conf_int")
     allocate(iv(n))
     iv = 1
     call h5ldelete_f(fid, "configuration", err)
     call h5ltmake_dataset_int_f(fid, "configuration", rank, dims(1:rank), iv, err)
  case default
     stop "h5mut: unknown op"
  end select
  if (err /= 0) stop "h5mut: HDF5 call failed"
  call h5fclose_f(fid, err)
  call h5close_f(err)
end program h5mut

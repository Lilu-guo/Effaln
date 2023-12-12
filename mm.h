#ifndef MM_H_
#define MM_H_
#define MM_READ(file, dest, sz) fread(dest, 1, sz, file)
#define MM_IS_IO_ERR(file_hd, ret, count) is_fread_err(file_hd, ret, count)
#endif

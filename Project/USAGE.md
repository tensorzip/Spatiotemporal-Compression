Operating system：linux x86

(BASIC: only employ time prediction model and last value prediction model )
Compress command:  ./compressor -c file.bin
Decompress command:  ./compressor -d bin_file.xor

Format:  
nz_num(unsigned int) t_num(unsigned int) 
nz_array_1 (double array)
nz_array_2 (double array)
...
nz_array_n (double array)

(Advanced version: employ time prediction model ,last value prediction model and stamping-based prediction model)
Developed on the basis of Xyce.

For the port of 3.2.2 to windows, these changes were made:

* All c files are cpp now.

* All the pthread methods were replaced with 'std::' methods.  This involved making a #define to do 
some logical task and defining it differently for pthread and std.  See mps.h.  
The program should still compile and run using pthread - you'll just have to provide access to 
the .h and .lib files.  My testing showed a weird fatal error when both pthread and 
std are in the same exe, so I had to do this conversion for my app which uses a lib version of mpsolve.

* ostream was being used to get a string equivalent of mpc variables, this is now done with 
a 'get_string' method, so ostream is no longer referenced either.

* mps_unistd.h was just copied from unistd.h for some types.  Now the program is completely windows.

* extra.cpp was added with some utility methods, but mostly it implements a 'tracked' mutex own/release 
so that no actual lock of a mutex is issued unless the thread does not already own it, and no actual 
unlock unless it does own it. pthread was doing nothing on a pthread_mutex_t except returning an error 
condition that is not checked. 'std::mutex' actually errors-out on such double-requests.  
The concept of a recursive mutex is not applicable here.

* malloc/free is now mps_malloc/mps_free, and new/delete is now mps_new_obj/mps_new_array_obj/
mps_delete_obj/mps_delete_array_obj.  structs are now obtained with new instead of malloc so that the 
std:mutex variables can be initialized.

* Variable-name improvements and standardization, such as 'ctx' everywhere for a ctx pointer, and such as 
'next_vertex' instead of just 'next'.

* _MPS_PRIVATE is always defined.

* a few changes to accomodate VS2022 needs, such as renaming the variable named 'free' to 'freeres'.

* the flex and bison code was regenerated using recent versions.  Not much change from flex, a lot of 
changes from bison. Oddly, lines 858 and 861 of tokenizer had to be changed to get a clean compile 
even before the regen.

* The phrase '// that may be a memory leak' appears in 4 places, as discovered with a leak-detector.  
They did not appear to be high-volume losses so were not investigated.

* Added 'mps_mutex_destroy(ctx->data_prec_max.value_mutex);' in data.cpp before line 272

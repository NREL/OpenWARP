/* This code is governed by the X11 license. See LICENSE for details. */
/* The lines below have little spacing to ease the line-length requirement */
#define log(level,format) if(logp(level))write(logu,format)trim(logl(level,__FILE__,__LINE__))," ",
#define log_root(level,format) if(logp(level,0))write(logu,format)trim(logl(level,__FILE__,__LINE__))," ",

/*
Tandem Repeats Finder 
Copyright (C) 1999-2020 Gary Benson

This file is part of the Tandem Repeats Finder (TRF) program.

TRF is free software: you can redistribute it and/or modify 
it under the terms of the GNU Affero General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

TRF is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See 
the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public 
License along with TRF.  If not, see <https://www.gnu.org/licenses/>.
*/

/* Platform selection and platform-dependent limits. */

#ifndef TRF_CORE_PLATFORM_H
#define TRF_CORE_PLATFORM_H

/* One platform identifier is in effect: WINDOWSGUI, WINDOWSCONSOLE,
   UNIXCONSOLE, or UNIXGUI. The console platform of the compiler's target is
   selected here. */
#if ( defined( WINDOWSGUI ) + defined( WINDOWSCONSOLE ) + defined( UNIXGUI ) + \
      defined( UNIXCONSOLE ) ) > 1
#error Only one platform can be defined
#endif

#if __unix__
#define UNIXCONSOLE
#undef WINDOWSGUI
#undef WINDOWSCONSOLE
#elif _WIN32
#define WINDOWSCONSOLE
#undef UNIXCONSOLE
#undef UNIXGUI
#endif

#ifdef UNIXCONSOLE
#define _MAX_PATH 260
#elif defined( UNIXGUI )

/* it is possible to do a windows GTK+ version */
#if !defined( _WIN32 )
#define _MAX_PATH 260  /* max. length of full pathname */
#define _MAX_DIR 260   /* max. length of path component */
#define _MAX_FNAME 260 /* max. length of file name component */
#define _MAX_EXT 260   /* max. length of extension component */
#endif

#endif

#ifndef _MAX_PATH
#define _MAX_PATH 260
#endif

#define farcalloc calloc

#ifndef UNIXGUI
/* alredy defined in glib, used in unix graphical version */
#define TRUE 1
#define FALSE 0
#endif

#endif

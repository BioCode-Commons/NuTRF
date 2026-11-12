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

/* Messages to the alignment file, the error stream, and the progress indicator. */

#include "core/globals.h"
#include "core/report.h"

/* explicit headers for the windows GUI definitions */
#ifdef WINDOWSGUI
 #include <windows.h>
 #include <windef.h>
 #include <winuser.h>
#endif

void trf_message( char *format, ... ) {

    va_list argp;

    if ( format == NULL )
        return;

    va_start( argp, format );

    if ( !paramset.HTMLoff )
        vfprintf( Fptxt, format, argp );

    va_end( argp );
}

void PrintError(char* errortext)
{
	/* the code here depends on the environment defined above */
#if defined(WINDOWSGUI)
	MessageBox((HWND)paramset.guihandle,errortext,"TRF",
			MB_ICONERROR|MB_OK);
#elif defined(WINDOWSCONSOLE)
	fprintf(stderr,"Error: %s\n",errortext);
#elif defined(UNIXGUI)
	unix_gui_error("Error: %s\n",errortext);
#elif defined(UNIXCONSOLE)
	fprintf(stderr,"Error: %s\n",errortext);
#endif

	return;
}

void PrintProgress(char* progresstext)
{
	/* platform dependent code */
#if defined(WINDOWSGUI)
	/* send message to status bar */
	SendDlgItemMessage((HWND)paramset.guihandle,IDC_STATUSBAR,
			SB_SETTEXT, (WPARAM) 0,(LPARAM)progresstext);
	SendDlgItemMessage((HWND)paramset.guihandle,IDC_STATUSBAR,
			WM_PAINT, (WPARAM) 0,(LPARAM)0);
#elif defined(WINDOWSCONSOLE)
	fprintf(stderr,"\n%s",progresstext);
#elif defined(UNIXGUI)
	print_progress("\n%s",progresstext);
#elif defined(UNIXCONSOLE)
	fprintf(stderr,"\n%s",progresstext);
#endif

	return;
}

void SetProgressBar(void)
{
	static int ready=0;

#if defined(WINDOWSGUI)
	static HWND barhwnd;
	static HDC hdc;
	static RECT rect,rect2;
	static char buffer[100];
#elif defined(WINDOWSCONSOLE)

#elif defined(UNIXGUI)

#elif defined(UNIXCONSOLE)

#endif

	/* if percent is minus one then un-ready the progress indicator */
	if(paramset.percent<0)
	{
		if(!ready) return;
		ready = 0;
		/* hide the progress indicator */
#if defined(WINDOWSGUI)
		InvalidateRect(barhwnd,NULL,TRUE);
		ReleaseDC(barhwnd,hdc);
#elif defined(WINDOWSCONSOLE)
		putchar('\n');
#elif defined(UNIXGUI)
		set_progress_bar( 0.0 );
#elif defined(UNIXCONSOLE)
		putchar('\n');
#endif

		return;
	}


	/* if not ready make ready */
	if(!ready)
	{
		ready = 1;

		/* position the progress progress indicator */
#if defined(WINDOWSGUI)
		/* get the rectangle of the second part of status bar */
		barhwnd = GetDlgItem((HWND)paramset.guihandle,IDC_STATUSBAR);
		SendMessage(barhwnd, SB_GETRECT,(WPARAM) 1,
				(LPARAM)(LPRECT) &rect);
		/* adjust for borders */
		rect.left+=2;
		rect.right-=2;
		rect.top+=2;
		rect.bottom-=2;
		/* set up progress rectangle */
		rect2 = rect;
		rect2.right = rect2.left;

		/* obtain a dc and draw on second panel */
		hdc = GetDC(barhwnd);
		FillRect(hdc,&rect,(HBRUSH)(COLOR_WINDOW+1));
		FillRect(hdc,&rect2,(HBRUSH)(COLOR_ACTIVECAPTION+1));

#elif defined(WINDOWSCONSOLE)
		putchar('\n');
#elif defined(UNIXGUI)

#elif defined(UNIXCONSOLE)
		putchar('\n');
#endif
	}

	/* set progress indicator to percent value */
	if(paramset.percent>100) paramset.percent=100;

#if defined(WINDOWSGUI)
	FillRect(hdc,&rect,(HBRUSH)(COLOR_WINDOW+1));
	rect2.right = rect2.left+(rect.right-rect.left)*paramset.percent/100;
	FillRect(hdc,&rect2,(HBRUSH)(COLOR_ACTIVECAPTION+1));
#elif defined(WINDOWSCONSOLE)
	if(!(paramset.percent%5)) putchar('.');
#elif defined(UNIXGUI)
	set_progress_bar( paramset.percent / 100.0 );
#elif defined(UNIXCONSOLE)
	if(!(paramset.percent%5)) putchar('.');
#endif

	return;
}

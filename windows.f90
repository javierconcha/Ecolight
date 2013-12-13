!******************************************************************************
!                                                                             *
!         WINDOWS.F90 -   Windows functions, types, and definitions           *
!                                                                             *
!                       Lahey Computer Systems, Inc.                          *
!                                PO Box 6091                                  *
!                         Incline Village, NV 89450                           *
!                            voice: 702-831-2500                              *
!                              fax: 702-831-8123                              *
!                              bbs: 702-831-8023                              *
!                        e-mail: support@lahey.com                            *
!                                                                             *
!              Copyright(c)  1996 - 1997, Lahey Computer Systems, Inc.               *
!                                                                             *
!******************************************************************************

module windows_h

!***** System Information *************************************************

INTEGER, PARAMETER :: NULL = 0

integer, parameter :: MAX_NEG = -2147483647
integer, parameter :: GFSR_SYSTEMRESOURCES =             0
integer, parameter :: GFSR_GDIRESOURCES =             1
integer, parameter :: GFSR_USERRESOURCES =             2

integer, parameter :: WF_PMODE =             1
integer, parameter :: WF_CPU286 =             2
integer, parameter :: WF_CPU386 =             4
integer, parameter :: WF_CPU486 =             8
integer, parameter :: WF_STANDARD =            16
integer, parameter :: WF_WIN286 =            16
integer, parameter :: WF_ENHANCED =            32
integer, parameter :: WF_WIN386 =            32
integer, parameter :: WF_CPU086 =            64
integer, parameter :: WF_CPU186 =           128
integer, parameter :: WF_LARGEFRAME =           256
integer, parameter :: WF_SMALLFRAME =           512
integer, parameter :: WF_80x87 =          1024
integer, parameter :: WF_PAGING =          2048
integer, parameter :: WF_WLO =         32768

!***** Error handling *****************************************************

!***** LogParamError/LogError values

! Error modifier bits

integer, parameter :: ERR_WARNING =         32768
integer, parameter :: ERR_PARAM =         16384

integer, parameter :: ERR_SIZE_MASK =         12288
integer, parameter :: ERR_BYTE =          4096
integer, parameter :: ERR_WORD =          8192
integer, parameter :: ERR_DWORD =         12288

!***** LogParamError() values

! Generic parameter values
integer, parameter :: ERR_BAD_VALUE =         24577
integer, parameter :: ERR_BAD_FLAGS =         24578
integer, parameter :: ERR_BAD_INDEX =         24579
integer, parameter :: ERR_BAD_DVALUE =         28676
integer, parameter :: ERR_BAD_DFLAGS =         28677
integer, parameter :: ERR_BAD_DINDEX =         28678
integer, parameter :: ERR_BAD_PTR =         28679
integer, parameter :: ERR_BAD_FUNC_PTR =         28680
integer, parameter :: ERR_BAD_SELECTOR =         24585
integer, parameter :: ERR_BAD_STRING_PTR =         28682
integer, parameter :: ERR_BAD_HANDLE =         24587

! KERNEL parameter errors
integer, parameter :: ERR_BAD_HINSTANCE =         24608
integer, parameter :: ERR_BAD_HMODULE =         24609
integer, parameter :: ERR_BAD_GLOBAL_HANDLE =         24610
integer, parameter :: ERR_BAD_LOCAL_HANDLE =         24611
integer, parameter :: ERR_BAD_ATOM =         24612
integer, parameter :: ERR_BAD_HFILE =         24613

! USER parameter errors
integer, parameter :: ERR_BAD_HWND =         24640
integer, parameter :: ERR_BAD_HMENU =         24641
integer, parameter :: ERR_BAD_HCURSOR =         24642
integer, parameter :: ERR_BAD_HICON =         24643
integer, parameter :: ERR_BAD_HDWP =         24644
integer, parameter :: ERR_BAD_CID =         24645
integer, parameter :: ERR_BAD_HDRVR =         24646

! GDI parameter errors
integer, parameter :: ERR_BAD_COORDS =         28768
integer, parameter :: ERR_BAD_GDI_OBJECT =         24673
integer, parameter :: ERR_BAD_HDC =         24674
integer, parameter :: ERR_BAD_HPEN =         24675
integer, parameter :: ERR_BAD_HFONT =         24676
integer, parameter :: ERR_BAD_HBRUSH =         24677
integer, parameter :: ERR_BAD_HBITMAP =         24678
integer, parameter :: ERR_BAD_HRGN =         24679
integer, parameter :: ERR_BAD_HPALETTE =         24680
integer, parameter :: ERR_BAD_HMETAFILE =         24681


!*** LogError() values

! KERNEL errors
integer, parameter :: ERR_GALLOC =             1
integer, parameter :: ERR_GREALLOC =             2
integer, parameter :: ERR_GLOCK =             3
integer, parameter :: ERR_LALLOC =             4
integer, parameter :: ERR_LREALLOC =             5
integer, parameter :: ERR_LLOCK =             6
integer, parameter :: ERR_ALLOCRES =             7
integer, parameter :: ERR_LOCKRES =             8
integer, parameter :: ERR_LOADMODULE =             9

! USER errors
integer, parameter :: ERR_CREATEDLG =            64
integer, parameter :: ERR_CREATEDLG2 =            65
integer, parameter :: ERR_REGISTERCLASS =            66
integer, parameter :: ERR_DCBUSY =            67
integer, parameter :: ERR_CREATEWND =            68
integer, parameter :: ERR_STRUCEXTRA =            69
integer, parameter :: ERR_LOADSTR =            70
integer, parameter :: ERR_LOADMENU =            71
integer, parameter :: ERR_NESTEDBEGINPAINT =            72
integer, parameter :: ERR_BADINDEX =            73
integer, parameter :: ERR_CREATEMENU =            74

! GDI errors
integer, parameter :: ERR_CREATEDC =           128
integer, parameter :: ERR_CREATEMETA =           129
integer, parameter :: ERR_DELOBJSELECTED =           130
integer, parameter :: ERR_SELBITMAP =           131

! Debugging support (DEBUG SYSTEM ONLY)

! WINDEBUGINFO flags values
integer, parameter :: WDI_OPTIONS =             1
integer, parameter :: WDI_FILTER =             2
integer, parameter :: WDI_ALLOCBREAK =             4

! dwOptions values
integer, parameter :: DBO_CHECKHEAP =             1
integer, parameter :: DBO_BUFFERFILL =             4
integer, parameter :: DBO_DISABLEGPTRAPPING =            16
integer, parameter :: DBO_CHECKFREE =            32

integer, parameter :: DBO_SILENT =         32768

integer, parameter :: DBO_TRACEBREAK =          8192
integer, parameter :: DBO_WARNINGBREAK =          4096
integer, parameter :: DBO_NOERRORBREAK =          2048
integer, parameter :: DBO_NOFATALBREAK =          1024
integer, parameter :: DBO_INT3BREAK =           256

! DebugOutput flags values
integer, parameter :: DBF_TRACE =             0
integer, parameter :: DBF_WARNING =         16384
integer, parameter :: DBF_ERROR =         32768
integer, parameter :: DBF_FATAL =         49152

! dwFilter values
integer, parameter :: DBF_KERNEL =          4096
integer, parameter :: DBF_KRN_MEMMAN =             1
integer, parameter :: DBF_KRN_LOADMODULE =             2
integer, parameter :: DBF_KRN_SEGMENTLOAD =             4
integer, parameter :: DBF_USER =          2048
integer, parameter :: DBF_GDI =          1024
integer, parameter :: DBF_MMSYSTEM =            64
integer, parameter :: DBF_PENWIN =            32
integer, parameter :: DBF_APPLICATION =             8
integer, parameter :: DBF_DRIVER =            16


integer, parameter :: EW_RESTARTWINDOWS =            66
integer, parameter :: EW_REBOOTSYSTEM =            67

! SetErrorMode() constants
integer, parameter :: SEM_FAILCRITICALERRORS =             1
integer, parameter :: SEM_NOGPFAULTERRORBOX =             2
integer, parameter :: SEM_NOOPENFILEERRORBOX =         32768


!***** Module Management **************************************************

integer, parameter :: HINSTANCE_ERROR = 32

! Windows Exit Procedure flag values
integer, parameter :: WEP_SYSTEM_EXIT = 1
integer, parameter :: WEP_FREE_DLL = 0


!***** Global memory management *******************************************

! Global Memory Flags

integer, parameter :: GMEM_FIXED =             0
integer, parameter :: GMEM_MOVEABLE =             2
integer, parameter :: GMEM_NOCOMPACT =            16
integer, parameter :: GMEM_NODISCARD =            32
integer, parameter :: GMEM_ZEROINIT =            64
integer, parameter :: GMEM_MODIFY =           128
integer, parameter :: GMEM_DISCARDABLE =           256
integer, parameter :: GMEM_NOT_BANKED =          4096
integer, parameter :: GMEM_SHARE =          8192
integer, parameter :: GMEM_DDESHARE =          8192
integer, parameter :: GMEM_NOTIFY =         16384
integer, parameter :: GMEM_LOWER = GMEM_NOT_BANKED

integer, parameter :: GHND = ior(GMEM_MOVEABLE, GMEM_ZEROINIT)
integer, parameter :: GPTR = ior(GMEM_FIXED, GMEM_ZEROINIT)

! GlobalFlags return flags (in addition to GMEM_DISCARDABLE)
integer, parameter :: GMEM_DISCARDED =         16384
integer, parameter :: GMEM_LOCKCOUNT =           255

! Low system memory notification message
integer, parameter :: WM_COMPACTING =            65

!**** Local Memory Management

! Local Memory Flags
integer, parameter :: LMEM_FIXED =             0
integer, parameter :: LMEM_MOVEABLE =             2
integer, parameter :: LMEM_NOCOMPACT =            16
integer, parameter :: LMEM_NODISCARD =            32
integer, parameter :: LMEM_ZEROINIT =            64
integer, parameter :: LMEM_MODIFY =           128
integer, parameter :: LMEM_DISCARDABLE =          3840

integer, parameter :: LHND = ior(LMEM_MOVEABLE, LMEM_ZEROINIT)
integer, parameter :: LPTR = ior(LMEM_FIXED, LMEM_ZEROINIT)

integer, parameter :: NONZEROLHND = (LMEM_MOVEABLE)
integer, parameter :: NONZEROLPTR = (LMEM_FIXED)

! LocalFlags return flags (in addition to LMEM_DISCARDABLE)
integer, parameter :: LMEM_DISCARDED =         16384
integer, parameter :: LMEM_LOCKCOUNT =           255


!***** File I/O ***********************************************************

integer, parameter :: HFILE_ERROR = -1

! OpenFile() Flags
integer, parameter :: OF_READ =             0
integer, parameter :: OF_WRITE =             1
integer, parameter :: OF_READWRITE =             2
integer, parameter :: OF_SHARE_COMPAT =             0
integer, parameter :: OF_SHARE_EXCLUSIVE =            16
integer, parameter :: OF_SHARE_DENY_WRITE =            32
integer, parameter :: OF_SHARE_DENY_READ =            48
integer, parameter :: OF_SHARE_DENY_NONE =            64
integer, parameter :: OF_PARSE =           256
integer, parameter :: OF_DELETE =           512
integer, parameter :: OF_VERIFY =          1024      ! Used with OF_REOPEN
integer, parameter :: OF_SEARCH =          1024      ! Used without OF_REOPEN
integer, parameter :: OF_CANCEL =          2048
integer, parameter :: OF_CREATE =          4096
integer, parameter :: OF_PROMPT =          8192
integer, parameter :: OF_EXIST =         16384
integer, parameter :: OF_REOPEN =         32768

! _lopen() flags
integer, parameter :: READ = 0
integer, parameter :: WRITE = 1
integer, parameter :: READ_WRITE = 2

! _llseek origin values
integer, parameter :: SEEK_SET = 0
integer, parameter :: SEEK_CUR = 1
integer, parameter :: SEEK_END = 2

! GetTempFileName() Flags
integer, parameter :: TF_FORCEDRIVE =           128

! GetDriveType return values
integer, parameter :: DRIVE_REMOVABLE = 2
integer, parameter :: DRIVE_FIXED = 3
integer, parameter :: DRIVE_REMOTE = 4

!***** Network support ****************************************************
! Errors
integer, parameter :: WN_SUCCESS =             0
integer, parameter :: WN_NOT_SUPPORTED =             1
integer, parameter :: WN_NET_ERROR =             2
integer, parameter :: WN_MORE_DATA =             3
integer, parameter :: WN_BAD_POINTER =             4
integer, parameter :: WN_BAD_VALUE =             5
integer, parameter :: WN_BAD_PASSWORD =             6
integer, parameter :: WN_ACCESS_DENIED =             7
integer, parameter :: WN_FUNCTION_BUSY =             8
integer, parameter :: WN_WINDOWS_ERROR =             9
integer, parameter :: WN_BAD_USER =            10
integer, parameter :: WN_OUT_OF_MEMORY =            11
integer, parameter :: WN_CANCEL =            12
integer, parameter :: WN_CONTINUE =            13

! Connection errors
integer, parameter :: WN_NOT_CONNECTED =            48
integer, parameter :: WN_OPEN_FILES =            49
integer, parameter :: WN_BAD_NETNAME =            50
integer, parameter :: WN_BAD_LOCALNAME =            51
integer, parameter :: WN_ALREADY_CONNECTED =            52
integer, parameter :: WN_DEVICE_ERROR =            53
integer, parameter :: WN_CONNECTION_CLOSED =            54

!***** Resource Management ************************************************
! Predefined Resource Types
integer, parameter :: RT_CURSOR            = 1
integer, parameter :: RT_BITMAP            = 2
integer, parameter :: RT_ICON              = 3
integer, parameter :: RT_MENU              = 4
integer, parameter :: RT_DIALOG            = 5
integer, parameter :: RT_STRING            = 6
integer, parameter :: RT_FONTDIR           = 7
integer, parameter :: RT_FONT              = 8
integer, parameter :: RT_ACCELERATOR       = 9
integer, parameter :: RT_RCDATA            = 10

integer, parameter :: RT_GROUP_CURSOR      = 12
integer, parameter :: RT_GROUP_ICON        = 14

! Standard cursor resource IDs
integer, parameter :: IDC_ARROW            = 32512
integer, parameter :: IDC_IBEAM            = 32513
integer, parameter :: IDC_WAIT             = 32514
integer, parameter :: IDC_CROSS            = 32515
integer, parameter :: IDC_UPARROW          = 32516
integer, parameter :: IDC_SIZE             = 32640
integer, parameter :: IDC_ICON             = 32641
integer, parameter :: IDC_SIZENWSE         = 32642
integer, parameter :: IDC_SIZENESW         = 32643
integer, parameter :: IDC_SIZEWE           = 32644
integer, parameter :: IDC_SIZENS           = 32645

! Standard icon resource IDs
integer, parameter :: IDI_APPLICATION      = 32512
integer, parameter :: IDI_HAND             = 32513
integer, parameter :: IDI_QUESTION         = 32514
integer, parameter :: IDI_EXCLAMATION      = 32515
integer, parameter :: IDI_ASTERISK         = 32516


! OEM Resource Ordinal Numbers
integer, parameter :: OBM_CLOSE = 32754
integer, parameter :: OBM_UPARROW = 32753
integer, parameter :: OBM_DNARROW = 32752
integer, parameter :: OBM_RGARROW = 32751
integer, parameter :: OBM_LFARROW = 32750
integer, parameter :: OBM_REDUCE = 32749
integer, parameter :: OBM_ZOOM = 32748
integer, parameter :: OBM_RESTORE = 32747
integer, parameter :: OBM_REDUCED = 32746
integer, parameter :: OBM_ZOOMD = 32745
integer, parameter :: OBM_RESTORED = 32744
integer, parameter :: OBM_UPARROWD = 32743
integer, parameter :: OBM_DNARROWD = 32742
integer, parameter :: OBM_RGARROWD = 32741
integer, parameter :: OBM_LFARROWD = 32740
integer, parameter :: OBM_MNARROW = 32739
integer, parameter :: OBM_COMBO = 32738
integer, parameter :: OBM_UPARROWI = 32737
integer, parameter :: OBM_DNARROWI = 32736
integer, parameter :: OBM_RGARROWI = 32735
integer, parameter :: OBM_LFARROWI = 32734

integer, parameter :: OBM_OLD_CLOSE = 32767
integer, parameter :: OBM_SIZE = 32766
integer, parameter :: OBM_OLD_UPARROW = 32765
integer, parameter :: OBM_OLD_DNARROW = 32764
integer, parameter :: OBM_OLD_RGARROW = 32763
integer, parameter :: OBM_OLD_LFARROW = 32762
integer, parameter :: OBM_BTSIZE = 32761
integer, parameter :: OBM_CHECK = 32760
integer, parameter :: OBM_CHECKBOXES = 32759
integer, parameter :: OBM_BTNCORNERS = 32758
integer, parameter :: OBM_OLD_REDUCE = 32757
integer, parameter :: OBM_OLD_ZOOM = 32756
integer, parameter :: OBM_OLD_RESTORE = 32755

integer, parameter :: OCR_NORMAL = 32512
integer, parameter :: OCR_IBEAM = 32513
integer, parameter :: OCR_WAIT = 32514
integer, parameter :: OCR_CROSS = 32515
integer, parameter :: OCR_UP = 32516
integer, parameter :: OCR_SIZE = 32640
integer, parameter :: OCR_ICON = 32641
integer, parameter :: OCR_SIZENWSE = 32642
integer, parameter :: OCR_SIZENESW = 32643
integer, parameter :: OCR_SIZEWE = 32644
integer, parameter :: OCR_SIZENS = 32645
integer, parameter :: OCR_SIZEALL = 32646
integer, parameter :: OCR_ICOCUR = 32647

integer, parameter :: OIC_SAMPLE = 32512
integer, parameter :: OIC_HAND = 32513
integer, parameter :: OIC_QUES = 32514
integer, parameter :: OIC_BANG = 32515
integer, parameter :: OIC_NOTE = 32516

!***** WIN.INI Support ****************************************************

integer, parameter :: WM_WININICHANGE =            26


!***** DC Management ******************************************************

! Drawing bounds accumulation APIs
integer, parameter :: DCB_RESET =             1
integer, parameter :: DCB_ACCUMULATE =             2
integer, parameter :: DCB_DIRTY = DCB_ACCUMULATE
integer, parameter :: DCB_SET = ior(DCB_RESET, DCB_ACCUMULATE)
integer, parameter :: DCB_ENABLE =             4
integer, parameter :: DCB_DISABLE =             8

!***** Device Capabilities ************************************************

! Device Parameters for GetDeviceCaps()
integer, parameter :: DRIVERVERSION = 0
integer, parameter :: TECHNOLOGY = 2
integer, parameter :: HORZSIZE = 4
integer, parameter :: VERTSIZE = 6
integer, parameter :: HORZRES = 8
integer, parameter :: VERTRES = 10
integer, parameter :: BITSPIXEL = 12
integer, parameter :: PLANES = 14
integer, parameter :: NUMBRUSHES = 16
integer, parameter :: NUMPENS = 18
integer, parameter :: NUMMARKERS = 20
integer, parameter :: NUMFONTS = 22
integer, parameter :: NUMCOLORS = 24
integer, parameter :: PDEVICESIZE = 26
integer, parameter :: CURVECAPS = 28
integer, parameter :: LINECAPS = 30
integer, parameter :: POLYGONALCAPS = 32
integer, parameter :: TEXTCAPS = 34
integer, parameter :: CLIPCAPS = 36
integer, parameter :: RASTERCAPS = 38
integer, parameter :: ASPECTX = 40
integer, parameter :: ASPECTY = 42
integer, parameter :: ASPECTXY = 44

integer, parameter :: LOGPIXELSX = 88
integer, parameter :: LOGPIXELSY = 90

integer, parameter :: SIZEPALETTE = 104
integer, parameter :: NUMRESERVED = 106
integer, parameter :: COLORRES = 108


! GetDeviceCaps() return value masks

! TECHNOLOGY
integer, parameter :: DT_PLOTTER = 0
integer, parameter :: DT_RASDISPLAY = 1
integer, parameter :: DT_RASPRINTER = 2
integer, parameter :: DT_RASCAMERA = 3
integer, parameter :: DT_CHARSTREAM = 4
integer, parameter :: DT_METAFILE = 5
integer, parameter :: DT_DISPFILE = 6

! CURVECAPS
integer, parameter :: CC_NONE =             0
integer, parameter :: CC_CIRCLES =             1
integer, parameter :: CC_PIE =             2
integer, parameter :: CC_CHORD =             4
integer, parameter :: CC_ELLIPSES =             8
integer, parameter :: CC_WIDE =            16
integer, parameter :: CC_STYLED =            32
integer, parameter :: CC_WIDESTYLED =            64
integer, parameter :: CC_INTERIORS =           128
integer, parameter :: CC_ROUNDRECT =           256

! LINECAPS
integer, parameter :: LC_NONE =             0
integer, parameter :: LC_POLYLINE =             2
integer, parameter :: LC_MARKER =             4
integer, parameter :: LC_POLYMARKER =             8
integer, parameter :: LC_WIDE =            16
integer, parameter :: LC_STYLED =            32
integer, parameter :: LC_WIDESTYLED =            64
integer, parameter :: LC_INTERIORS =           128

! POLYGONALCAPS
integer, parameter :: PC_NONE =             0
integer, parameter :: PC_POLYGON =             1
integer, parameter :: PC_RECTANGLE =             2
integer, parameter :: PC_WINDPOLYGON =             4
integer, parameter :: PC_SCANLINE =             8
integer, parameter :: PC_WIDE =            16
integer, parameter :: PC_STYLED =            32
integer, parameter :: PC_WIDESTYLED =            64
integer, parameter :: PC_INTERIORS =           128

! TEXTCAPS
integer, parameter :: TC_OP_CHARACTER =             1
integer, parameter :: TC_OP_STROKE =             2
integer, parameter :: TC_CP_STROKE =             4
integer, parameter :: TC_CR_90 =             8
integer, parameter :: TC_CR_ANY =            16
integer, parameter :: TC_SF_X_YINDEP =            32
integer, parameter :: TC_SA_DOUBLE =            64
integer, parameter :: TC_SA_INTEGER =           128
integer, parameter :: TC_SA_CONTIN =           256
integer, parameter :: TC_EA_DOUBLE =           512
integer, parameter :: TC_IA_ABLE =          1024
integer, parameter :: TC_UA_ABLE =          2048
integer, parameter :: TC_SO_ABLE =          4096
integer, parameter :: TC_RA_ABLE =          8192
integer, parameter :: TC_VA_ABLE =         16384
integer, parameter :: TC_RESERVED =         32768

! CLIPCAPS
integer, parameter :: CP_NONE =             0
integer, parameter :: CP_RECTANGLE =             1
integer, parameter :: CP_REGION =             2

! RASTERCAPS
integer, parameter :: RC_NONE = 0
integer, parameter :: RC_BITBLT =             1
integer, parameter :: RC_BANDING =             2
integer, parameter :: RC_SCALING =             4
integer, parameter :: RC_BITMAP64 =             8
integer, parameter :: RC_GDI20_OUTPUT =            16
integer, parameter :: RC_GDI20_STATE =            32
integer, parameter :: RC_SAVEBITMAP =            64
integer, parameter :: RC_DI_BITMAP =           128
integer, parameter :: RC_PALETTE =           256
integer, parameter :: RC_DIBTODEV =           512
integer, parameter :: RC_BIGFONT =          1024
integer, parameter :: RC_STRETCHBLT =          2048
integer, parameter :: RC_FLOODFILL =          4096
integer, parameter :: RC_STRETCHDIB =          8192
integer, parameter :: RC_OP_DX_OUTPUT =         16384
integer, parameter :: RC_DEVBITS =         32768


!***** Coordinate transformation support **********************************

! Map modes
integer, parameter :: MM_TEXT = 1
integer, parameter :: MM_LOMETRIC = 2
integer, parameter :: MM_HIMETRIC = 3
integer, parameter :: MM_LOENGLISH = 4
integer, parameter :: MM_HIENGLISH = 5
integer, parameter :: MM_TWIPS = 6
integer, parameter :: MM_ISOTROPIC = 7
integer, parameter :: MM_ANISOTROPIC = 8

! Coordinate Modes
integer, parameter :: ABSOLUTE = 1
integer, parameter :: RELATIVE = 2

!***** Color support ******************************************************

integer, parameter :: COLOR_SCROLLBAR = 0
integer, parameter :: COLOR_BACKGROUND = 1
integer, parameter :: COLOR_ACTIVECAPTION = 2
integer, parameter :: COLOR_INACTIVECAPTION = 3
integer, parameter :: COLOR_MENU = 4
integer, parameter :: COLOR_WINDOW = 5
integer, parameter :: COLOR_WINDOWFRAME = 6
integer, parameter :: COLOR_MENUTEXT = 7
integer, parameter :: COLOR_WINDOWTEXT = 8
integer, parameter :: COLOR_CAPTIONTEXT = 9
integer, parameter :: COLOR_ACTIVEBORDER = 10
integer, parameter :: COLOR_INACTIVEBORDER = 11
integer, parameter :: COLOR_APPWORKSPACE = 12
integer, parameter :: COLOR_HIGHLIGHT = 13
integer, parameter :: COLOR_HIGHLIGHTTEXT = 14
integer, parameter :: COLOR_BTNFACE = 15
integer, parameter :: COLOR_BTNSHADOW = 16
integer, parameter :: COLOR_GRAYTEXT = 17
integer, parameter :: COLOR_BTNTEXT = 18
integer, parameter :: COLOR_INACTIVECAPTIONTEXT = 19
integer, parameter :: COLOR_BTNHIGHLIGHT = 20

integer, parameter :: WM_SYSCOLORCHANGE =            21

!***** GDI Object Support *************************************************

! Object types for EnumObjects()
integer, parameter :: OBJ_PEN = 1
integer, parameter :: OBJ_BRUSH = 2

!***** Pen support ********************************************************

! Pen Styles
integer, parameter :: PS_SOLID = 0
integer, parameter :: PS_DASH = 1
integer, parameter :: PS_DOT = 2
integer, parameter :: PS_DASHDOT = 3
integer, parameter :: PS_DASHDOTDOT = 4
integer, parameter :: PS_NULL = 5
integer, parameter :: PS_INSIDEFRAME = 6

! Stock pens for use with GetStockObject();
integer, parameter :: WHITE_PEN = 6
integer, parameter :: BLACK_PEN = 7
integer, parameter :: NULL_PEN = 8

!***** Brush support ******************************************************

! Brush Styles
integer, parameter :: BS_SOLID = 0
integer, parameter :: BS_NULL = 1
integer, parameter :: BS_HOLLOW =          BS_NULL
integer, parameter :: BS_HATCHED = 2
integer, parameter :: BS_PATTERN = 3
integer, parameter :: BS_INDEXED = 4
integer, parameter :: BS_DIBPATTERN = 5

! Hatch Styles
integer, parameter :: HS_HORIZONTAL = 0
integer, parameter :: HS_VERTICAL = 1
integer, parameter :: HS_FDIAGONAL = 2
integer, parameter :: HS_BDIAGONAL = 3
integer, parameter :: HS_CROSS = 4
integer, parameter :: HS_DIAGCROSS = 5

! Stock brushes for use with GetStockObject()
integer, parameter :: WHITE_BRUSH = 0
integer, parameter :: LTGRAY_BRUSH = 1
integer, parameter :: GRAY_BRUSH = 2
integer, parameter :: DKGRAY_BRUSH = 3
integer, parameter :: BLACK_BRUSH = 4
integer, parameter :: NULL_BRUSH = 5
integer, parameter :: HOLLOW_BRUSH =       NULL_BRUSH

!***** Region support *****************************************************

! Region type flags
integer, parameter :: ERROR = 0
integer, parameter :: NULLREGION = 1
integer, parameter :: SIMPLEREGION = 2
integer, parameter :: COMPLEXREGION = 3

! CombineRgn() command values
integer, parameter :: RGN_AND = 1
integer, parameter :: RGN_OR = 2
integer, parameter :: RGN_XOR = 3
integer, parameter :: RGN_DIFF = 4
integer, parameter :: RGN_COPY = 5

!***** Color palette Support ***********************************************


! Palette entry flags
integer, parameter :: PC_RESERVED =             1    ! palette index used for animation
integer, parameter :: PC_EXPLICIT =             2    ! palette index is explicit to device
integer, parameter :: PC_NOCOLLAPSE =             4    ! do not match color to system palette

! Logical Palette

! Get/SetSystemPaletteUse() values
integer, parameter :: SYSPAL_STATIC = 1
integer, parameter :: SYSPAL_NOSTATIC = 2

! Palette window messages
integer, parameter :: WM_QUERYNEWPALETTE =           783
integer, parameter :: WM_PALETTEISCHANGING =           784
integer, parameter :: WM_PALETTECHANGED =           785

!***** General drawing support *******************************************

! PolyFill Modes
integer, parameter :: ALTERNATE = 1
integer, parameter :: WINDING = 2

! ExtFloodFill style flags
integer, parameter ::  FLOODFILLBORDER = 0
integer, parameter ::  FLOODFILLSURFACE = 1

!***** Text support *******************************************************

integer, parameter :: ETO_GRAYED =             1
integer, parameter :: ETO_OPAQUE =             2
integer, parameter :: ETO_CLIPPED =             4

! DrawText() Format Flags
integer, parameter :: DT_TOP =             0
integer, parameter :: DT_LEFT =             0
integer, parameter :: DT_CENTER =             1
integer, parameter :: DT_RIGHT =             2
integer, parameter :: DT_VCENTER =             4
integer, parameter :: DT_BOTTOM =             8
integer, parameter :: DT_WORDBREAK =            16
integer, parameter :: DT_SINGLELINE =            32
integer, parameter :: DT_EXPANDTABS =            64
integer, parameter :: DT_TABSTOP =           128
integer, parameter :: DT_NOCLIP =           256
integer, parameter :: DT_EXTERNALLEADING =           512
integer, parameter :: DT_CALCRECT =          1024
integer, parameter :: DT_NOPREFIX =          2048
integer, parameter :: DT_INTERNAL =          4096


! Background Modes
integer, parameter :: TRANSPARENT = 1
integer, parameter :: OPAQUE = 2

! Text Alignment Options
integer, parameter :: TA_NOUPDATECP =             0
integer, parameter :: TA_UPDATECP =             1
integer, parameter :: TA_LEFT =             0
integer, parameter :: TA_RIGHT =             2
integer, parameter :: TA_CENTER =             6
integer, parameter :: TA_TOP =             0
integer, parameter :: TA_BOTTOM =             8
integer, parameter :: TA_BASELINE =            24

!***** Font support *******************************************************

! Logical Font
integer, parameter :: LF_FACESIZE = 32

! weight values
integer, parameter :: FW_DONTCARE = 0
integer, parameter :: FW_THIN = 100
integer, parameter :: FW_EXTRALIGHT = 200
integer, parameter :: FW_LIGHT = 300
integer, parameter :: FW_NORMAL = 400
integer, parameter :: FW_MEDIUM = 500
integer, parameter :: FW_SEMIBOLD = 600
integer, parameter :: FW_BOLD = 700
integer, parameter :: FW_EXTRABOLD = 800
integer, parameter :: FW_HEAVY = 900

integer, parameter :: FW_ULTRALIGHT =      FW_EXTRALIGHT
integer, parameter :: FW_REGULAR =         FW_NORMAL
integer, parameter :: FW_DEMIBOLD =        FW_SEMIBOLD
integer, parameter :: FW_ULTRABOLD =       FW_EXTRABOLD
integer, parameter :: FW_BLACK =           FW_HEAVY

! CharSet values
integer, parameter :: ANSI_CHARSET = 0
integer, parameter :: DEFAULT_CHARSET = 1
integer, parameter :: SYMBOL_CHARSET = 2
integer, parameter :: SHIFTJIS_CHARSET = 128
integer, parameter :: HANGEUL_CHARSET = 129
integer, parameter :: CHINESEBIG5_CHARSET = 136
integer, parameter :: OEM_CHARSET = 255

! OutPrecision values
integer, parameter :: OUT_DEFAULT_PRECIS = 0
integer, parameter :: OUT_STRING_PRECIS = 1
integer, parameter :: OUT_CHARACTER_PRECIS = 2
integer, parameter :: OUT_STROKE_PRECIS = 3
integer, parameter :: OUT_TT_PRECIS = 4
integer, parameter :: OUT_DEVICE_PRECIS = 5
integer, parameter :: OUT_RASTER_PRECIS = 6
integer, parameter :: OUT_TT_ONLY_PRECIS = 7

! ClipPrecision values
integer, parameter :: CLIP_DEFAULT_PRECIS =             0
integer, parameter :: CLIP_CHARACTER_PRECIS =             1
integer, parameter :: CLIP_STROKE_PRECIS =             2
integer, parameter :: CLIP_MASK =            15
integer, parameter :: CLIP_LH_ANGLES =            16
integer, parameter :: CLIP_TT_ALWAYS =            32
integer, parameter :: CLIP_EMBEDDED =           128

! Quality values
integer, parameter :: DEFAULT_QUALITY = 0
integer, parameter :: DRAFT_QUALITY = 1
integer, parameter :: PROOF_QUALITY = 2

! PitchAndFamily pitch values (low 4 bits)
integer, parameter :: DEFAULT_PITCH =             0
integer, parameter :: FIXED_PITCH =             1
integer, parameter :: VARIABLE_PITCH =             2

! PitchAndFamily family values (high 4 bits)
integer, parameter :: FF_DONTCARE =             0
integer, parameter :: FF_ROMAN =            16
integer, parameter :: FF_SWISS =            32
integer, parameter :: FF_MODERN =            48
integer, parameter :: FF_SCRIPT =            64
integer, parameter :: FF_DECORATIVE =            80

! Stock fonts for use with GetStockObject()
integer, parameter :: OEM_FIXED_FONT = 10
integer, parameter :: ANSI_FIXED_FONT = 11
integer, parameter :: ANSI_VAR_FONT = 12
integer, parameter :: SYSTEM_FONT = 13
integer, parameter :: DEVICE_DEFAULT_FONT = 14
integer, parameter :: DEFAULT_PALETTE = 15
integer, parameter :: SYSTEM_FIXED_FONT = 16

integer, parameter :: ASPECT_FILTERING =             1

integer, parameter :: WM_FONTCHANGE =            29

! tmPitchAndFamily values
integer, parameter :: TMPF_FIXED_PITCH =             1
integer, parameter :: TMPF_VECTOR =             2
integer, parameter :: TMPF_DEVICE =             8
integer, parameter :: TMPF_TRUETYPE =             4

! ntmFlags field flags
integer, parameter :: NTM_REGULAR =            64
integer, parameter :: NTM_BOLD =            32
integer, parameter :: NTM_ITALIC =             1

integer, parameter :: LF_FULLFACESIZE = 64

! EnumFonts font type values
integer, parameter :: RASTER_FONTTYPE =             1
integer, parameter :: DEVICE_FONTTYPE =             2
integer, parameter :: TRUETYPE_FONTTYPE =             4

! GetGlyphOutline constants
integer, parameter :: GGO_METRICS = 0
integer, parameter :: GGO_BITMAP = 1
integer, parameter :: GGO_NATIVE = 2

integer, parameter :: TT_POLYGON_TYPE = 24

integer, parameter :: TT_PRIM_LINE = 1
integer, parameter :: TT_PRIM_QSPLINE = 2

! bits defined in wFlags of RASTERIZER_STATUS
integer, parameter :: TT_AVAILABLE =             1
integer, parameter :: TT_ENABLED =             2

! constants for the biCompression field
integer, parameter :: BI_RGB = 0
integer, parameter :: BI_RLE8 = 1
integer, parameter :: BI_RLE4 = 2

! DIB color table identifiers
integer, parameter :: DIB_RGB_COLORS = 0
integer, parameter :: DIB_PAL_COLORS = 1

! constants for CreateDIBitmap
integer, parameter :: CBM_INIT =             4

! Binary raster ops
integer, parameter :: R2_BLACK = 1
integer, parameter :: R2_NOTMERGEPEN = 2
integer, parameter :: R2_MASKNOTPEN = 3
integer, parameter :: R2_NOTCOPYPEN = 4
integer, parameter :: R2_MASKPENNOT = 5
integer, parameter :: R2_NOT = 6
integer, parameter :: R2_XORPEN = 7
integer, parameter :: R2_NOTMASKPEN = 8
integer, parameter :: R2_MASKPEN = 9
integer, parameter :: R2_NOTXORPEN = 10
integer, parameter :: R2_NOP = 11
integer, parameter :: R2_MERGENOTPEN = 12
integer, parameter :: R2_COPYPEN = 13
integer, parameter :: R2_MERGEPENNOT = 14
integer, parameter :: R2_MERGEPEN = 15
integer, parameter :: R2_WHITE = 16

! Ternary raster operations
integer, parameter :: SRCCOPY =      13369376
integer, parameter :: SRCPAINT =      15597702
integer, parameter :: SRCAND =       8913094
integer, parameter :: SRCINVERT =       6684742
integer, parameter :: SRCERASE =       4457256
integer, parameter :: NOTSRCCOPY =       3342344
integer, parameter :: NOTSRCERASE =       1114278
integer, parameter :: MERGECOPY =      12583114
integer, parameter :: MERGEPAINT =      12255782
integer, parameter :: PATCOPY =      15728673
integer, parameter :: PATPAINT =      16452105
integer, parameter :: PATINVERT =       5898313
integer, parameter :: DSTINVERT =       5570569
integer, parameter :: BLACKNESS =            66
integer, parameter :: WHITENESS =      16711778

! StretchBlt() Modes
integer, parameter :: BLACKONWHITE = 1
integer, parameter :: WHITEONBLACK = 2
integer, parameter :: COLORONCOLOR = 3

! new StretchBlt() Modes (simpler names)
integer, parameter :: STRETCH_ANDSCANS = 1
integer, parameter :: STRETCH_ORSCANS = 2
integer, parameter :: STRETCH_DELETESCANS = 3

! Metafile Functions
integer, parameter :: META_SETBKCOLOR =           513
integer, parameter :: META_SETBKMODE =           258
integer, parameter :: META_SETMAPMODE =           259
integer, parameter :: META_SETROP2 =           260
integer, parameter :: META_SETRELABS =           261
integer, parameter :: META_SETPOLYFILLMODE =           262
integer, parameter :: META_SETSTRETCHBLTMODE =           263
integer, parameter :: META_SETTEXTCHAREXTRA =           264
integer, parameter :: META_SETTEXTCOLOR =           521
integer, parameter :: META_SETTEXTJUSTIFICATION =           522
integer, parameter :: META_SETWINDOWORG =           523
integer, parameter :: META_SETWINDOWEXT =           524
integer, parameter :: META_SETVIEWPORTORG =           525
integer, parameter :: META_SETVIEWPORTEXT =           526
integer, parameter :: META_OFFSETWINDOWORG =           527
integer, parameter :: META_SCALEWINDOWEXT =          1040
integer, parameter :: META_OFFSETVIEWPORTORG =           529
integer, parameter :: META_SCALEVIEWPORTEXT =          1042
integer, parameter :: META_LINETO =           531
integer, parameter :: META_MOVETO =           532
integer, parameter :: META_EXCLUDECLIPRECT =          1045
integer, parameter :: META_INTERSECTCLIPRECT =          1046
integer, parameter :: META_ARC =          2071
integer, parameter :: META_ELLIPSE =          1048
integer, parameter :: META_FLOODFILL =          1049
integer, parameter :: META_PIE =          2074
integer, parameter :: META_RECTANGLE =          1051
integer, parameter :: META_ROUNDRECT =          1564
integer, parameter :: META_PATBLT =          1565
integer, parameter :: META_SAVEDC =            30
integer, parameter :: META_SETPIXEL =          1055
integer, parameter :: META_OFFSETCLIPRGN =           544
integer, parameter :: META_TEXTOUT =          1313
integer, parameter :: META_BITBLT =          2338
integer, parameter :: META_STRETCHBLT =          2851
integer, parameter :: META_POLYGON =           804
integer, parameter :: META_POLYLINE =           805
integer, parameter :: META_ESCAPE =          1574
integer, parameter :: META_RESTOREDC =           295
integer, parameter :: META_FILLREGION =           552
integer, parameter :: META_FRAMEREGION =          1065
integer, parameter :: META_INVERTREGION =           298
integer, parameter :: META_PAINTREGION =           299
integer, parameter :: META_SELECTCLIPREGION =           300
integer, parameter :: META_SELECTOBJECT =           301
integer, parameter :: META_SETTEXTALIGN =           302
integer, parameter :: META_DRAWTEXT =          1583

integer, parameter :: META_CHORD =          2096
integer, parameter :: META_SETMAPPERFLAGS =           561
integer, parameter :: META_EXTTEXTOUT =          2610
integer, parameter :: META_SETDIBTODEV =          3379
integer, parameter :: META_SELECTPALETTE =           564
integer, parameter :: META_REALIZEPALETTE =            53
integer, parameter :: META_ANIMATEPALETTE =          1078
integer, parameter :: META_SETPALENTRIES =            55
integer, parameter :: META_POLYPOLYGON =          1336
integer, parameter :: META_RESIZEPALETTE =           313

integer, parameter :: META_DIBBITBLT =          2368
integer, parameter :: META_DIBSTRETCHBLT =          2881
integer, parameter :: META_DIBCREATEPATTERNBRUSH =           322
integer, parameter :: META_STRETCHDIB =          3907

integer, parameter :: META_EXTFLOODFILL =          1352

integer, parameter :: META_RESETDC =           332
integer, parameter :: META_STARTDOC =           333
integer, parameter :: META_STARTPAGE =            79
integer, parameter :: META_ENDPAGE =            80
integer, parameter :: META_ABORTDOC =            82
integer, parameter :: META_ENDDOC =            94

integer, parameter :: META_DELETEOBJECT =           496

integer, parameter :: META_CREATEPALETTE =           247
integer, parameter :: META_CREATEBRUSH =           248
integer, parameter :: META_CREATEPATTERNBRUSH =           505
integer, parameter :: META_CREATEPENINDIRECT =           762
integer, parameter :: META_CREATEFONTINDIRECT =           763
integer, parameter :: META_CREATEBRUSHINDIRECT =           764
integer, parameter :: META_CREATEBITMAPINDIRECT =           765
integer, parameter :: META_CREATEBITMAP =          1790
integer, parameter :: META_CREATEREGION =          1791


!***** Printing support ***************************************************

! Spooler Error Codes
integer, parameter :: SP_NOTREPORTED =         16384
integer, parameter :: SP_ERROR = (-1)
integer, parameter :: SP_APPABORT = (-2)
integer, parameter :: SP_USERABORT = (-3)
integer, parameter :: SP_OUTOFDISK = (-4)
integer, parameter :: SP_OUTOFMEMORY = (-5)

integer, parameter :: PR_JOBSTATUS =             0

! Spooler status notification message
integer, parameter :: WM_SPOOLERSTATUS =            42

!****** GDI Escape support ************************************************

! GDI Escapes
integer, parameter :: NEWFRAME = 1
integer, parameter :: ABORTDOC = 2
integer, parameter :: NEXTBAND = 3
integer, parameter :: SETCOLORTABLE = 4
integer, parameter :: GETCOLORTABLE = 5
integer, parameter :: FLUSHOUTPUT = 6
integer, parameter :: DRAFTMODE = 7
integer, parameter :: QUERYESCSUPPORT = 8
integer, parameter :: SETABORTPROC = 9
integer, parameter :: STARTDOC = 10
integer, parameter :: ENDDOC = 11
integer, parameter :: GETPHYSPAGESIZE = 12
integer, parameter :: GETPRINTINGOFFSET = 13
integer, parameter :: GETSCALINGFACTOR = 14
integer, parameter :: MFCOMMENT = 15
integer, parameter :: GETPENWIDTH = 16
integer, parameter :: SETCOPYCOUNT = 17
integer, parameter :: SELECTPAPERSOURCE = 18
integer, parameter :: DEVICEDATA = 19
integer, parameter :: PASSTHROUGH = 19
integer, parameter :: GETTECHNOLGY = 20
integer, parameter :: GETTECHNOLOGY = 20
integer, parameter :: SETLINECAP = 21
integer, parameter :: SETLINEJOIN = 22
integer, parameter :: SETMITERLIMIT = 23
integer, parameter :: BANDINFO = 24
integer, parameter :: DRAWPATTERNRECT = 25
integer, parameter :: GETVECTORPENSIZE = 26
integer, parameter :: GETVECTORBRUSHSIZE = 27
integer, parameter :: ENABLEDUPLEX = 28
integer, parameter :: GETSETPAPERBINS = 29
integer, parameter :: GETSETPRINTORIENT = 30
integer, parameter :: ENUMPAPERBINS = 31
integer, parameter :: SETDIBSCALING = 32
integer, parameter :: EPSPRINTING = 33
integer, parameter :: ENUMPAPERMETRICS = 34
integer, parameter :: GETSETPAPERMETRICS = 35
integer, parameter :: POSTSCRIPT_DATA = 37
integer, parameter :: POSTSCRIPT_IGNORE = 38
integer, parameter :: MOUSETRAILS = 39

integer, parameter :: GETEXTENDEDTEXTMETRICS = 256
integer, parameter :: GETEXTENTTABLE = 257
integer, parameter :: GETPAIRKERNTABLE = 258
integer, parameter :: GETTRACKKERNTABLE = 259
integer, parameter :: EXTTEXTOUT = 512
integer, parameter :: GETFACENAME = 513
integer, parameter :: ENABLERELATIVEWIDTHS = 768
integer, parameter :: ENABLEPAIRKERNING = 769
integer, parameter :: SETKERNTRACK = 770
integer, parameter :: SETALLJUSTVALUES = 771
integer, parameter :: SETCHARSET = 772

integer, parameter :: STRETCHBLT = 2048

integer, parameter :: GETSETSCREENPARAMS = 3072

integer, parameter :: BEGIN_PATH = 4096
integer, parameter :: CLIP_TO_PATH = 4097
integer, parameter :: END_PATH = 4098
integer, parameter :: EXT_DEVICE_CAPS = 4099
integer, parameter :: RESTORE_CTM = 4100
integer, parameter :: SAVE_CTM = 4101
integer, parameter :: SET_ARC_DIRECTION = 4102
integer, parameter :: SET_BACKGROUND_COLOR = 4103
integer, parameter :: SET_POLY_MODE = 4104
integer, parameter :: SET_SCREEN_ANGLE = 4105
integer, parameter :: SET_SPREAD = 4106
integer, parameter :: TRANSFORM_CTM = 4107
integer, parameter :: SET_CLIP_BOX = 4108
integer, parameter :: SET_BOUNDS = 4109

!***** USER typedefs, structures, and functions ****************************

! GetSystemMetrics() codes
integer, parameter :: SM_CXSCREEN = 0
integer, parameter :: SM_CYSCREEN = 1
integer, parameter :: SM_CXVSCROLL = 2
integer, parameter :: SM_CYHSCROLL = 3
integer, parameter :: SM_CYCAPTION = 4
integer, parameter :: SM_CXBORDER = 5
integer, parameter :: SM_CYBORDER = 6
integer, parameter :: SM_CXDLGFRAME = 7
integer, parameter :: SM_CYDLGFRAME = 8
integer, parameter :: SM_CYVTHUMB = 9
integer, parameter :: SM_CXHTHUMB = 10
integer, parameter :: SM_CXICON = 11
integer, parameter :: SM_CYICON = 12
integer, parameter :: SM_CXCURSOR = 13
integer, parameter :: SM_CYCURSOR = 14
integer, parameter :: SM_CYMENU = 15
integer, parameter :: SM_CXFULLSCREEN = 16
integer, parameter :: SM_CYFULLSCREEN = 17
integer, parameter :: SM_CYKANJIWINDOW = 18
integer, parameter :: SM_MOUSEPRESENT = 19
integer, parameter :: SM_CYVSCROLL = 20
integer, parameter :: SM_CXHSCROLL = 21
integer, parameter :: SM_DEBUG = 22
integer, parameter :: SM_SWAPBUTTON = 23
integer, parameter :: SM_RESERVED1 = 24
integer, parameter :: SM_RESERVED2 = 25
integer, parameter :: SM_RESERVED3 = 26
integer, parameter :: SM_RESERVED4 = 27
integer, parameter :: SM_CXMIN = 28
integer, parameter :: SM_CYMIN = 29
integer, parameter :: SM_CXSIZE = 30
integer, parameter :: SM_CYSIZE = 31
integer, parameter :: SM_CXFRAME = 32
integer, parameter :: SM_CYFRAME = 33
integer, parameter :: SM_CXMINTRACK = 34
integer, parameter :: SM_CYMINTRACK = 35

integer, parameter :: SM_CXDOUBLECLK = 36
integer, parameter :: SM_CYDOUBLECLK = 37
integer, parameter :: SM_CXICONSPACING = 38
integer, parameter :: SM_CYICONSPACING = 39
integer, parameter :: SM_MENUDROPALIGNMENT = 40
integer, parameter :: SM_PENWINDOWS = 41
integer, parameter :: SM_DBCSENABLED = 42

integer, parameter :: SM_CMETRICS = 43

integer, parameter :: WM_DEVMODECHANGE =            27
integer, parameter :: WM_TIMECHANGE =            30

!***** System Parameters support *******************************************

integer, parameter :: SPI_GETBEEP = 1
integer, parameter :: SPI_SETBEEP = 2
integer, parameter :: SPI_GETMOUSE = 3
integer, parameter :: SPI_SETMOUSE = 4
integer, parameter :: SPI_GETBORDER = 5
integer, parameter :: SPI_SETBORDER = 6
integer, parameter :: SPI_GETKEYBOARDSPEED = 10
integer, parameter :: SPI_SETKEYBOARDSPEED = 11
integer, parameter :: SPI_LANGDRIVER = 12
integer, parameter :: SPI_ICONHORIZONTALSPACING = 13
integer, parameter :: SPI_GETSCREENSAVETIMEOUT = 14
integer, parameter :: SPI_SETSCREENSAVETIMEOUT = 15
integer, parameter :: SPI_GETSCREENSAVEACTIVE = 16
integer, parameter :: SPI_SETSCREENSAVEACTIVE = 17
integer, parameter :: SPI_GETGRIDGRANULARITY = 18
integer, parameter :: SPI_SETGRIDGRANULARITY = 19
integer, parameter :: SPI_SETDESKWALLPAPER = 20
integer, parameter :: SPI_SETDESKPATTERN = 21
integer, parameter :: SPI_GETKEYBOARDDELAY = 22
integer, parameter :: SPI_SETKEYBOARDDELAY = 23
integer, parameter :: SPI_ICONVERTICALSPACING = 24
integer, parameter :: SPI_GETICONTITLEWRAP = 25
integer, parameter :: SPI_SETICONTITLEWRAP = 26
integer, parameter :: SPI_GETMENUDROPALIGNMENT = 27
integer, parameter :: SPI_SETMENUDROPALIGNMENT = 28
integer, parameter :: SPI_SETDOUBLECLKWIDTH = 29
integer, parameter :: SPI_SETDOUBLECLKHEIGHT = 30
integer, parameter :: SPI_GETICONTITLELOGFONT = 31
integer, parameter :: SPI_SETDOUBLECLICKTIME = 32
integer, parameter :: SPI_SETMOUSEBUTTONSWAP = 33
integer, parameter :: SPI_SETICONTITLELOGFONT = 34
integer, parameter :: SPI_GETFASTTASKSWITCH = 35
integer, parameter :: SPI_SETFASTTASKSWITCH = 36

! SystemParametersInfo flags
integer, parameter :: SPIF_UPDATEINIFILE =             1
integer, parameter :: SPIF_SENDWININICHANGE =             2

!***** Window message support **********************************************

integer, parameter :: WM_NULL =             0

! NOTE: All messages below =          1024 are RESERVED by Windows
integer, parameter :: WM_USER =          1024

! PeekMessage() options
integer, parameter :: PM_NOREMOVE =             0
integer, parameter :: PM_REMOVE =             1
integer, parameter :: PM_NOYIELD =             2

! GetQueueStatus flags
integer, parameter :: QS_KEY =             1
integer, parameter :: QS_MOUSEMOVE =             2
integer, parameter :: QS_MOUSEBUTTON =             4
integer, parameter :: QS_MOUSE =       ior(QS_MOUSEMOVE, QS_MOUSEBUTTON)
integer, parameter :: QS_POSTMESSAGE =             8
integer, parameter :: QS_TIMER =            16
integer, parameter :: QS_PAINT =            32
integer, parameter :: QS_SENDMESSAGE =            64

integer, parameter :: QS_ALLINPUT =           127

! Special HWND value for use with PostMessage() and SendMessage()
integer, parameter :: HWND_BROADCAST =         65535

integer, parameter :: WH_GETMESSAGE = 3

integer, parameter :: WH_CALLWNDPROC = 4

integer, parameter :: WH_MSGFILTER = (-1)
integer, parameter :: WH_SYSMSGFILTER = 6

! CallMsgFilter() and WH_SYS/MSGFILTER context codes
integer, parameter :: MSGF_DIALOGBOX = 0
integer, parameter :: MSGF_MENU = 2
integer, parameter :: MSGF_MOVE = 3
integer, parameter :: MSGF_SIZE = 4
integer, parameter :: MSGF_SCROLLBAR = 5
integer, parameter :: MSGF_NEXTWINDOW = 6
integer, parameter :: MSGF_MAINLOOP = 8
integer, parameter :: MSGF_USER = 4096

! Standard window messages
! PenWindows specific messages
integer, parameter :: WM_PENWINFIRST =           896
integer, parameter :: WM_PENWINLAST =           911

! Coalescing messages
integer, parameter :: WM_COALESCE_FIRST =           912
integer, parameter :: WM_COALESCE_LAST =           927

!***** Power management ***************************************************
integer, parameter :: WM_POWER =            72

! wParam for WM_POWER window message and DRV_POWER driver notification
integer, parameter :: PWR_OK = 1
integer, parameter :: PWR_FAIL = (-1)
integer, parameter :: PWR_SUSPENDREQUEST = 1
integer, parameter :: PWR_SUSPENDRESUME = 2
integer, parameter :: PWR_CRITICALRESUME = 3

!***** Application termination ********************************************

integer, parameter :: WM_QUERYENDSESSION =            17
integer, parameter :: WM_ENDSESSION =            22

integer, parameter :: WM_QUIT =            18

integer, parameter :: WM_SYSTEMERROR =            23

!***** Window class management ********************************************

! Class styles
integer, parameter :: CS_VREDRAW =             1
integer, parameter :: CS_HREDRAW =             2

integer, parameter :: CS_OWNDC =            32
integer, parameter :: CS_CLASSDC =            64
integer, parameter :: CS_PARENTDC =           128

integer, parameter :: CS_SAVEBITS =          2048

integer, parameter :: CS_DBLCLKS =             8

integer, parameter :: CS_BYTEALIGNCLIENT =          4096
integer, parameter :: CS_BYTEALIGNWINDOW =          8192

integer, parameter :: CS_NOCLOSE =           512

integer, parameter :: CS_KEYCVTWINDOW =             4
integer, parameter :: CS_NOKEYCVT =           256

integer, parameter :: CS_GLOBALCLASS =         16384

! Class field offsets for GetClassLong() and GetClassWord()
integer, parameter :: GCL_MENUNAME = (-8)
integer, parameter :: GCW_HBRBACKGROUND = (-10)
integer, parameter :: GCW_HCURSOR = (-12)
integer, parameter :: GCW_HICON = (-14)
integer, parameter :: GCW_HMODULE = (-16)
integer, parameter :: GCW_CBWNDEXTRA = (-18)
integer, parameter :: GCW_CBCLSEXTRA = (-20)
integer, parameter :: GCL_WNDPROC = (-24)
integer, parameter :: GCW_STYLE = (-26)

integer, parameter :: GCW_ATOM = (-32)

!***** Window creation/destroy ********************************************

! Window Styles

! Basic window types
integer, parameter :: WS_OVERLAPPED =             0
integer, parameter :: WS_POPUP = MAX_NEG - 1  ! -2147483648
integer, parameter :: WS_CHILD =    1073741824

! Clipping styles
integer, parameter :: WS_CLIPSIBLINGS =      67108864
integer, parameter :: WS_CLIPCHILDREN =      33554432

! Generic window states
integer, parameter :: WS_VISIBLE =     268435456
integer, parameter :: WS_DISABLED =     134217728

! Main window states
integer, parameter :: WS_MINIMIZE =     536870912
integer, parameter :: WS_MAXIMIZE =      16777216

! Main window styles
integer, parameter :: WS_CAPTION =      12582912     ! WS_BORDER | WS_DLGFRAME
integer, parameter :: WS_BORDER =       8388608
integer, parameter :: WS_DLGFRAME =       4194304
integer, parameter :: WS_VSCROLL =       2097152
integer, parameter :: WS_HSCROLL =       1048576
integer, parameter :: WS_SYSMENU =        524288
integer, parameter :: WS_THICKFRAME =        262144
integer, parameter :: WS_MINIMIZEBOX =        131072
integer, parameter :: WS_MAXIMIZEBOX =         65536

! Control window styles
integer, parameter :: WS_GROUP =        131072
integer, parameter :: WS_TABSTOP =         65536

! Common Window Styles
integer, parameter :: WS_OVERLAPPEDWINDOW = ior(WS_OVERLAPPED, ior(WS_CAPTION, ior(WS_SYSMENU, ior(WS_THICKFRAME, &
                                            ior(WS_MINIMIZEBOX, WS_MAXIMIZEBOX)))))
integer, parameter :: WS_POPUPWINDOW = ior(WS_POPUP,  ior(WS_BORDER, WS_SYSMENU))
integer, parameter :: WS_CHILDWINDOW = (WS_CHILD)

! Extended Window Styles
integer, parameter :: WS_EX_DLGMODALFRAME =             1
integer, parameter :: WS_EX_NOPARENTNOTIFY =             4

integer, parameter :: WS_EX_TOPMOST =             8
integer, parameter :: WS_EX_ACCEPTFILES =            16
integer, parameter :: WS_EX_TRANSPARENT =            32

! Special value for CreateWindow, et al.
integer, parameter :: HWND_DESKTOP = 0

integer, parameter :: WM_CREATE =             1
integer, parameter :: WM_NCCREATE =           129

integer, parameter :: WM_DESTROY =             2
integer, parameter :: WM_NCDESTROY =           130

! Basic window attributes

integer, parameter :: SW_HIDE = 0
integer, parameter :: SW_SHOWNORMAL = 1
integer, parameter :: SW_NORMAL = 1
integer, parameter :: SW_SHOWMINIMIZED = 2
integer, parameter :: SW_SHOWMAXIMIZED = 3
integer, parameter :: SW_MAXIMIZE = 3
integer, parameter :: SW_SHOWNOACTIVATE = 4
integer, parameter :: SW_SHOW = 5
integer, parameter :: SW_MINIMIZE = 6
integer, parameter :: SW_SHOWMINNOACTIVE = 7
integer, parameter :: SW_SHOWNA = 8
integer, parameter :: SW_RESTORE = 9

integer, parameter :: WM_SHOWWINDOW =            24

! WM_SHOWWINDOW wParam codes
integer, parameter :: SW_PARENTCLOSING = 1
integer, parameter :: SW_OTHERMAXIMIZED = 2
integer, parameter :: SW_PARENTOPENING = 3
integer, parameter :: SW_OTHERRESTORED = 4

integer, parameter :: WM_SETREDRAW =            11

integer, parameter :: WM_ENABLE =            10

! Window text
integer, parameter :: WM_SETTEXT =            12
integer, parameter :: WM_GETTEXT =            13
integer, parameter :: WM_GETTEXTLENGTH =            14

! Window field offsets for GetWindowLong() and GetWindowWord()
integer, parameter :: GWL_WNDPROC = (-4)
integer, parameter :: GWW_HINSTANCE = (-6)
integer, parameter :: GWW_HWNDPARENT = (-8)
integer, parameter :: GWW_ID = (-12)
integer, parameter :: GWL_STYLE = (-16)
integer, parameter :: GWL_EXSTYLE = (-20)

!***** Window size, position, Z-order, and visibility *********************

integer, parameter :: CW_USEDEFAULT =         32768

integer, parameter :: WPF_SETMINPOSITION =             1
integer, parameter :: WPF_RESTORETOMAXIMIZED =             2

! SetWindowPos() and WINDOWPOS flags
integer, parameter :: SWP_NOSIZE =             1
integer, parameter :: SWP_NOMOVE =             2
integer, parameter :: SWP_NOZORDER =             4
integer, parameter :: SWP_NOREDRAW =             8
integer, parameter :: SWP_NOACTIVATE =            16
integer, parameter :: SWP_FRAMECHANGED =            32  ! The frame changed: send WM_NCCALCSIZE
integer, parameter :: SWP_SHOWWINDOW =            64
integer, parameter :: SWP_HIDEWINDOW =           128
integer, parameter :: SWP_NOCOPYBITS =           256
integer, parameter :: SWP_NOOWNERZORDER =           512  ! Don't do owner Z ordering

integer, parameter :: SWP_DRAWFRAME =      SWP_FRAMECHANGED
integer, parameter :: SWP_NOREPOSITION =    SWP_NOOWNERZORDER

integer, parameter :: SWP_NOSENDCHANGING =          1024
integer, parameter :: SWP_DEFERERASE =          8192

! SetWindowPos() hwndInsertAfter field values
integer, parameter :: HWND_TOP = (0)
integer, parameter :: HWND_BOTTOM = (1)
integer, parameter :: HWND_TOPMOST = (-1)
integer, parameter :: HWND_NOTOPMOST = (-2)

integer, parameter :: WM_WINDOWPOSCHANGING =            70
integer, parameter :: WM_WINDOWPOSCHANGED =            71

! WM_WINDOWPOSCHANGING/CHANGED struct pointed to by lParam
integer, parameter :: WM_MOVE =             3
integer, parameter :: WM_SIZE =             5

! WM_SIZE message wParam values
integer, parameter :: SIZE_RESTORED = 0
integer, parameter :: SIZE_MINIMIZED = 1
integer, parameter :: SIZE_MAXIMIZED = 2
integer, parameter :: SIZE_MAXSHOW = 3
integer, parameter :: SIZE_MAXHIDE = 4

!***** Main window support *************************************************

integer, parameter :: WM_QUERYOPEN =            19
integer, parameter :: WM_CLOSE =            16

integer, parameter :: WM_GETMINMAXINFO =            36


!***** Window query and enumeration ***************************************

! GetWindow() constants
integer, parameter :: GW_HWNDFIRST = 0
integer, parameter :: GW_HWNDLAST = 1
integer, parameter :: GW_HWNDNEXT = 2
integer, parameter :: GW_HWNDPREV = 3
integer, parameter :: GW_OWNER = 4
integer, parameter :: GW_CHILD = 5

!***** Window drawing support *********************************************

integer, parameter :: DCX_WINDOW =             1
integer, parameter :: DCX_CACHE =             2
integer, parameter :: DCX_CLIPCHILDREN =             8
integer, parameter :: DCX_CLIPSIBLINGS =            16
integer, parameter :: DCX_PARENTCLIP =            32

integer, parameter :: DCX_EXCLUDERGN =            64
integer, parameter :: DCX_INTERSECTRGN =           128


integer, parameter :: DCX_LOCKWINDOWUPDATE =          1024

integer, parameter :: DCX_USESTYLE =         65536

!***** Window repainting **************************************************

integer, parameter :: WM_PAINT =            15
integer, parameter :: WM_ERASEBKGND =            20
integer, parameter :: WM_ICONERASEBKGND =            39

integer, parameter :: RDW_INVALIDATE =             1
integer, parameter :: RDW_INTERNALPAINT =             2
integer, parameter :: RDW_ERASE =             4

integer, parameter :: RDW_VALIDATE =             8
integer, parameter :: RDW_NOINTERNALPAINT =            16
integer, parameter :: RDW_NOERASE =            32

integer, parameter :: RDW_NOCHILDREN =            64
integer, parameter :: RDW_ALLCHILDREN =           128

integer, parameter :: RDW_UPDATENOW =           256
integer, parameter :: RDW_ERASENOW =           512

integer, parameter :: RDW_FRAME =          1024
integer, parameter :: RDW_NOFRAME =          2048


!***** Window scrolling ***************************************************

integer, parameter :: SW_SCROLLCHILDREN =             1
integer, parameter :: SW_INVALIDATE =             2
integer, parameter :: SW_ERASE =             4

!***** Non-client window area management ***********************************

integer, parameter :: WM_NCPAINT =           133

integer, parameter :: WM_NCCALCSIZE =           131

! WM_NCCALCSIZE return flags
integer, parameter :: WVR_ALIGNTOP =            16
integer, parameter :: WVR_ALIGNLEFT =            32
integer, parameter :: WVR_ALIGNBOTTOM =            64
integer, parameter :: WVR_ALIGNRIGHT =           128
integer, parameter :: WVR_HREDRAW =           256
integer, parameter :: WVR_VREDRAW =           512
integer, parameter :: WVR_REDRAW =         ior(WVR_HREDRAW, WVR_VREDRAW)
integer, parameter :: WVR_VALIDRECTS =          1024

integer, parameter :: WM_NCHITTEST =           132

! WM_NCHITTEST return codes
integer, parameter :: HTERROR = (-2)
integer, parameter :: HTTRANSPARENT = (-1)
integer, parameter :: HTNOWHERE = 0
integer, parameter :: HTCLIENT = 1
integer, parameter :: HTCAPTION = 2
integer, parameter :: HTSYSMENU = 3
integer, parameter :: HTSIZE = 4
integer, parameter :: HTMENU = 5
integer, parameter :: HTHSCROLL = 6
integer, parameter :: HTVSCROLL = 7
integer, parameter :: HTMINBUTTON = 8
integer, parameter :: HTMAXBUTTON = 9
integer, parameter :: HTLEFT = 10
integer, parameter :: HTRIGHT = 11
integer, parameter :: HTTOP = 12
integer, parameter :: HTTOPLEFT = 13
integer, parameter :: HTTOPRIGHT = 14
integer, parameter :: HTBOTTOM = 15
integer, parameter :: HTBOTTOMLEFT = 16
integer, parameter :: HTBOTTOMRIGHT = 17
integer, parameter :: HTBORDER = 18
integer, parameter :: HTGROWBOX =          HTSIZE
integer, parameter :: HTREDUCE =           HTMINBUTTON
integer, parameter :: HTZOOM =             HTMAXBUTTON

!***** Drag-and-drop support **********************************************

integer, parameter :: WM_QUERYDRAGICON =            55
integer, parameter :: WM_DROPFILES =           563

!***** Window activation **************************************************

! WM_ACTIVATE state values
integer, parameter :: WA_INACTIVE = 0
integer, parameter :: WA_ACTIVE = 1
integer, parameter :: WA_CLICKACTIVE = 2

integer, parameter :: WM_ACTIVATE =             6
integer, parameter :: WM_ACTIVATEAPP =            28
integer, parameter :: WM_NCACTIVATE =           134

!***** Keyboard input support *********************************************
integer, parameter :: WM_SETFOCUS =             7
integer, parameter :: WM_KILLFOCUS =             8

integer, parameter :: WM_KEYDOWN =           256
integer, parameter :: WM_KEYUP =           257

integer, parameter :: WM_CHAR =           258
integer, parameter :: WM_DEADCHAR =           259

integer, parameter :: WM_SYSKEYDOWN =           260
integer, parameter :: WM_SYSKEYUP =           261

integer, parameter :: WM_SYSCHAR =           262
integer, parameter :: WM_SYSDEADCHAR =           263


! Keyboard message range
integer, parameter :: WM_KEYFIRST =           256
integer, parameter :: WM_KEYLAST =           264

! WM_KEYUP/DOWN/CHAR HIWORD(lParam) flags
integer, parameter :: KF_EXTENDED =           256
integer, parameter :: KF_DLGMODE =          2048
integer, parameter :: KF_MENUMODE =          4096
integer, parameter :: KF_ALTDOWN =          8192
integer, parameter :: KF_REPEAT =         16384
integer, parameter :: KF_UP =         32768

! Virtual key codes
integer, parameter :: VK_LBUTTON =             1
integer, parameter :: VK_RBUTTON =             2
integer, parameter :: VK_CANCEL =             3
integer, parameter :: VK_MBUTTON =             4
integer, parameter :: VK_BACK =             8
integer, parameter :: VK_TAB =             9
integer, parameter :: VK_CLEAR =            12
integer, parameter :: VK_RETURN =            13
integer, parameter :: VK_SHIFT =            16
integer, parameter :: VK_CONTROL =            17
integer, parameter :: VK_MENU =            18
integer, parameter :: VK_PAUSE =            19
integer, parameter :: VK_CAPITAL =            20
integer, parameter :: VK_ESCAPE =            27
integer, parameter :: VK_SPACE =            32
integer, parameter :: VK_PRIOR =            33
integer, parameter :: VK_NEXT =            34
integer, parameter :: VK_END =            35
integer, parameter :: VK_HOME =            36
integer, parameter :: VK_LEFT =            37
integer, parameter :: VK_UP =            38
integer, parameter :: VK_RIGHT =            39
integer, parameter :: VK_DOWN =            40
integer, parameter :: VK_SELECT =            41
integer, parameter :: VK_PRINT =            42
integer, parameter :: VK_EXECUTE =            43
integer, parameter :: VK_SNAPSHOT =            44
integer, parameter :: VK_INSERT =            45
integer, parameter :: VK_DELETE =            46
integer, parameter :: VK_HELP =            47
integer, parameter :: VK_NUMPAD0 =            96
integer, parameter :: VK_NUMPAD1 =            97
integer, parameter :: VK_NUMPAD2 =            98
integer, parameter :: VK_NUMPAD3 =            99
integer, parameter :: VK_NUMPAD4 =           100
integer, parameter :: VK_NUMPAD5 =           101
integer, parameter :: VK_NUMPAD6 =           102
integer, parameter :: VK_NUMPAD7 =           103
integer, parameter :: VK_NUMPAD8 =           104
integer, parameter :: VK_NUMPAD9 =           105
integer, parameter :: VK_MULTIPLY =           106
integer, parameter :: VK_ADD =           107
integer, parameter :: VK_SEPARATOR =           108
integer, parameter :: VK_SUBTRACT =           109
integer, parameter :: VK_DECIMAL =           110
integer, parameter :: VK_DIVIDE =           111
integer, parameter :: VK_F1 =           112
integer, parameter :: VK_F2 =           113
integer, parameter :: VK_F3 =           114
integer, parameter :: VK_F4 =           115
integer, parameter :: VK_F5 =           116
integer, parameter :: VK_F6 =           117
integer, parameter :: VK_F7 =           118
integer, parameter :: VK_F8 =           119
integer, parameter :: VK_F9 =           120
integer, parameter :: VK_F10 =           121
integer, parameter :: VK_F11 =           122
integer, parameter :: VK_F12 =           123
integer, parameter :: VK_F13 =           124
integer, parameter :: VK_F14 =           125
integer, parameter :: VK_F15 =           126
integer, parameter :: VK_F16 =           127
integer, parameter :: VK_F17 =           128
integer, parameter :: VK_F18 =           129
integer, parameter :: VK_F19 =           130
integer, parameter :: VK_F20 =           131
integer, parameter :: VK_F21 =           132
integer, parameter :: VK_F22 =           133
integer, parameter :: VK_F23 =           134
integer, parameter :: VK_F24 =           135
integer, parameter :: VK_NUMLOCK =           144
integer, parameter :: VK_SCROLL =           145

! VK_A thru VK_Z are the same as their ASCII equivalents: 'A' thru 'Z'
! VK_0 thru VK_9 are the same as their ASCII equivalents: '0' thru '0'


! Defines for the fVirt field of the Accelerator table structure.
integer, parameter :: FVIRTKEY  = 1
integer, parameter :: FNOINVERT = 02
integer, parameter :: FSHIFT    = 04
integer, parameter :: FCONTROL  = 08
integer, parameter :: FALT      = 16


! SetWindowsHook() keyboard hook
integer, parameter :: WH_KEYBOARD = 2

!***** Mouse input support ************************************************

! Mouse input messages
integer, parameter :: WM_MOUSEMOVE =           512
integer, parameter :: WM_LBUTTONDOWN =           513
integer, parameter :: WM_LBUTTONUP =           514
integer, parameter :: WM_LBUTTONDBLCLK =           515
integer, parameter :: WM_RBUTTONDOWN =           516
integer, parameter :: WM_RBUTTONUP =           517
integer, parameter :: WM_RBUTTONDBLCLK =           518
integer, parameter :: WM_MBUTTONDOWN =           519
integer, parameter :: WM_MBUTTONUP =           520
integer, parameter :: WM_MBUTTONDBLCLK =           521

! Mouse input message range
integer, parameter :: WM_MOUSEFIRST =           512
integer, parameter :: WM_MOUSELAST =           521

! Mouse message wParam key states
integer, parameter :: MK_LBUTTON =             1
integer, parameter :: MK_RBUTTON =             2
integer, parameter :: MK_SHIFT =             4
integer, parameter :: MK_CONTROL =             8
integer, parameter :: MK_MBUTTON =            16

! Non-client mouse messages
integer, parameter :: WM_NCMOUSEMOVE =           160
integer, parameter :: WM_NCLBUTTONDOWN =           161
integer, parameter :: WM_NCLBUTTONUP =           162
integer, parameter :: WM_NCLBUTTONDBLCLK =           163
integer, parameter :: WM_NCRBUTTONDOWN =           164
integer, parameter :: WM_NCRBUTTONUP =           165
integer, parameter :: WM_NCRBUTTONDBLCLK =           166
integer, parameter :: WM_NCMBUTTONDOWN =           167
integer, parameter :: WM_NCMBUTTONUP =           168
integer, parameter :: WM_NCMBUTTONDBLCLK =           169

! Mouse click activation support
integer, parameter :: WM_MOUSEACTIVATE =            33

! WM_MOUSEACTIVATE return codes
integer, parameter :: MA_ACTIVATE = 1
integer, parameter :: MA_ACTIVATEANDEAT = 2
integer, parameter :: MA_NOACTIVATE = 3
integer, parameter :: MA_NOACTIVATEANDEAT = 4

! SetWindowsHook() mouse hook
integer, parameter :: WH_MOUSE = 7

!***** Mode control *******************************************************

integer, parameter :: WM_CANCELMODE =            31


!***** Timer support ******************************************************

integer, parameter :: WM_TIMER =           275

!***** Menu support *******************************************************

integer, parameter :: MF_INSERT =             0
integer, parameter :: MF_CHANGE =           128
integer, parameter :: MF_APPEND =           256
integer, parameter :: MF_DELETE =           512
integer, parameter :: MF_REMOVE =          4096

! Menu flags for Add/Check/EnableMenuItem()
integer, parameter :: MF_BYCOMMAND =             0
integer, parameter :: MF_BYPOSITION =          1024

integer, parameter :: MF_SEPARATOR =          2048

integer, parameter :: MF_ENABLED =             0
integer, parameter :: MF_GRAYED =             1
integer, parameter :: MF_DISABLED =             2

integer, parameter :: MF_UNCHECKED =             0
integer, parameter :: MF_CHECKED =             8
integer, parameter :: MF_USECHECKBITMAPS =           512

integer, parameter :: MF_STRING =             0
integer, parameter :: MF_BITMAP =             4
integer, parameter :: MF_OWNERDRAW =           256

integer, parameter :: MF_POPUP =            16
integer, parameter :: MF_MENUBARBREAK =            32
integer, parameter :: MF_MENUBREAK =            64

integer, parameter :: MF_UNHILITE =             0
integer, parameter :: MF_HILITE =           128

integer, parameter :: MF_SYSMENU =          8192
integer, parameter :: MF_HELP =         16384
integer, parameter :: MF_MOUSESELECT =         32768


integer, parameter :: MF_END =           128  ! Only valid in menu resource templates

! Flags for TrackPopupMenu
integer, parameter :: TPM_LEFTBUTTON =             0
integer, parameter :: TPM_RIGHTBUTTON =             2
integer, parameter :: TPM_LEFTALIGN =             0
integer, parameter :: TPM_CENTERALIGN =             4
integer, parameter :: TPM_RIGHTALIGN =             8

! Menu messages
integer, parameter :: WM_INITMENU =           278
integer, parameter :: WM_INITMENUPOPUP =           279

integer, parameter :: WM_MENUSELECT =           287
integer, parameter :: WM_MENUCHAR =           288

! Menu and control command messages
integer, parameter :: WM_COMMAND =           273

!***** Scroll bar support *************************************************

integer, parameter :: WM_HSCROLL =           276
integer, parameter :: WM_VSCROLL =           277

! WM_H/VSCROLL commands
integer, parameter :: SB_LINEUP = 0
integer, parameter :: SB_LINELEFT = 0
integer, parameter :: SB_LINEDOWN = 1
integer, parameter :: SB_LINERIGHT = 1
integer, parameter :: SB_PAGEUP = 2
integer, parameter :: SB_PAGELEFT = 2
integer, parameter :: SB_PAGEDOWN = 3
integer, parameter :: SB_PAGERIGHT = 3
integer, parameter :: SB_THUMBPOSITION = 4
integer, parameter :: SB_THUMBTRACK = 5
integer, parameter :: SB_TOP = 6
integer, parameter :: SB_LEFT = 6
integer, parameter :: SB_BOTTOM = 7
integer, parameter :: SB_RIGHT = 7
integer, parameter :: SB_ENDSCROLL = 8

! Scroll bar selection constants
integer, parameter :: SB_HORZ = 0
integer, parameter :: SB_VERT = 1
integer, parameter :: SB_CTL = 2
integer, parameter :: SB_BOTH = 3

! EnableScrollBar() flags
integer, parameter :: ESB_ENABLE_BOTH =             0
integer, parameter :: ESB_DISABLE_BOTH =             3

integer, parameter :: ESB_DISABLE_LEFT =             1
integer, parameter :: ESB_DISABLE_RIGHT =             2

integer, parameter :: ESB_DISABLE_UP =             1
integer, parameter :: ESB_DISABLE_DOWN =             2

integer, parameter :: ESB_DISABLE_LTUP =   ESB_DISABLE_LEFT
integer, parameter :: ESB_DISABLE_RTDN =   ESB_DISABLE_RIGHT


!****** Clipboard manager *************************************************

! Predefined Clipboard Formats
integer, parameter :: CF_TEXT = 1
integer, parameter :: CF_BITMAP = 2
integer, parameter :: CF_METAFILEPICT = 3
integer, parameter :: CF_SYLK = 4
integer, parameter :: CF_DIF = 5
integer, parameter :: CF_TIFF = 6
integer, parameter :: CF_OEMTEXT = 7
integer, parameter :: CF_DIB = 8
integer, parameter :: CF_PALETTE = 9
integer, parameter :: CF_PENDATA = 10
integer, parameter :: CF_RIFF = 11
integer, parameter :: CF_WAVE = 12

integer, parameter :: CF_OWNERDISPLAY =           128
integer, parameter :: CF_DSPTEXT =           129
integer, parameter :: CF_DSPBITMAP =           130
integer, parameter :: CF_DSPMETAFILEPICT =           131

! "Private" formats don't get GlobalFree()'d
integer, parameter :: CF_PRIVATEFIRST =           512
integer, parameter :: CF_PRIVATELAST =           767

! "GDIOBJ" formats do get DeleteObject()'d
integer, parameter :: CF_GDIOBJFIRST =           768
integer, parameter :: CF_GDIOBJLAST =          1023

! Clipboard command messages
integer, parameter :: WM_CUT =           768
integer, parameter :: WM_COPY =           769
integer, parameter :: WM_PASTE =           770
integer, parameter :: WM_CLEAR =           771
integer, parameter :: WM_UNDO =           772

! Clipboard owner messages
integer, parameter :: WM_RENDERFORMAT =           773
integer, parameter :: WM_RENDERALLFORMATS =           774
integer, parameter :: WM_DESTROYCLIPBOARD =           775

! Clipboard viewer messages
integer, parameter :: WM_DRAWCLIPBOARD =           776
integer, parameter :: WM_PAINTCLIPBOARD =           777
integer, parameter :: WM_SIZECLIPBOARD =           779
integer, parameter :: WM_VSCROLLCLIPBOARD =           778
integer, parameter :: WM_HSCROLLCLIPBOARD =           782
integer, parameter :: WM_ASKCBFORMATNAME =           780
integer, parameter :: WM_CHANGECBCHAIN =           781


!***** Mouse cursor support ************************************************

integer, parameter :: WM_SETCURSOR =            32


!***** Message Box support ************************************************

integer, parameter :: MB_OK =                    0
integer, parameter :: MB_OKCANCEL =              1
integer, parameter :: MB_ABORTRETRYIGNORE =      2
integer, parameter :: MB_YESNOCANCEL =           3
integer, parameter :: MB_YESNO =                 4
integer, parameter :: MB_RETRYCANCEL =           5
integer, parameter :: MB_TYPEMASK =             15

integer, parameter :: MB_ICONHAND =             16
integer, parameter :: MB_ICONQUESTION =         32
integer, parameter :: MB_ICONEXCLAMATION =      48
integer, parameter :: MB_ICONASTERISK =         64
integer, parameter :: MB_ICONMASK =            240

integer, parameter :: MB_ICONINFORMATION = MB_ICONASTERISK
integer, parameter :: MB_ICONSTOP =        MB_ICONHAND

integer, parameter :: MB_DEFBUTTON1 =            0
integer, parameter :: MB_DEFBUTTON2 =          256
integer, parameter :: MB_DEFBUTTON3 =          512
integer, parameter :: MB_DEFMASK =            3840

integer, parameter :: MB_APPLMODAL =             0
integer, parameter :: MB_SYSTEMMODAL =        4096
integer, parameter :: MB_TASKMODAL =          8192
integer, parameter :: MB_SETFOREGROUND =     65536

integer, parameter :: MB_NOFOCUS =           32768


!***** WM_SYSCOMMAND support **********************************************

integer, parameter :: WM_SYSCOMMAND =           274

! System Menu Command Values
integer, parameter :: SC_SIZE =         61440
integer, parameter :: SC_MOVE =         61456
integer, parameter :: SC_MINIMIZE =         61472
integer, parameter :: SC_MAXIMIZE =         61488
integer, parameter :: SC_NEXTWINDOW =         61504
integer, parameter :: SC_PREVWINDOW =         61520
integer, parameter :: SC_CLOSE =         61536
integer, parameter :: SC_VSCROLL =         61552
integer, parameter :: SC_HSCROLL =         61568
integer, parameter :: SC_MOUSEMENU =         61584
integer, parameter :: SC_KEYMENU =         61696
integer, parameter :: SC_ARRANGE =         61712
integer, parameter :: SC_RESTORE =         61728
integer, parameter :: SC_TASKLIST =         61744
integer, parameter :: SC_SCREENSAVE =         61760
integer, parameter :: SC_HOTKEY =         61776


!***** MDI Support ********************************************************

! MDI client style bits
integer, parameter :: MDIS_ALLCHILDSTYLES =             1

! MDI messages
integer, parameter :: WM_MDICREATE =           544
integer, parameter :: WM_MDIDESTROY =           545
integer, parameter :: WM_MDIACTIVATE =           546
integer, parameter :: WM_MDIRESTORE =           547
integer, parameter :: WM_MDINEXT =           548
integer, parameter :: WM_MDIMAXIMIZE =           549
integer, parameter :: WM_MDITILE =           550
integer, parameter :: WM_MDICASCADE =           551
integer, parameter :: WM_MDIICONARRANGE =           552
integer, parameter :: WM_MDIGETACTIVE =           553
integer, parameter :: WM_MDISETMENU =           560

! wParam values for WM_MDITILE and WM_MDICASCADE messages.
integer, parameter :: MDITILE_VERTICAL =             0
integer, parameter :: MDITILE_HORIZONTAL =             1
integer, parameter :: MDITILE_SKIPDISABLED =             2

integer, parameter :: WM_CHILDACTIVATE =            34


!***** Dialog and Control Management **************************************

! Dialog window class
integer, parameter :: WC_DIALOG =         32770

! cbWndExtra bytes needed by dialog manager for dialog classes
integer, parameter :: DLGWINDOWEXTRA = 30

! Dialog styles
integer, parameter :: DS_ABSALIGN =             1
integer, parameter :: DS_SYSMODAL =             2
integer, parameter :: DS_LOCALEDIT =            32
integer, parameter :: DS_SETFONT =            64
integer, parameter :: DS_MODALFRAME =           128
integer, parameter :: DS_NOIDLEMSG =           256

! Dialog messages
integer, parameter :: DM_GETDEFID = (WM_USER+0)
integer, parameter :: DM_SETDEFID = (WM_USER+1)

! Returned in HIWORD() of DM_GETDEFID result if msg is supported
integer, parameter :: DC_HASDEFID =         21323

! Dialog notification messages
integer, parameter :: WM_INITDIALOG =           272
integer, parameter :: WM_NEXTDLGCTL =            40

integer, parameter :: WM_PARENTNOTIFY =           528

integer, parameter :: WM_ENTERIDLE =           289

! Get/SetWindowWord/Long offsets for use with WC_DIALOG windows
integer, parameter :: DWL_MSGRESULT = 0
integer, parameter :: DWL_DLGPROC = 4
integer, parameter :: DWL_USER = 8

integer, parameter :: WM_GETDLGCODE =           135

! dialog codes
integer, parameter :: DLGC_WANTARROWS =             1
integer, parameter :: DLGC_WANTTAB =             2
integer, parameter :: DLGC_WANTALLKEYS =             4
integer, parameter :: DLGC_WANTMESSAGE =             4
integer, parameter :: DLGC_HASSETSEL =             8
integer, parameter :: DLGC_DEFPUSHBUTTON =            16
integer, parameter :: DLGC_UNDEFPUSHBUTTON =            32
integer, parameter :: DLGC_RADIOBUTTON =            64
integer, parameter :: DLGC_WANTCHARS =           128
integer, parameter :: DLGC_STATIC =           256
integer, parameter :: DLGC_BUTTON =          8192

integer, parameter :: WM_CTLCOLOR =            25

! WM_CTLCOLOR control IDs
integer, parameter :: CTLCOLOR_MSGBOX = 0
integer, parameter :: CTLCOLOR_EDIT = 1
integer, parameter :: CTLCOLOR_LISTBOX = 2
integer, parameter :: CTLCOLOR_BTN = 3
integer, parameter :: CTLCOLOR_DLG = 4
integer, parameter :: CTLCOLOR_SCROLLBAR = 5
integer, parameter :: CTLCOLOR_STATIC = 6

integer, parameter :: WM_SETFONT =            48
integer, parameter :: WM_GETFONT =            49

! Standard dialog button IDs
integer, parameter :: IDOK = 1
integer, parameter :: IDCANCEL = 2
integer, parameter :: IDABORT = 3
integer, parameter :: IDRETRY = 4
integer, parameter :: IDIGNORE = 5
integer, parameter :: IDYES = 6
integer, parameter :: IDNO = 7

!***** Owner draw control support *****************************************

! Owner draw control types
integer, parameter :: ODT_MENU = 1
integer, parameter :: ODT_LISTBOX = 2
integer, parameter :: ODT_COMBOBOX = 3
integer, parameter :: ODT_BUTTON = 4

! Owner draw actions
integer, parameter :: ODA_DRAWENTIRE =             1
integer, parameter :: ODA_SELECT =             2
integer, parameter :: ODA_FOCUS =             4

! Owner draw state
integer, parameter :: ODS_SELECTED =             1
integer, parameter :: ODS_GRAYED =             2
integer, parameter :: ODS_DISABLED =             4
integer, parameter :: ODS_CHECKED =             8
integer, parameter :: ODS_FOCUS =            16

integer, parameter :: WM_DRAWITEM =            43

integer, parameter :: WM_MEASUREITEM =            44

integer, parameter :: WM_DELETEITEM =            45

integer, parameter :: WM_COMPAREITEM =            57

!***** Static control *****************************************************

! Static Control Styles
integer, parameter :: SS_LEFT =             0
integer, parameter :: SS_CENTER =             1
integer, parameter :: SS_RIGHT =             2
integer, parameter :: SS_ICON =             3
integer, parameter :: SS_BLACKRECT =             4
integer, parameter :: SS_GRAYRECT =             5
integer, parameter :: SS_WHITERECT =             6
integer, parameter :: SS_BLACKFRAME =             7
integer, parameter :: SS_GRAYFRAME =             8
integer, parameter :: SS_WHITEFRAME =             9
integer, parameter :: SS_SIMPLE =            11
integer, parameter :: SS_LEFTNOWORDWRAP =            12
integer, parameter :: SS_NOPREFIX =           128

! Static Control Mesages
integer, parameter :: STM_SETICON = (WM_USER+0)
integer, parameter :: STM_GETICON = (WM_USER+1)

!***** Button control ****************************************************

! Button Control Styles
integer, parameter :: BS_PUSHBUTTON =             0
integer, parameter :: BS_DEFPUSHBUTTON =             1
integer, parameter :: BS_CHECKBOX =             2
integer, parameter :: BS_AUTOCHECKBOX =             3
integer, parameter :: BS_RADIOBUTTON =             4
integer, parameter :: BS_3STATE =             5
integer, parameter :: BS_AUTO3STATE =             6
integer, parameter :: BS_GROUPBOX =             7
integer, parameter :: BS_USERBUTTON =             8
integer, parameter :: BS_AUTORADIOBUTTON =             9
integer, parameter :: BS_OWNERDRAW =            11
integer, parameter :: BS_LEFTTEXT =            32

! Button Control Messages
integer, parameter :: BM_GETCHECK = (WM_USER+0)
integer, parameter :: BM_SETCHECK = (WM_USER+1)
integer, parameter :: BM_GETSTATE = (WM_USER+2)
integer, parameter :: BM_SETSTATE = (WM_USER+3)
integer, parameter :: BM_SETSTYLE = (WM_USER+4)

! User Button Notification Codes
integer, parameter :: BN_CLICKED = 0
integer, parameter :: BN_PAINT = 1
integer, parameter :: BN_HILITE = 2
integer, parameter :: BN_UNHILITE = 3
integer, parameter :: BN_DISABLE = 4
integer, parameter :: BN_DOUBLECLICKED = 5


!***** Edit control ******************************************************


! Edit control styles
integer, parameter :: ES_LEFT =             0
integer, parameter :: ES_CENTER =             1
integer, parameter :: ES_RIGHT =             2
integer, parameter :: ES_MULTILINE =             4
integer, parameter :: ES_UPPERCASE =             8
integer, parameter :: ES_LOWERCASE =            16
integer, parameter :: ES_PASSWORD =            32
integer, parameter :: ES_AUTOVSCROLL =            64
integer, parameter :: ES_AUTOHSCROLL =           128
integer, parameter :: ES_NOHIDESEL =           256
integer, parameter :: ES_OEMCONVERT =          1024
integer, parameter :: ES_READONLY =          2048
integer, parameter :: ES_WANTRETURN =          4096

! Edit control messages
integer, parameter :: EM_GETSEL = (WM_USER+0)
integer, parameter :: EM_SETSEL = (WM_USER+1)
integer, parameter :: EM_GETRECT = (WM_USER+2)
integer, parameter :: EM_SETRECT = (WM_USER+3)
integer, parameter :: EM_SETRECTNP = (WM_USER+4)
integer, parameter :: EM_LINESCROLL = (WM_USER+6)
integer, parameter :: EM_GETMODIFY = (WM_USER+8)
integer, parameter :: EM_SETMODIFY = (WM_USER+9)
integer, parameter :: EM_GETLINECOUNT = (WM_USER+10)
integer, parameter :: EM_LINEINDEX = (WM_USER+11)
integer, parameter :: EM_SETHANDLE = (WM_USER+12)
integer, parameter :: EM_GETHANDLE = (WM_USER+13)
integer, parameter :: EM_LINELENGTH = (WM_USER+17)
integer, parameter :: EM_REPLACESEL = (WM_USER+18)
integer, parameter :: EM_SETFONT = (WM_USER+19)    ! NOT IMPLEMENTED: use WM_SETFONT
integer, parameter :: EM_GETLINE = (WM_USER+20)
integer, parameter :: EM_LIMITTEXT = (WM_USER+21)
integer, parameter :: EM_CANUNDO = (WM_USER+22)
integer, parameter :: EM_UNDO = (WM_USER+23)
integer, parameter :: EM_FMTLINES = (WM_USER+24)
integer, parameter :: EM_LINEFROMCHAR = (WM_USER+25)
integer, parameter :: EM_SETWORDBREAK = (WM_USER+26)    ! NOT IMPLEMENTED: use EM_SETWORDBREAK
integer, parameter :: EM_SETTABSTOPS = (WM_USER+27)
integer, parameter :: EM_SETPASSWORDCHAR = (WM_USER+28)
integer, parameter :: EM_EMPTYUNDOBUFFER = (WM_USER+29)
integer, parameter :: EM_GETFIRSTVISIBLELINE = (WM_USER+30)
integer, parameter :: EM_SETREADONLY = (WM_USER+31)
integer, parameter :: EM_SETWORDBREAKPROC = (WM_USER+32)
integer, parameter :: EM_GETWORDBREAKPROC = (WM_USER+33)
integer, parameter :: EM_GETPASSWORDCHAR = (WM_USER+34)

! EDITWORDBREAKPROC code values
integer, parameter :: WB_LEFT = 0
integer, parameter :: WB_RIGHT = 1
integer, parameter :: WB_ISDELIMITER = 2

! Edit control notification codes
integer, parameter :: EN_SETFOCUS =           256
integer, parameter :: EN_KILLFOCUS =           512
integer, parameter :: EN_CHANGE =           768
integer, parameter :: EN_UPDATE =          1024
integer, parameter :: EN_ERRSPACE =          1280
integer, parameter :: EN_MAXTEXT =          1281
integer, parameter :: EN_HSCROLL =          1537
integer, parameter :: EN_VSCROLL =          1538


!***** Scroll bar control ************************************************
! Also see scrolling support

! Scroll bar styles
integer, parameter :: SBS_HORZ =             0
integer, parameter :: SBS_VERT =             1
integer, parameter :: SBS_TOPALIGN =             2
integer, parameter :: SBS_LEFTALIGN =             2
integer, parameter :: SBS_BOTTOMALIGN =             4
integer, parameter :: SBS_RIGHTALIGN =             4
integer, parameter :: SBS_SIZEBOXTOPLEFTALIGN =             2
integer, parameter :: SBS_SIZEBOXBOTTOMRIGHTALIGN =             4
integer, parameter :: SBS_SIZEBOX =             8

!***** Listbox control ***************************************************

! Listbox styles
integer, parameter :: LBS_NOTIFY =             1
integer, parameter :: LBS_SORT =             2
integer, parameter :: LBS_NOREDRAW =             4
integer, parameter :: LBS_MULTIPLESEL =             8
integer, parameter :: LBS_OWNERDRAWFIXED =            16
integer, parameter :: LBS_OWNERDRAWVARIABLE =            32
integer, parameter :: LBS_HASSTRINGS =            64
integer, parameter :: LBS_USETABSTOPS =           128
integer, parameter :: LBS_NOINTEGRALHEIGHT =           256
integer, parameter :: LBS_MULTICOLUMN =           512
integer, parameter :: LBS_WANTKEYBOARDINPUT =          1024
integer, parameter :: LBS_EXTENDEDSEL =          2048
integer, parameter :: LBS_DISABLENOSCROLL =          4096
integer, parameter :: LBS_STANDARD = ior(LBS_NOTIFY, ior(LBS_SORT, ior(WS_VSCROLL, WS_BORDER)))

! Listbox messages
integer, parameter :: LB_ADDSTRING = 384
integer, parameter :: LB_INSERTSTRING = 385
integer, parameter :: LB_DELETESTRING = 386
integer, parameter :: LB_SELITEMRANGEEX = 387
integer, parameter :: LB_RESETCONTENT = 388
integer, parameter :: LB_SETSEL = 389
integer, parameter :: LB_SETCURSEL = 390
integer, parameter :: LB_GETSEL = 391
integer, parameter :: LB_GETCURSEL = 392
integer, parameter :: LB_GETTEXT = 393
integer, parameter :: LB_GETTEXTLEN = 394
integer, parameter :: LB_GETCOUNT = 395
integer, parameter :: LB_SELECTSTRING = 396
integer, parameter :: LB_DIR = 397
integer, parameter :: LB_GETTOPINDEX = 398
integer, parameter :: LB_FINDSTRING = 399
integer, parameter :: LB_GETSELCOUNT = 400
integer, parameter :: LB_GETSELITEMS = 401
integer, parameter :: LB_SETTABSTOPS = 402
integer, parameter :: LB_GETHORIZONTALEXTENT = 403
integer, parameter :: LB_SETHORIZONTALEXTENT = 404
integer, parameter :: LB_SETCOLUMNWIDTH = 405
integer, parameter :: LB_ADDFILE = 406
integer, parameter :: LB_SETTOPINDEX = 407
integer, parameter :: LB_GETITEMRECT = 408
integer, parameter :: LB_GETITEMDATA = 409
integer, parameter :: LB_SETITEMDATA = 410
integer, parameter :: LB_SELITEMRANGE = 411
integer, parameter :: LB_SETANCHORINDEX = 412
integer, parameter :: LB_GETANCHORINDEX = 413
integer, parameter :: LB_SETCARETINDEX = 414
integer, parameter :: LB_GETCARETINDEX = 415
integer, parameter :: LB_SETITEMHEIGHT = 416
integer, parameter :: LB_GETITEMHEIGHT = 417
integer, parameter :: LB_FINDSTRINGEXACT = 418
integer, parameter :: LB_SETLOCALE = 421
integer, parameter :: LB_GETLOCALE = 422
integer, parameter :: LB_SETCOUNT = 423
integer, parameter :: LB_INITSTORAGE = 424
integer, parameter :: LB_ITEMFROMPOINT = 425

! Listbox notification codes
integer, parameter :: LBN_ERRSPACE = (-2)
integer, parameter :: LBN_SELCHANGE = 1
integer, parameter :: LBN_DBLCLK = 2
integer, parameter :: LBN_SELCANCEL = 3
integer, parameter :: LBN_SETFOCUS = 4
integer, parameter :: LBN_KILLFOCUS = 5

! Listbox notification messages
integer, parameter :: WM_VKEYTOITEM =            46
integer, parameter :: WM_CHARTOITEM =            47

! Listbox message return values
integer, parameter :: LB_OKAY = 0
integer, parameter :: LB_ERR = (-1)
integer, parameter :: LB_ERRSPACE = (-2)

integer, parameter :: LB_CTLCODE = 0

!***** Dialog directory support *******************************************

! DlgDirList, DlgDirListComboBox flags values
integer, parameter :: DDL_READWRITE =             0
integer, parameter :: DDL_READONLY =             1
integer, parameter :: DDL_HIDDEN =             2
integer, parameter :: DDL_SYSTEM =             4
integer, parameter :: DDL_DIRECTORY =            16
integer, parameter :: DDL_ARCHIVE =            32

integer, parameter :: DDL_POSTMSGS =          8192
integer, parameter :: DDL_DRIVES =         16384
integer, parameter :: DDL_EXCLUSIVE =         32768

!***** Combo box control *************************************************

! Combo box styles
integer, parameter :: CBS_SIMPLE =             1
integer, parameter :: CBS_DROPDOWN =             2
integer, parameter :: CBS_DROPDOWNLIST =             3
integer, parameter :: CBS_OWNERDRAWFIXED =            16
integer, parameter :: CBS_OWNERDRAWVARIABLE =            32
integer, parameter :: CBS_AUTOHSCROLL =            64
integer, parameter :: CBS_OEMCONVERT =           128
integer, parameter :: CBS_SORT =           256
integer, parameter :: CBS_HASSTRINGS =           512
integer, parameter :: CBS_NOINTEGRALHEIGHT =          1024
integer, parameter :: CBS_DISABLENOSCROLL =          2048

! Combo box messages
integer, parameter :: CB_GETEDITSEL = (WM_USER+0)
integer, parameter :: CB_LIMITTEXT = (WM_USER+1)
integer, parameter :: CB_SETEDITSEL = (WM_USER+2)
integer, parameter :: CB_ADDSTRING = (WM_USER+3)
integer, parameter :: CB_DELETESTRING = (WM_USER+4)
integer, parameter :: CB_DIR = (WM_USER+5)
integer, parameter :: CB_GETCOUNT = (WM_USER+6)
integer, parameter :: CB_GETCURSEL = (WM_USER+7)
integer, parameter :: CB_GETLBTEXT = (WM_USER+8)
integer, parameter :: CB_GETLBTEXTLEN = (WM_USER+9)
integer, parameter :: CB_INSERTSTRING = (WM_USER+10)
integer, parameter :: CB_RESETCONTENT = (WM_USER+11)
integer, parameter :: CB_FINDSTRING = (WM_USER+12)
integer, parameter :: CB_SELECTSTRING = (WM_USER+13)
integer, parameter :: CB_SETCURSEL = (WM_USER+14)
integer, parameter :: CB_SHOWDROPDOWN = (WM_USER+15)
integer, parameter :: CB_GETITEMDATA = (WM_USER+16)
integer, parameter :: CB_SETITEMDATA = (WM_USER+17)
integer, parameter :: CB_GETDROPPEDCONTROLRECT = (WM_USER+18)
integer, parameter :: CB_SETITEMHEIGHT = (WM_USER+19)
integer, parameter :: CB_GETITEMHEIGHT = (WM_USER+20)
integer, parameter :: CB_SETEXTENDEDUI = (WM_USER+21)
integer, parameter :: CB_GETEXTENDEDUI = (WM_USER+22)
integer, parameter :: CB_GETDROPPEDSTATE = (WM_USER+23)
integer, parameter :: CB_FINDSTRINGEXACT = (WM_USER+24)

! Combo box notification codes
integer, parameter :: CBN_ERRSPACE = (-1)
integer, parameter :: CBN_SELCHANGE = 1
integer, parameter :: CBN_DBLCLK = 2
integer, parameter :: CBN_SETFOCUS = 3
integer, parameter :: CBN_KILLFOCUS = 4
integer, parameter :: CBN_EDITCHANGE = 5
integer, parameter :: CBN_EDITUPDATE = 6
integer, parameter :: CBN_DROPDOWN = 7
integer, parameter :: CBN_CLOSEUP = 8
integer, parameter :: CBN_SELENDOK = 9
integer, parameter :: CBN_SELENDCANCEL = 10

! Combo box message return values
integer, parameter :: CB_OKAY = 0
integer, parameter :: CB_ERR = (-1)
integer, parameter :: CB_ERRSPACE = (-2)


!****** Windows hook support *********************************************

! Standard hook code
integer, parameter :: HC_ACTION = 0


!***** Computer-based-training (CBT) support ******************************

integer, parameter :: WM_QUEUESYNC =            35

! SetWindowsHook() code
integer, parameter :: WH_CBT = 5

integer, parameter :: HCBT_MOVESIZE = 0
integer, parameter :: HCBT_MINMAX = 1
integer, parameter :: HCBT_QS = 2
integer, parameter :: HCBT_CREATEWND = 3
integer, parameter :: HCBT_DESTROYWND = 4
integer, parameter :: HCBT_ACTIVATE = 5
integer, parameter :: HCBT_CLICKSKIPPED = 6
integer, parameter :: HCBT_KEYSKIPPED = 7
integer, parameter :: HCBT_SYSCOMMAND = 8
integer, parameter :: HCBT_SETFOCUS = 9

!***** Hardware hook support **********************************************

integer, parameter :: WH_HARDWARE = 8

!***** Shell support ******************************************************

! SetWindowsHook() Shell hook code
integer, parameter :: WH_SHELL = 10

integer, parameter :: HSHELL_WINDOWCREATED = 1
integer, parameter :: HSHELL_WINDOWDESTROYED = 2
integer, parameter :: HSHELL_ACTIVATESHELLWINDOW = 3


!***** Journalling support ************************************************

integer, parameter :: WH_JOURNALRECORD = 0
integer, parameter :: WH_JOURNALPLAYBACK = 1

! Journalling hook codes
integer, parameter :: HC_GETNEXT = 1
integer, parameter :: HC_SKIP = 2
integer, parameter :: HC_NOREMOVE = 3
integer, parameter :: HC_NOREM =           HC_NOREMOVE
integer, parameter :: HC_SYSMODALON = 4
integer, parameter :: HC_SYSMODALOFF = 5


!***** Debugger support ***************************************************

! SetWindowsHook debug hook support
integer, parameter :: WH_DEBUG = 9

! Flags returned by GetSystemDebugState.

integer, parameter :: SDS_MENU =             1
integer, parameter :: SDS_SYSMODAL =             2
integer, parameter :: SDS_NOTASKQUEUE =             4
integer, parameter :: SDS_DIALOG =             8
integer, parameter :: SDS_TASKLOCKED =            16

!***** Help support *******************************************************

! WinHelp() commands
integer, parameter :: HELP_CONTEXT =             1
integer, parameter :: HELP_QUIT =             2
integer, parameter :: HELP_INDEX =             3
integer, parameter :: HELP_CONTENTS =             3
integer, parameter :: HELP_HELPONHELP =             4
integer, parameter :: HELP_SETINDEX =             5
integer, parameter :: HELP_SETCONTENTS =             5
integer, parameter :: HELP_CONTEXTPOPUP =             8
integer, parameter :: HELP_FORCEFILE =             9
integer, parameter :: HELP_KEY =           257
integer, parameter :: HELP_COMMAND =           258
integer, parameter :: HELP_PARTIALKEY =           261
integer, parameter :: HELP_MULTIKEY =           513
integer, parameter :: HELP_SETWINPOS =           515

!***** Sound support *****************************************************

! SetSoundNoise() Sources
integer, parameter :: S_PERIOD512 = 0
integer, parameter :: S_PERIOD1024 = 1
integer, parameter :: S_PERIOD2048 = 2
integer, parameter :: S_PERIODVOICE = 3
integer, parameter :: S_WHITE512 = 4
integer, parameter :: S_WHITE1024 = 5
integer, parameter :: S_WHITE2048 = 6
integer, parameter :: S_WHITEVOICE = 7

! WaitSoundState() constants
integer, parameter :: S_QUEUEEMPTY = 0
integer, parameter :: S_THRESHOLD = 1
integer, parameter :: S_ALLTHRESHOLD = 2

! Accent Modes
integer, parameter :: S_NORMAL = 0
integer, parameter :: S_LEGATO = 1
integer, parameter :: S_STACCATO = 2

! Error return values
integer, parameter :: S_SERDVNA = (-1)
integer, parameter :: S_SEROFM = (-2)
integer, parameter :: S_SERMACT = (-3)
integer, parameter :: S_SERQFUL = (-4)
integer, parameter :: S_SERBDNT = (-5)
integer, parameter :: S_SERDLN = (-6)
integer, parameter :: S_SERDCC = (-7)
integer, parameter :: S_SERDTP = (-8)
integer, parameter :: S_SERDVL = (-9)
integer, parameter :: S_SERDMD = (-10)
integer, parameter :: S_SERDSH = (-11)
integer, parameter :: S_SERDPT = (-12)
integer, parameter :: S_SERDFQ = (-13)
integer, parameter :: S_SERDDR = (-14)
integer, parameter :: S_SERDSR = (-15)
integer, parameter :: S_SERDST = (-16)


!***** Comm support *****************************************************

integer, parameter :: NOPARITY = 0
integer, parameter :: ODDPARITY = 1
integer, parameter :: EVENPARITY = 2
integer, parameter :: MARKPARITY = 3
integer, parameter :: SPACEPARITY = 4

integer, parameter :: ONESTOPBIT = 0
integer, parameter :: ONE5STOPBITS = 1
integer, parameter :: TWOSTOPBITS = 2

integer, parameter :: IGNORE = 0
integer, parameter :: INFINITE =         65535

! Error Flags
integer, parameter :: CE_RXOVER =             1
integer, parameter :: CE_OVERRUN =             2
integer, parameter :: CE_RXPARITY =             4
integer, parameter :: CE_FRAME =             8
integer, parameter :: CE_BREAK =            16
integer, parameter :: CE_CTSTO =            32
integer, parameter :: CE_DSRTO =            64
integer, parameter :: CE_RLSDTO =           128
integer, parameter :: CE_TXFULL =           256
integer, parameter :: CE_PTO =           512
integer, parameter :: CE_IOE =          1024
integer, parameter :: CE_DNS =          2048
integer, parameter :: CE_OOP =          4096
integer, parameter :: CE_MODE =         32768

integer, parameter :: IE_BADID = (-1)
integer, parameter :: IE_OPEN = (-2)
integer, parameter :: IE_NOPEN = (-3)
integer, parameter :: IE_MEMORY = (-4)
integer, parameter :: IE_DEFAULT = (-5)
integer, parameter :: IE_HARDWARE = (-10)
integer, parameter :: IE_BYTESIZE = (-11)
integer, parameter :: IE_BAUDRATE = (-12)

! Events
integer, parameter :: EV_RXCHAR =             1
integer, parameter :: EV_RXFLAG =             2
integer, parameter :: EV_TXEMPTY =             4
integer, parameter :: EV_CTS =             8
integer, parameter :: EV_DSR =            16
integer, parameter :: EV_RLSD =            32
integer, parameter :: EV_BREAK =            64
integer, parameter :: EV_ERR =           128
integer, parameter :: EV_RING =           256
integer, parameter :: EV_PERR =           512
integer, parameter :: EV_CTSS =          1024
integer, parameter :: EV_DSRS =          2048
integer, parameter :: EV_RLSDS =          4096
integer, parameter :: EV_RINGTE =          8192

! Escape Functions
integer, parameter :: SETXOFF = 1
integer, parameter :: SETXON = 2
integer, parameter :: SETRTS = 3
integer, parameter :: CLRRTS = 4
integer, parameter :: SETDTR = 5
integer, parameter :: CLRDTR = 6
integer, parameter :: RESETDEV = 7

integer, parameter :: LPTx =           128

! new escape functions
integer, parameter :: GETMAXLPT = 8
integer, parameter :: GETMAXCOM = 9
integer, parameter :: GETBASEIRQ = 10

! Comm Baud Rate indices
integer, parameter :: CBR_110 =         65296
integer, parameter :: CBR_300 =         65297
integer, parameter :: CBR_600 =         65298
integer, parameter :: CBR_1200 =         65299
integer, parameter :: CBR_2400 =         65300
integer, parameter :: CBR_4800 =         65301
integer, parameter :: CBR_9600 =         65302
integer, parameter :: CBR_14400 =         65303
integer, parameter :: CBR_19200 =         65304
integer, parameter :: CBR_38400 =         65307
integer, parameter :: CBR_56000 =         65311
integer, parameter :: CBR_128000 =         65315
integer, parameter :: CBR_256000 =         65319

! notifications passed in low word of lParam on WM_COMMNOTIFY messages
integer, parameter :: CN_RECEIVE =             1
integer, parameter :: CN_TRANSMIT =             2
integer, parameter :: CN_EVENT =             4

integer, parameter :: CSTF_CTSHOLD =             1
integer, parameter :: CSTF_DSRHOLD =             2
integer, parameter :: CSTF_RLSDHOLD =             4
integer, parameter :: CSTF_XOFFHOLD =             8
integer, parameter :: CSTF_XOFFSENT =            16
integer, parameter :: CSTF_EOF =            32
integer, parameter :: CSTF_TXIM =            64

integer, parameter :: WM_COMMNOTIFY =            68

!***** Driver support *****************************************************

! Driver messages
integer, parameter :: DRV_LOAD =             1
integer, parameter :: DRV_ENABLE =             2
integer, parameter :: DRV_OPEN =             3
integer, parameter :: DRV_CLOSE =             4
integer, parameter :: DRV_DISABLE =             5
integer, parameter :: DRV_FREE =             6
integer, parameter :: DRV_CONFIGURE =             7
integer, parameter :: DRV_QUERYCONFIGURE =             8
integer, parameter :: DRV_INSTALL =             9
integer, parameter :: DRV_REMOVE =            10
integer, parameter :: DRV_EXITSESSION =            11
integer, parameter :: DRV_EXITAPPLICATION =            12
integer, parameter :: DRV_POWER =            15

integer, parameter :: DRV_RESERVED =          2048
integer, parameter :: DRV_USER =         16384

! Supported return values for DRV_CONFIGURE message
integer, parameter :: DRVCNF_CANCEL =             0
integer, parameter :: DRVCNF_OK =             1
integer, parameter :: DRVCNF_RESTART =             2

! Supported lParam1 of DRV_EXITAPPLICATION notification
integer, parameter :: DRVEA_NORMALEXIT =             1
integer, parameter :: DRVEA_ABNORMALEXIT =             2

! GetNextDriver flags
integer, parameter :: GND_FIRSTINSTANCEONLY =             1

integer, parameter :: GND_FORWARD =             0
integer, parameter :: GND_REVERSE =             2

end module

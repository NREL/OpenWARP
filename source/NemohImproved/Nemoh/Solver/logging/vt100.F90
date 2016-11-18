!> This module sets terminal colors, boldness and other settings
!! Using ANSI/VT100 control sequences
!! See http://misc.flogisoft.com/bash/tip_colors_and_formatting for a list of sequences
!! This code is governed by the X11 license. See LICENSE for details.
module vt100
  implicit none
  ! Control start character
  character(len=*), parameter :: start = achar(27)
  character(len=*), parameter :: reset = "0"
  ! Styles
  character(len=*), parameter :: bold = "1", dimmed = "2", &
      underline = "4", blink = "5", invert = "7", hidden = "8"
  contains
    subroutine tput(lu, code)
      implicit none
      character(len=*), intent(in) :: code
      integer, intent(in) :: lu
      write(lu, '(a,"[",a,"m")', advance="no") start, code
    end subroutine tput
    subroutine stput(str, code)
      implicit none
      character(len=*), intent(inout) :: str
      character(len=*), intent(in)    :: code
      str = trim(str) // start // "[" // trim(code) // "m"
    end subroutine stput
end module vt100

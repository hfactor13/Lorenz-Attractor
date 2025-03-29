module linspace_module
! This module provides linearly spaced data points
implicit none
public :: linspace

contains

function linspace(start, end, len) result(x)
    real, intent(in) :: start, end
    integer, intent(in) :: len
    real :: dx
    integer :: i
    real, dimension(1:len) :: x


    dx = (end - start) / (len - 1)
    x(1:len) = [(start + (i - 1) * dx, i = 1, len)]


end function

end module
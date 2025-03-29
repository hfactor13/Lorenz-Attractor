program lorenz_attractor
    ! The purpose of this script is to calculate the state variables of the Lorenz attractor using the Runge-Kutta method order 4 (RK4)
    use linspace_module
    implicit none

    ! Initialize parameters
    real, parameter :: rho = 28.0
    real, parameter :: sigma = 10.0
    real, parameter :: beta = 8.0 / 3.0
    real, parameter :: T = 10.0
    real, parameter :: dt = 0.01
    integer, parameter :: n = int(T / dt)
    
    ! Declare variables
    real, dimension(n, 3) :: y
    integer :: i

    ! Generate time vector
    real, dimension(n) :: time
    time = linspace(0.0, T, n)

    ! Open file for writing data
    open(unit=1, file="./output/lorenz_data.dat", status = "replace")

    ! Provide inital guess
    y(1,:) = [-8.0, 8.0, 27.0]

    ! Apply RK4 algorithm at each step
    do i = 1, n-1
        y(i+1,:) = rk4_singlestep(lorenz, dt, time(i), y(i,:))
    end do

    ! Write results to a file
    do i = 1, n
        write(1, *) y(i, 1), y(i, 2), y(i, 3)
    end do

    close(1)

contains

    function lorenz(t, y) result(dy)
        real, intent(in) :: t
        real, dimension(3), intent(in) :: y
        real, dimension(3) :: dy

        dy(1) = sigma * (y(2) - y(1))
        dy(2) = y(1) * (rho - y(3)) - y(2)
        dy(3) = y(1) * y(2) - beta * y(3)
    end function

    function rk4_singlestep(fun, dt, t0, y0) result(yout)
        interface
            function fun(t, y) result(dy)
                real, intent(in) :: t
                real, dimension(3), intent(in) :: y
                real, dimension(3) :: dy
            end function
        end interface
        real, intent(in) :: dt, t0
        real, dimension(3), intent(in) :: y0
        real, dimension(3) :: k1, k2, k3, k4
        real, dimension(3) :: yout

        k1 = fun(t0, y0)
        k2 = fun(t0 + dt / 2, y0 + dt / 2 * k1)
        k3 = fun(t0 + dt / 2, y0 + dt / 2 * k2)
        k4 = fun(t0 + dt, y0 + dt * k3)
        yout = y0 + dt / 6 * (k1 + 2*k2 + 2*k3 + k4)

    end function

end program

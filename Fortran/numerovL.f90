module Constants
! This module contains all the physical constants (it define your world)
real(8)               ::  mass
real(8)               ::  hbaron2m
real(8),    parameter ::  pi       = 3.141592653589793238462643383279
logical,    parameter ::  debug    = .true.
real(8),save          :: LEC(10)
real(8),save          :: LEC2(10)
real(8),save          :: test
real(8), allocatable  :: Identity(:)   
end module Constants





module solvers
use Constants
  implicit none

! -----------------
!
! This module contains all the solvers I need
!
! Distance converts a integer to a distance
! f is an auxilliary function for numerov
! der_7p performs derivatives using 7 points
! simpson performs the simpson integration
! NLO correction apply the second order of EFT potentials
! 

contains



function Distance(Step,i) result (D)
! ---------------------------- !
! Step is a real number that   !
! indicates the spatial        !
! resolution                   !
! i is the index of the spatial!
! vector                       !
! ---------------------------- !
  real(8)                 :: D
  real(8),    intent(in)  :: Step
  integer, intent(in)     :: i
  D= (i-1)*Step
end function Distance


   function f(Pot, E,  Step, i) result(v)
! ---------------------------- !
! f is used by numerov to get  !
! the next wave function point !
! Pot is the pointer to the    !
! potential function, E is the !
! energy, step the spatial res !
! and i the index of the acual !
! point.                       !
! ---------------------------- !
    interface
        function Pot(r) result(v)
            real(8), intent (in)  :: r
            real(8)              :: v
        end function Pot
    end interface

      integer, intent(in)  :: i
      real(8), intent(in)  :: E, Step
      real(8)              :: v
        v = (Pot(Distance(Step,i))-E)/hbaron2m
   end function f



subroutine der_7p(xr,f,fp,fpp)
! ---------------------------- !
! xr : Vectror of coordinates  !
! ---------------------------- !
   integer              :: i,nr,nstep
   real(8)              :: hr
   real(8)              :: xr(:),f(:),fp(:),fpp(:)
   
   nstep = size(xr)
   nr    = size(xr)
   hr=xr(2)-xr(1)
! 1st derivative
   fp(1)=(-10.0 *f(7)+72.0 *f(6)-225.0 *f(5)+400.0 *f(4)-450.0 *f(3)+360.0 *f(2)-147.0 *f(1))/(60.0 *hr)
   fp(2)=(  2.0 *f(7)-15.0 *f(6)+ 50.0 *f(5)-100.0 *f(4)+150.0 *f(3)- 77.0 *f(2)- 10.0 *f(1))/(60.0 *hr)
   fp(3)=(- 1.0 *f(7)+ 8.0 *f(6)- 30.0 *f(5)+ 80.0 *f(4)- 35.0 *f(3)- 24.0 *f(2)+  2.0 *f(1))/(60.0 *hr)
   do i=4,nstep-3
      fp(i)=(f(i+3)-9.0 *f(i+2)+45.0 *f(i+1)-45.0 *f(i-1)+9.0 *f(i-2)-f(i-3))/(60.0 *hr)
   enddo
   fp(nstep-2)=(-  2.0 *f(nstep)+ 24.0 *f(nstep-1)+ 35.0 *f(nstep-2)- 80.0 *f(nstep-3)+ 30.0 *f(nstep-4) &
   & -  8.0 *f(nstep-5)+  1.0 *f(nstep-6))/(60.0 *hr)
   fp(nstep-1)=(  10.0 *f(nstep)+ 77.0 *f(nstep-1)-150.0 *f(nstep-2)+100.0 *f(nstep-3)- 50.0 *f(nstep-4) &
   & + 15.0 *f(nstep-5)-  2.0 *f(nstep-6))/(60.0 *hr)
   fp(nstep)  =( 147.0 *f(nstep)-360.0 *f(nstep-1)+450.0 *f(nstep-2)-400.0 *f(nstep-3)+225.0 *f(nstep-4) &
   & - 72.0 *f(nstep-5)+ 10.0 *f(nstep-6))/(60.0 *hr)
! 2nd derivative
   fpp(1)=( 137.0 *f(7)-972.0 *f(6)+2970.0 *f(5)-5080.0 *f(4)+5265.0 *f(3)-3132.0 *f(2)+812.0 *f(1))/(180.0 *hr**2)
   fpp(2)=(- 13.0 *f(7)+ 93.0 *f(6)- 285.0 *f(5)+ 470.0 *f(4)- 255.0 *f(3)- 147.0 *f(2)+137.0 *f(1))/(180.0 *hr**2)
   fpp(3)=(   2.0 *f(7)- 12.0 *f(6)+  15.0 *f(5)+ 200.0 *f(4)- 420.0 *f(3)+ 228.0 *f(2)- 13.0 *f(1))/(180.0 *hr**2)
   do i=4,nr-3
      fpp(i)=(2.0 *f(i+3)-27.0 *f(i+2)+270.0 *f(i+1)-490 *f(i)+270.0 *f(i-1)-27.0 *f(i-2)+2.0 *f(i-3))/(180.0 *hr**2)
   enddo
   fpp(nstep-2)=(- 13.0 *f(nstep)+ 228.0 *f(nstep-1)- 420.0 *f(nstep-2)+ 200.0 *f(nstep-3)+  15.0 *f(nstep-4) &
   & - 12.0 *f(nstep-5)+  2.0 *f(nstep-6))/(180.0 *hr**2)
   fpp(nstep-1)=( 137.0 *f(nstep)- 147.0 *f(nstep-1)- 255.0 *f(nstep-2)+ 470.0 *f(nstep-3)- 285.0 *f(nstep-4)& 
   & + 93.0 *f(nstep-5)- 13.0 *f(nstep-6))/(180.0 *hr**2)
   fpp(nstep)  =( 812.0 *f(nstep)-3132.0 *f(nstep-1)+5265.0 *f(nstep-2)-5080.0 *f(nstep-3)+2970.0 *f(nstep-4)& 
   & -972.0 *f(nstep-5)+137.0 *f(nstep-6))/(180.0 *hr**2)
   return
end subroutine der_7p



function simpson(x,f) result (Inte)
! ---------------------------- !
! Performs Simpson integration !
! in x on the function f       !
! ---------------------------- !
   real(8)              :: inte, x(:),f(:)
   integer              :: i
   if ((mod(size(x),2) .ne. 0).or.(size(x) .ne. size(f))) stop "error in simpson integral!"

   Inte = 0
   do i=1, int(size(x)/2)
     Inte = Inte + 4. * pi * (x(2) - x(1))/3. * ( &
&      (x(2*i-1))**2    * f(2*i-1)  +&
&  4 * (x(2*i))  **2    * f(2*i)    +&
&      (x(2*i+1))**2    * f(2*i+1)  )
    enddo

end function simpson



 subroutine NLO_correction(np,r,u,v2) !result (v2)
! ---------------------------- !
! Still testing different      !
! approaches. Do not use this  !
! fucntion if you dont know    !
! what you are doing or expect !
! wrong results                !
! ---------------------------- !

  !  USE GAUSS_QUADRATURE
  integer                  :: np
  real(8)                  ::  u(:),   r(:)
  real(8),allocatable      :: du(:), ddu(:) 
  real(8)                  :: v2
  real(8)    :: alpha
  integer    :: n_point_quadrature, j, i
  real(8),allocatable      :: xg(:), wgs(:)
  integer,parameter        :: kind = 5
  
  
  allocate(du(np), ddu(np))
  call     der_7p(r,u,du,ddu)
  
  test=  LEC2(1) * (r(2)-r(1)) * sum(  exp(-.25 *( LEC(2)*r(:))**2)      * u(:) *   u(:) ) &
  &+LEC2(2) * (r(2)-r(1)) * sum(  exp(-.25 *( LEC(2)*r(:))**2)  * 2 * u(:) * ddu(:) )

  deallocate(du, ddu)
 end subroutine NLO_correction 


 subroutine Numerov_scattering_zeroE(Pot, Max_range, Binni, l, pot_DWBA, ERE) 
!----------------------------------------------------------------------------!
! Numerov for scattering lenght (ERE at 0 Energy)                            !
! The function returns (#of_boundstates, scattering_lenght, effective range) !
!----------------------------------------------------------------------------!

  integer                :: Binni, l
  real(8)                :: Max_range, E_step, Step, Numerov_precision, adwba
  real(8)                :: max_WF, R1, R2, A0, B0, ERE(4), Kappa, KR, R,  E
  real(8), allocatable   :: WF(:), erre(:)
  integer                :: nodi, i, Binnied
  real(8)                :: w0, w1, w2
  real(8)                :: a1
  

    interface
        function Pot(r) result(v)
            real(8), intent (in)  :: r
            real(8)               :: v
        end function Pot
    end interface

    interface
        function Pot_DWBA(r) result(v)
            real(8), intent (in)  :: r
            real(8)               :: v
        end function Pot_DWBA
    end interface

  allocate(WF(Binni), erre(Binni))

    E           = 0.                 ! Scattering at ZERO energy...
    Binnied     = Binni*0.9          ! Place where to calculate the asimptotic slope
    Step        = Max_range/Binni    ! Just the positional step
    max_WF      = WF(1)              ! For renormalization porpuse
    nodi        = 0                  ! Nodes of the function

    !From 0 to Binni (match)
    WF(1) = 0.
    WF(2) = Step**(l+1)
    w0    = 0
    w1    = (1-Step*Step/12._8*  f(Pot,E, Step, 0) ) *  WF(2)

    do i=3,Binni
          w2 = 2*w1 - w0 + Step*Step* f(Pot,E, Step, i-1) *WF(i-1)
          WF(i) = w2/(1-Step*Step/12._8* f(Pot,E, Step, i) )
          w0=w1
          w1=w2
      if (WF(i)*WF(i-1) < 0)         nodi   = nodi+1
      if (abs(WF(i))>max_WF)         max_WF = abs(WF(i))
    enddo

    !Calcolo phaseshift
    R1          = Distance(Step,Binni)
    R2          = Distance(Step,Binnied)
    B0          = (WF(Binnied)-WF(Binni)) / (R2-R1)
    A0          = WF(Binni) - B0*R1
    ERE(1)      = nodi+1
    ERE(2)      = -A0/B0

    !print result
   OPEN(9, FILE='Wave_function_sacattering.dat')
   write(9,*) "# r      psi(r)     V(r)     Assintotics_psi(r)"
   ERE(3) = 0.
   ERE(4) = 0.
   do i=1,Binni
    erre(i) = Step * (i-1)
    ERE(3)  = ERE(3) + Step * (2./(A0**2)) * (WF(i)**2 - (A0 + B0*Distance(Step,i))**2)
    ERE(4)  = ( step * WF(i) * Pot(Distance(Step,i))) 
    write(9,*) Distance(Step,i), WF(i), Pot(Distance(Step,i)),(A0 + B0*Distance(Step,i))
   enddo
   CLOSE(9)
   
   
   call NLO_correction(Binni, erre, WF, adwba)     ! DWBA correction (we are working on it)
   ERE(4) = adwba / ERE(4)                         ! renormaization
end subroutine Numerov_scattering_zeroE



subroutine Numerov_bound(Pot, energy, precision_, n_, l_, rmatch_, rmax_, r_, u_, v_)
!-------------------------------------------------------------------!
! Kernel of numerov algorithm. Just for boundstates (use scattering !
! function otherwise) but for all the spectrum.                     !
! Pot is the potential function pointer.                            !
! Energy is the initial energy and precision is the final precison  !
! n is the number of nodes, l the angular momentum                  !
! rmatch is the meeting point for the match in physical units       !
! rmax is the maximum box lenght (should be big enought)            !
! r is the distance array, u is the wave function and v will contain!
! the final potential.                                              !
!-------------------------------------------------------------------!
implicit none
real(8),parameter :: coef_(7)=(/49._8/20._8 ,- 6._8, 15._8/2._8, &
- 20._8/3._8, 15._8/4._8 ,- 6._8/5._8,  &
1._8/6._8 /) ! parameters for derivation
real(8), intent(inout) :: energy
real(8), allocatable   :: f(:)
real(8) :: w0,w1,w2,h,dlogu_out,dlogu_in, dlogu, old_dlogu
real(8) :: old_energy, new_energy,de,E_pert
integer :: s,i,iter, nodes
real(8) :: precision_, rmin_,rmax_,rmatch_
integer :: m_                       ! match point
integer :: n_                       ! number of grid points
integer :: l_                       ! angular momentum
real(8) :: r_(:),v_(:)              ! potential
real(8) :: u_(:)                    ! wavefunction


interface
function Pot(r) result(v)
real(8), intent (in)  :: r
real(8)               :: v
end function Pot
end interface
m_    =  ceiling( (rmatch_-rmin_)/(rmax_-rmin_)*n_)     ! matching point
de    =  1d-6                                           ! energy error
rmin_ =  0                                              ! Minimum r
h     = (rmax_-rmin_)/(n_-1)
allocate(f(n_))

! Creation of r(:) and v(:)
do i=1,n_
    r_(i) = h*(i-1)+rmin_
    v_(i) = Pot(r_(i))
enddo

! initial conditions
s = 3
u_(s-2) = 0
u_(s-1) = h**(l_+1)
iter=0
iteration: do while( (iter .lt. 2).or.(abs(de/energy) .gt. precision_) )
    nodes = 0
    ! set up of f function
    f =  (v_-energy)/hbaron2m
    w0 = u_(s-2)
    w1 = (1-h*h/12._8* f(s-1)) *  u_(s-1)
    do i=s,m_
        w2 = 2*w1 - w0 + h*h*f(i-1)*u_(i-1)
        u_(i) = w2/(1-h*h/12._8*f(i))
        w0=w1
        w1=w2
        if (u_(i)*u_(i-1) < 0) nodes = nodes +1
    end do

    ! now we can calculate the derivative at the matching point
    !        dlogu_out = dot_product(coef_, u(m_,m_-6:-1))
    dlogu_out = 49._8/20._8 * u_(m_) - 6._8 * u_(m_-1) + 15._8/2._8 * u_(m_-2)       &
    - 20._8/3._8 * u_(m_-3) + 15._8/4._8 * u_(m_-4) - 6._8/5._8*u_(m_-5) &
    + 1._8/6._8*u_(m_-6)
    dlogu_out = dlogu_out/h/u_(m_)


    ! inWard
    u_(n_) = exp(-r_(n_)*sqrt(abs(energy/hbaron2m)))      ! asymptotic bound state
    u_(n_-1) = u_(n_) * exp(h*sqrt(abs(energy/hbaron2m)))

    ! from inf to mathc
    w2 = (1-h*h/12._8* f(n_))   *  u_(n_)
    w1 = (1-h*h/12._8* f(n_-1)) *  u_(n_-1)
    do i=n_-2,m_,-1
        w0 = 2*w1 - w2 + h*h*f(i+1)*u_(i+1)
        u_(i) = w0/(1-h*h/12._8*f(i))
        w2 = w1
        w1 = w0
        if (u_(i)*u_(i+1) < 0) nodes = nodes +1
    end do

    ! now we can calculate the derivative at the matching point
    dlogu_in = -dot_product(coef_,u_(m_:m_+6))
    dlogu_in = dlogu_in/h/u_(m_)

    ! check for conditions
    dlogu = dlogu_out - dlogu_in


    if((nodes > 0).and.(energy<0)) then  !I want only groundstates no exitations now
      new_energy =  2*energy
      iter = -1
    else
      if(iter .eq. 0) then
          new_energy = energy + 1.d-5 * energy
      else
          new_energy = old_energy  - old_dlogu* (energy-old_energy)/(dlogu-old_dlogu)
      end if
    end if
    
    de = new_energy-energy
    old_energy = energy
    energy = new_energy
    old_dlogu =dlogu
    iter = iter + 1
end do iteration

! reconstruct the wf
f =  (v_-energy)/hbaron2m
w0 = u_(s-2)
w1 = (1-h*h/12._8* f(s-1)) *  u_(s-1)

! outWard
do i=s,n_
    w2 = 2*w1 - w0 + h*h*f(i-1)*u_(i-1)
    u_(i) = w2/(1-h*h/12._8*f(i))
    w0=w1
    w1=w2
end do

! time to renormalize the full wave function
u_(:)=u_(:)/r_(:)
u_(1)=u_(2)
u_(:)=u_(:)/sqrt(simpson(r_,u_**2)) 


deallocate(f)
end subroutine Numerov_bound


end module solvers








module Potential
use Constants
  implicit none
!-------------------------------------------------------------------!
! this module contains few useful potentials                        !
!-------------------------------------------------------------------!

   abstract interface
   function Type_Potential(r) result(v)
      real(8), intent(in)   :: r
      real(8)               :: v
   end function Type_Potential
   end interface
   procedure (Type_Potential), pointer ::  Pot_point_NLO => null (), Pot_point => null ()

   real(8), save      :: Pot_alpha
   contains

  function Pot_MTV(r) result(v)
      real(8) :: A   = 1458.05
      real(8) :: B   = 700.
      real(8) :: mua = 3.11
      real(8) :: mub = 1.55
      real(8), intent(in)  :: r
      real(8)              :: v
        v = A*exp(-mua*r)/r - B*exp(-mub*r)/r
   end function Pot_MTV

   function Pot_Volkov_inrho(rho,alpha) result(v)
      real(8) :: A   = 144.86
      real(8) :: B   = -83.34
      real(8) :: mua = 0.82
      real(8) :: mub = 1.60
      real(8), intent(in)  :: rho,alpha
      real(8)              :: v,r
        r = rho * cos(alpha)
        v = A*exp(-(r/mua)**2) + B*exp(-(r/mub)**2)
   end function Pot_Volkov_inrho

   function Pot_Volkov_rho(rho) result(v)
      real(8) :: A   = 144.86
      real(8) :: B   = -83.34
      real(8) :: mua = 0.82
      real(8) :: mub = 1.60
      real(8), intent(in)  :: rho
      real(8)              :: v,r
        r = rho * cos(Pot_alpha)
        v = A*exp(-(r/mua)**2) + B*exp(-(r/mub)**2)
   end function Pot_Volkov_rho

   function Buca_Q(r) result(v)
      real(8), intent(in)  :: r
      real(8)              :: v
        if(r<2.)  v=-50.!-23.3685!+25.*r**2
        if(r>=2.) v=0.
   end function Buca_Q

   function Pionless_LO(r) result(v)
      real(8), intent(in)  :: r
      real(8)              :: v
        v = LEC(1)*exp(-.25*(r*LEC(2))**2)
   end function Pionless_LO

   function Pionless_NLO(r) result(v)
      real(8), intent(in)  :: r
      real(8)              :: v
        v = .25 * LEC2(2) * ( LEC(2)**4 * r**2 - 6.*LEC(2)**2) * exp(-.25 *( LEC(2)*r)**2) &
                + LEC2(1) * exp(-.25 *( LEC(2)*r)**2)
   end function Pionless_NLO

   function Pionless_NLO_Batzalel(r) result(v)
      real(8), intent(in)  :: r
      real(8)              :: v
        v = .25 * LEC2(2) * ( LEC(2)**4 * r**2 - 6.*LEC(2)**2) * exp(-.25 *( LEC(2)*r)**2) &
                + LEC2(1) * exp(-.25 *( LEC(2)*r)**2)
   end function Pionless_NLO_Batzalel

   function No_pot(r) result(v)
      real(8), intent(in)  :: r
      real(8)              :: v
         v = 0.
   end function No_pot

end module Potential







program Numerov
  use solvers
  use Potential
  use Constants
  implicit none
  ! ------------------------------------------------------------------ !
  ! This program calculates the boundstate of a system with an         !
  ! interaction. Can be used both to fix lecs or to make predictions   ! 
  ! ------------------------------------------------------------------ !


  real(8)              :: C1lo, c2lo, C1nlo,C2nlo,C3nlo,C4nlo
  real(8)              :: Ei, Es, Inv, E_stima, Pot_Range, ERE(2), E_pert_c, E_pert_s, r
  integer              :: n_points_rho, i, caso_labda
  real                 :: LECmin, LECstep, BEnew, BEprime, BEtarg, energy, LECtemp, pionmass,R_square
  real(8), allocatable :: r_(:),v_(:)   ! potential
  real(8), allocatable :: u_(:)         ! wavefunction
  logical              :: def = .false.


  
! Get arguments from bash
integer :: icom, a
character(len=32) :: arg, strvec(4)
strvec = (/ "Energy   ", "Pion mass", "Λ       ", "LEC      " /)
do icom=1, 10
  call get_command_argument(icom, arg)
  if (len_trim(arg) == 0) then
   def = .true.
  endif 
 enddo
! ------------------------
write(6,*) "---Numerov:---"



!Standard potential choice (more choices in mod potential):
Pot_point      => Pionless_LO
Pot_point_NLO  => Pionless_NLO

   
if (def) then !Standard debug code:
    write(6,*) "Default system"
    E_stima               =  -11.5             ! Energy estimation
    n_points_rho          =  100000            ! Points rho
    mass         =  938.3                      ! Nucleon Mass
    hbaron2m     =  197.3269718**2 / (mass)    ! kinetic them
    Ei           =  E_stima                    ! Initial energy
    Es           =  0.0001                     ! Precision on the energy
    Inv          =  .1                         ! Matching point for numerov
    Pot_Range    =  4.                         ! Max range of the simulation
    LEC(2)       =  4.                         ! CUT-OFF
    LEC(1)       = -505.1643                   ! Main interaction
    LEC2(1)      =  0                          ! Perturbation
    LEC2(2)      =  0                          ! Second perturbation


else  ! If you are reading parameters from file.
    write(6,*) " System pionmass:", pionmass
    select case(int(pionmass)) ! What is the nucleon mass for each pion?
        case(140)
            mass = 938.0
        case(510)
            mass = 1320.0
        case(805)
            mass = 1634.0
        case default
            write(*,*)"Cannot use pion mass =", pionmass
            stop
      end select
    E_stima               =  -2.              ! Energy estimation
    n_points_rho          =  100000           ! Points rho
    hbaron2m              =  197.3269718**2 / (mass)       ! Kinetic therm
    Ei        =  E_stima                      ! Initial energy
    Es        =  0.001                        ! Precision on the energy
    Inv       =  0.5                          ! Matching point for numerov
    Pot_Range =  3.                           ! Max range of the simulation
    LEC(2)  =  4.                             ! CUT-OFF
    LEC(1)  = C1lo  + C2lo                    ! Channel coupling
    LEC2(1) = C1nlo + C2nlo
    LEC2(2) = C3nlo + C4nlo

endif





    write (6,*)   " Mass    = ", mass
    write (6,*)   " Cut-off = ", LEC(2)
    write (6,*)   " LEC     = ", LEC(1)
    write (6,*)   " Hbar/2m = ", hbaron2m

    allocate(r_(n_points_rho),u_(n_points_rho),v_(n_points_rho))
    ! kernel !
    call Numerov_bound(Pot_point, Ei, Es, n_points_rho, 0, Inv,Pot_Range, r_(:), u_(:), v_(:))

    write(6,*) " "
    write(6,*) " "
    write(6,*) " LO  Energy:", Ei
    write(6,*) " ---------------- "





deallocate(u_,v_)
end program Numerov



! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2023, Miguel A. Caro
! HND X
! HND X   soap_turbo is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   soap_turbo is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


module soap_turbo_radial_op

  contains

!**************************************************************************
!
! This function returns the normalization coefficients for the polynomial
! basis functions used to construct the orthonormal basis.
!
  function N_a(rcut, a)
    implicit none

    integer :: a, b
    real*8 :: rcut, N_a

    b = 2*a + 5

!   **************** New basis ******************
!    N_a = dsqrt( rcut**b / dfloat(b) )
    N_a = dsqrt( rcut / dfloat(b) )
!   *********************************************

  return
  end function
!**************************************************************************



!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and a polynomial piecewise representation for 
! the atomic density
!
  subroutine get_radial_expansion_coefficients_poly3operator_fast(n_sites, n_neigh, rjs_in, alpha_max, &
                                                             rcut_soft_in, rcut_hard_in, atom_sigma_in, &
                                                             atom_sigma_scaling, amplitude_scaling, W, &
                                                             scaling_mode, mask, radial_enhancement, &
                                                             do_derivatives, do_central, central_weight, &
                                                             exp_coeff, exp_coeff_der)

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling, central_weight
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
!
    integer :: i, j, k
    real*8 :: pi, rj, s2, atom_width, atom_width_scaling, filter_width, x
    real*8 :: lim_soft(1:3), lim_buffer(1:3)!, B(1:7) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real*8 :: M_left(1:7, 1:2), M_right(1:7, 1:2)
    real*8, allocatable :: A(:,:), I0(:,:), I_left(:), I_right(:)
    real*8 :: W(:,:)
!   derivatives
    real*8 :: amplitude_der
    real*8 :: poly_left(1:4), poly_right(1:4), poly_filter(1:4), B_der(1:7)
    real*8 :: M_left_der(1:7, 1:2), M_right_der(1:7, 1:2)
    real*8, allocatable :: I_left_der(:), I_right_der(:), upper_limit_der(:), lower_limit_der(:)
!   Results will be stored in exp_coeff, which is an array of dimension (alpha_max, n_atom_pairs)
    real*8 :: exp_coeff(:,:), exp_coeff_der(:,:)
    real*8, allocatable :: exp_coeff_soft(:), exp_coeff_buffer(:)
    real*8, allocatable :: exp_coeff_soft_der(:), exp_coeff_buffer_der(:) 
    logical, save :: print_basis = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, allocatable :: rjs(:), lim_soft_array(:,:), atom_sigma_scaleds(:), s2s(:), atom_widths(:), I0_array(:,:,:), &
                       g_aux_left_array(:,:,:), g_aux_right_array(:,:,:), M_left_array(:,:,:), M_right_array(:,:,:), &
                       exp_coeff_soft_array(:,:), amplitudes(:), amplitudes_der(:), I_left_array(:,:), &
                       I_right_array(:,:), B(:,:), exp_coeff_buffer_array(:,:), lim_buffer_array(:,:)
real*8 :: vect(1:4)
integer, allocatable :: rjs_idx(:), soft_weights(:), buffer_weights(:)
integer :: nn, k2
real*8, save :: elapsed_time = 0.d0
real*8 :: t1, t2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate( exp_coeff_soft(1:alpha_max) )
    allocate( exp_coeff_buffer(1:alpha_max) )
    allocate( A(1:alpha_max, 1:7) )
    allocate( I0(1:alpha_max + 4, 1:3) )
    allocate( I_left(1:alpha_max) )
    allocate( I_right(1:alpha_max) )

    allocate( exp_coeff_soft_der(1:alpha_max) )
    allocate( exp_coeff_buffer_der(1:alpha_max) )
    allocate( I_left_der(1:alpha_max) )
    allocate( I_right_der(1:alpha_max) )
    allocate( lower_limit_der(1:alpha_max) )
    allocate( upper_limit_der(1:alpha_max) )


!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      do j=1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,*) "0."
      end do
    end if

!   *************** New basis *******************
!   Normalized with rcut_hard to avoid instabilities
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
!   *********************************************
    pi = dacos(-1.d0)
    filter_width = 2.d0*sqrt(2.d0*log(2.d0))*(rcut_hard - rcut_soft)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if

    call get_constant_poly_coeff(alpha_max, rcut_hard, A)


! NOTE: the vectorization can potentially be done for all atoms at one, where nn = size(rjs_in), then only
! obvious thing to be careful about I can think of is handling of the central atoms in the flattened array
    k = 0
    do i = 1, n_sites


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! check number of atoms within soap cutoff
!nn = count( rjs_in(k+1:k+n_neigh(i)) < rcut_hard_in .and. mask(k+1:k+n_neigh(i)) )
!if( nn < 1 )then
!  cycle
!end if





! SOFT REGION
! Temporarily allocate this as oversized array:
allocate( atom_widths(1:n_neigh(i)) )
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*(atom_sigma_in + atom_sigma_scaling*rjs_in(k+1:k+n_neigh(i)))

! count number of atoms fully within soft region
nn = count( rjs_in(k+1:k+n_neigh(i)) - atom_widths < rcut_soft_in .and. mask(k+1:k+n_neigh(i)) )

if( nn > 0 )then

! These need to be allocated immediately
allocate( rjs(1:nn) )
allocate( rjs_idx(1:nn) )

! make assignment between original neighbors and reduced neighbors
k2 = 1
do j = 1, n_neigh(i)
!  if( rjs_in(k+j) < rcut_hard_in .and. mask(k+j) )then
  if( rjs_in(k+j) - atom_widths(j) < rcut_soft_in .and. mask(k+j) )then
    rjs(k2) = rjs_in(k+j)
    rjs_idx(k2) = k+j
    k2 = k2 + 1
  end if
end do

deallocate( atom_widths )




! allocate arrays that depend on the number of neighbors inside soap cutoff
allocate( atom_sigma_scaleds(1:nn) )
allocate( s2s(1:nn) )
allocate( atom_widths(1:nn) )
allocate( lim_soft_array(1:nn, 1:3) )
allocate( soft_weights(1:nn) )
allocate( buffer_weights(1:nn) )
! we might be able to get some extra speed by swapping the order of some of these dimensions
allocate( I0_array(1:nn, 1:alpha_max, 1:3) )
allocate( g_aux_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( g_aux_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( M_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( M_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( I_left_array(1:nn, 1:alpha_max) )
allocate( I_right_array(1:nn, 1:alpha_max) )
allocate( exp_coeff_soft_array(1:nn, 1:alpha_max) )
allocate( exp_coeff_buffer_array(1:nn, 1:alpha_max) )
allocate( amplitudes(1:nn) )
allocate( amplitudes_der(1:nn) )
!allocate( B(1:nn, 1:7) )
allocate( B(1:7, 1:nn) )



! define atom parameters in vector form
!   **************** New basis ******************
rjs = rjs/rcut_hard_in
!   *********************************************
atom_sigma_scaleds = atom_sigma + atom_sigma_scaling*rjs
s2s = atom_sigma_scaleds**2
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaleds
atom_width_scaling = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaling
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitudes = 1.d0 / atom_sigma_scaleds
              amplitudes_der = - atom_sigma_scaling / s2s 
!            else if( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 <= 1.d-10 )then ! FIX THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              amplitudes = 0.d0
!              amplitude_ders = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )
                amplitudes_der = 6.d0 / atom_sigma_scaleds * (rjs**2 - rjs) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              else
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**amplitude_scaling
                amplitudes_der = 6.d0*amplitude_scaling / atom_sigma_scaleds * (rjs**2 - rjs) &
                                * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( size(amplitudes) > 0 )then
            if( rjs_idx(1) == k+1 )then
              amplitudes(1) = central_weight * amplitudes(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              amplitudes_der(1) = central_weight * amplitudes_der(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
          end if

!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitudes_der = amplitudes * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                             amplitudes_der * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
            amplitudes = amplitudes * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
          else if( radial_enhancement == 2 )then
            amplitudes_der = amplitudes*( 2.d0*rjs + 2.d0*atom_sigma_scaleds*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaleds + dsqrt(8.d0/pi)*rjs*atom_sigma_scaling ) + &
                            amplitudes_der*( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
            amplitudes = amplitudes * ( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
          end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! NOTE: this below is design to be used with arrays that contain all the atoms. In some cases, the left_weights or
! right_weights might contain lots of zeros (*), in which case a speed up of close to x2 could be achieved by not doing
! all atoms but only those for which the weights are 1. One could refactor this code with the nn construction to
! have, e.g., nn_soft and nn_buffer
! (*) we are allocating and computing contributions from soft and buffer zones for all atoms, even when they do not
!     have any


! define arrays that give contribution of atoms within left and right of rcut_soft
! note that these are integer arrays (0|1 = .false.|.true.)
!soft_weights = rjs - atom_widths < rcut_soft
!nn_soft = count(soft_weights)
!buffer_weights = rcut_soft < rcut_hard .and. rjs + atom_widths > rcut_soft
!nn_buffer = count(buffer_weights)




! handle the central atom contribution; NOTE: I'M NOT SURE ABOUT THE NEED FOR THIS DO_CENTRAL TAG HERE, IT SHOULD BE HANDLED VIA THE CENTRAL_WEIGHT
if( .not. do_central )then
  if( size(amplitudes) > 0 )then
    if( rjs_idx(1) == k+1 )then
      amplitudes(1) = 0.d0
    end if
  end if
end if



! integration limits inside r_soft
lim_soft_array(:, 1) = max( 0.d0, rjs - atom_widths ) ! lower limit left
!lim_soft_array(:, 2) = min( max(0.d0, rjs), rcut_soft )   ! upper limit left / lower limit right ! rj should always be > 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lim_soft_array(:, 2) = min( rjs, rcut_soft )               ! upper limit left / lower limit right
lim_soft_array(:, 3) = min( rcut_soft, rjs + atom_widths ) ! upper limit right



! dimensions of M_radial_poly_array are (1:nn, 1:alpha_max, 1:3)
I0_array = M_radial_poly_array(lim_soft_array, alpha_max + 4, rcut_hard)
g_aux_left_array(1:nn, 1:4, 1:2) = g_aux_array(lim_soft_array(:, 1:2), rjs(:), atom_widths(:), "left")
g_aux_right_array(1:nn, 1:4, 2:3) = g_aux_array(lim_soft_array(:, 2:3), rjs(:), atom_widths(:), "right")


vect = [-1.d0, -1.d0, -2.d0, -6.d0]
do j = 1, 4
  M_left_array(1:nn, j, 1) = I0_array(1:nn, j, 1) * g_aux_left_array(1:nn, j, 1) * vect(j)
end do
do j = 1, 4
  M_left_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_left_array(1:nn, j, 2) * vect(j)
end do
do j = 1, 4
  M_right_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_right_array(1:nn, j, 2) * vect(j)
end do
do j = 1, 4
  M_right_array(1:nn, j, 3) = I0_array(1:nn, j, 3) * g_aux_right_array(1:nn, j, 3) * vect(j)
end do

I_left_array(1:nn, 1:alpha_max) = matmul( M_left_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                  matmul( M_left_array(1:nn, 1:4, 1), transpose(A(1:alpha_max, 1:4)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 1)
I_right_array(1:nn, 1:alpha_max) = matmul( M_right_array(1:nn, 1:4, 3), transpose(A(1:alpha_max, 1:4)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                   matmul( M_right_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 2)

exp_coeff_soft_array = I_left_array + I_right_array

do j = 1, nn
  k2 = rjs_idx(j)
  exp_coeff(1:alpha_max, k2) = amplitudes(j) * exp_coeff_soft_array(j, 1:alpha_max)
end do


deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_soft_array, soft_weights, &
            buffer_weights, I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
            I_left_array, I_right_array, exp_coeff_soft_array, amplitudes, amplitudes_der, B, &
            exp_coeff_buffer_array )

else

deallocate( atom_widths )

end if













! BUFFER REGION
! Temporarily allocate this as oversized array:
allocate( atom_widths(1:n_neigh(i)) )
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*(atom_sigma_in + atom_sigma_scaling*rjs_in(k+1:k+n_neigh(i)))

! count number of atoms within buffer region
nn = count( rcut_soft_in < rcut_hard_in .and. rjs_in(k+1:k+n_neigh(i)) + atom_widths > rcut_soft_in .and. mask(k+1:k+n_neigh(i)) )

if( nn > 0 )then

! These need to be allocated immediately
allocate( rjs(1:nn) )
allocate( rjs_idx(1:nn) )

! make assignment between original neighbors and reduced neighbors
k2 = 1
do j = 1, n_neigh(i)
!  if( rjs_in(k+j) < rcut_hard_in .and. mask(k+j) )then
  if( rcut_soft_in < rcut_hard_in .and. rjs_in(k+j) + atom_widths(j) > rcut_soft_in .and. mask(k+j) )then
    rjs(k2) = rjs_in(k+j)
    rjs_idx(k2) = k+j
    k2 = k2 + 1
  end if
end do

deallocate( atom_widths )




! allocate arrays that depend on the number of neighbors inside soap cutoff
allocate( atom_sigma_scaleds(1:nn) )
allocate( s2s(1:nn) )
allocate( atom_widths(1:nn) )
allocate( lim_buffer_array(1:nn, 1:3) )
allocate( soft_weights(1:nn) )
allocate( buffer_weights(1:nn) )
! we might be able to get some extra speed by swapping the order of some of these dimensions
allocate( I0_array(1:nn, 1:alpha_max, 1:3) )
allocate( g_aux_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( g_aux_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( M_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( M_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( I_left_array(1:nn, 1:alpha_max) )
allocate( I_right_array(1:nn, 1:alpha_max) )
allocate( exp_coeff_soft_array(1:nn, 1:alpha_max) )
allocate( exp_coeff_buffer_array(1:nn, 1:alpha_max) )
allocate( amplitudes(1:nn) )
allocate( amplitudes_der(1:nn) )
!allocate( B(1:nn, 1:7) )
allocate( B(1:7, 1:nn) )



! define atom parameters in vector form
!   **************** New basis ******************
rjs = rjs/rcut_hard_in
!   *********************************************
atom_sigma_scaleds = atom_sigma + atom_sigma_scaling*rjs
s2s = atom_sigma_scaleds**2
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaleds
atom_width_scaling = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaling
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitudes = 1.d0 / atom_sigma_scaleds
              amplitudes_der = - atom_sigma_scaling / s2s 
!            else if( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 <= 1.d-10 )then ! FIX THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              amplitudes = 0.d0
!              amplitude_ders = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )
                amplitudes_der = 6.d0 / atom_sigma_scaleds * (rjs**2 - rjs) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              else
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**amplitude_scaling
                amplitudes_der = 6.d0*amplitude_scaling / atom_sigma_scaleds * (rjs**2 - rjs) &
                                * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( size(amplitudes) > 0 )then
            if( rjs_idx(1) == k+1 )then
              amplitudes(1) = central_weight * amplitudes(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              amplitudes_der(1) = central_weight * amplitudes_der(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
          end if

!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitudes_der = amplitudes * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                             amplitudes_der * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
            amplitudes = amplitudes * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
          else if( radial_enhancement == 2 )then
            amplitudes_der = amplitudes*( 2.d0*rjs + 2.d0*atom_sigma_scaleds*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaleds + dsqrt(8.d0/pi)*rjs*atom_sigma_scaling ) + &
                            amplitudes_der*( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
            amplitudes = amplitudes * ( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
          end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! NOTE: this below is design to be used with arrays that contain all the atoms. In some cases, the left_weights or
! right_weights might contain lots of zeros, in which case a speed up of close to x2 could be achieved by not doing
! all atoms but only those for which the weights are 1. One could refactor this code with the nn construction to
! have, e.g., nn_soft and nn_buffer


! define arrays that give contribution of atoms within left and right of rcut_soft
! note that these are integer arrays (0|1 = .false.|.true.)
!soft_weights = rjs - atom_widths < rcut_soft
!nn_soft = count(soft_weights)
!buffer_weights = rcut_soft < rcut_hard .and. rjs + atom_widths > rcut_soft
!nn_buffer = count(buffer_weights)




! handle the central atom contribution; NOTE: I'M NOT SURE ABOUT THE NEED FOR THIS DO_CENTRAL TAG HERE, IT SHOULD BE HANDLED VIA THE CENTRAL_WEIGHT
if( .not. do_central )then
  if( size(amplitudes) > 0 )then
    if( rjs_idx(1) == k+1 )then
      amplitudes(1) = 0.d0
    end if
  end if
end if


! Atoms within the buffer region

lim_buffer_array(:, 1) = max( rcut_soft, rjs - atom_widths ) ! lower limit left
lim_buffer_array(:, 2) = max( rjs, rcut_soft ) ! upper limit left / lower limit right
lim_buffer_array(:, 3) = min( rcut_hard, rjs + atom_widths ) ! upper limit right
I0_array = M_radial_poly_array(lim_buffer_array, max(7, alpha_max + 4), rcut_hard)

call get_constant_poly_filter_coeff_array(rjs, atom_widths, rcut_soft, filter_width, 'left', B)

! We should try to figure out a more "vectorized way" of doing this
do k2 = 1, nn
M_left_array(k2, 1:7, 1) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 1), 6) ) * I0_array(k2, 1:7, 1)
M_left_array(k2, 1:7, 2) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * I0_array(k2, 1:7, 2)
end do

I_left_array(1:nn, 1:alpha_max) = matmul( M_left_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                  matmul( M_left_array(1:nn, 1:7, 1), transpose(A(1:alpha_max, 1:7)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 1)

call get_constant_poly_filter_coeff_array(rjs, atom_widths, rcut_soft, filter_width, 'right', B)

! We should try to figure out a more "vectorized way" of doing this
do k2 = 1, nn
M_right_array(k2, 1:7, 2) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * I0_array(k2, 1:7, 2)
M_right_array(k2, 1:7, 3) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 3), 6) ) * I0_array(k2, 1:7, 3)
end do

I_right_array(1:nn, 1:alpha_max) = matmul( M_right_array(1:nn, 1:7, 3), transpose(A(1:alpha_max, 1:7)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                   matmul( M_right_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 2)

exp_coeff_buffer_array = I_left_array + I_right_array


do j = 1, nn
  k2 = rjs_idx(j)
  exp_coeff(1:alpha_max, k2) = exp_coeff(1:alpha_max, k2) + amplitudes(j) * exp_coeff_buffer_array(j, 1:alpha_max)
end do


deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_buffer_array, soft_weights, &
            buffer_weights, I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
            I_left_array, I_right_array, exp_coeff_soft_array, amplitudes, amplitudes_der, B, &
            exp_coeff_buffer_array )

else

deallocate( atom_widths )

end if








! Transform from g_alpha to g_n (the orthonormal basis)
!exp_coeff_soft_array(1:nn, 1:alpha_max) = matmul( exp_coeff_soft_array(1:nn, 1:alpha_max), W )

!exp_coeff_buffer_array(1:nn, 1:alpha_max) = matmul( exp_coeff_buffer_array(1:nn, 1:alpha_max), W )



!do j = 1, nn
!  k2 = rjs_idx(j)
!  exp_coeff(1:alpha_max, k2) = amplitudes(j) * soft_weights(j) * exp_coeff_soft_array(j, 1:alpha_max) + &
!                               amplitudes(j) * buffer_weights(j) * exp_coeff_buffer_array(j, 1:alpha_max)
!end do


! Transform from g_alpha to g_n (the orthonormal basis)
exp_coeff(1:alpha_max, k+1:k+n_neigh(i)) = matmul( W, exp_coeff(1:alpha_max, k+1:k+n_neigh(i)) )


k = k + n_neigh(i)

!deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_soft_array, soft_weights, &
!            buffer_weights, I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
!            I_left_array, I_right_array, exp_coeff_soft_array, amplitudes, amplitudes_der, B, &
!            exp_coeff_buffer_array )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do

!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
    if( do_derivatives )then
      exp_coeff_der = exp_coeff_der / dsqrt(rcut_hard_in)
    end if
!   *********************************************

!   This is for debugging
    if( .true. )then
      open(10, file="radial_expansion_coefficients_operator-fast.dat", status="unknown")
      do i=1, size(exp_coeff, 2)
        write(10,*) exp_coeff(1:alpha_max, i)
      end do
      close(10)
      if( do_derivatives )then
        open(10, file="radial_expansion_coefficients_der_operator-fast.dat", status="unknown")
        do i=1, size(exp_coeff, 2)
          write(10,*) exp_coeff_der(1:alpha_max, i)
        end do
        close(10)
      end if
    end if

    deallocate( exp_coeff_soft, exp_coeff_buffer, A, I0, I_left, I_right )
    deallocate( exp_coeff_soft_der, exp_coeff_buffer_der )

  return
  end subroutine
!**************************************************************************





!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and a polynomial piecewise representation for 
! the atomic density
!
  subroutine get_radial_expansion_coefficients_poly3operator_fast_der(n_sites, n_neigh, rjs_in, alpha_max, &
                                                             rcut_soft_in, rcut_hard_in, atom_sigma_in, &
                                                             atom_sigma_scaling, amplitude_scaling, W, &
                                                             scaling_mode, mask, radial_enhancement, &
                                                             do_derivatives, do_central, central_weight, &
                                                             exp_coeff, exp_coeff_der)

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling, central_weight
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
!
    integer :: i, j, k
    real*8 :: pi, rj, s2, atom_width, atom_width_scaling, filter_width, x
    real*8, allocatable :: A(:,:)
    real*8 :: W(:,:)
!   derivatives
    real*8 :: amplitude_der
!   Results will be stored in exp_coeff, which is an array of dimension (alpha_max, n_atom_pairs)
    real*8 :: exp_coeff(:,:), exp_coeff_der(:,:)
    logical, save :: print_basis = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8, allocatable :: rjs(:), lim_soft_array(:,:), atom_sigma_scaleds(:), s2s(:), atom_widths(:), I0_array(:,:,:), &
                       g_aux_left_array(:,:,:), g_aux_right_array(:,:,:), M_left_array(:,:,:), M_right_array(:,:,:), &
                       exp_coeff_soft_array(:,:), amplitudes(:), I_left_array(:,:), &
                       I_right_array(:,:), B(:,:), exp_coeff_buffer_array(:,:), lim_buffer_array(:,:)
!!!! derivatives
real*8, allocatable :: g_aux_left_der_array(:,:,:), g_aux_right_der_array(:,:,:), M_left_der_array(:,:,:), &
                       M_right_der_array(:,:,:), exp_coeff_soft_der_array(:,:), amplitudes_der(:), &
                       poly_left(:,:), poly_right(:,:), I0_block(:,:), B_der(:,:), &
                       I_left_der_array(:,:), I_right_der_array(:,:), exp_coeff_buffer_der_array(:,:)
logical, allocatable :: mask_der_array(:)
integer, allocatable :: upper_limit_der_array(:,:), lower_limit_der_array(:,:), I_ones(:,:)
integer :: nn_lower_lim_der, nn_upper_lim_der

real*8 :: vect(1:4)
integer, allocatable :: rjs_idx(:), soft_weights(:), buffer_weights(:)
integer :: nn, k2
real*8, save :: elapsed_time = 0.d0
real*8 :: t1, t2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      do j=1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,*) "0."
      end do
    end if

!   *************** New basis *******************
!   Normalized with rcut_hard to avoid instabilities
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
!   *********************************************
    pi = dacos(-1.d0)
    filter_width = 2.d0*sqrt(2.d0*log(2.d0))*(rcut_hard - rcut_soft)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if

    allocate( A(1:alpha_max, 1:7) )
    call get_constant_poly_coeff(alpha_max, rcut_hard, A)

    allocate( I_ones(1:1, 1:alpha_max) )
    I_ones = 1

! NOTE: the vectorization can potentially be done for all atoms at one, where nn = size(rjs_in), then only
! obvious thing to be careful about I can think of is handling of the central atoms in the flattened array
    k = 0
    do i = 1, n_sites


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! check number of atoms within soap cutoff
!nn = count( rjs_in(k+1:k+n_neigh(i)) < rcut_hard_in .and. mask(k+1:k+n_neigh(i)) )
!if( nn < 1 )then
!  cycle
!end if





! SOFT REGION
! Temporarily allocate this as oversized array:
allocate( atom_widths(1:n_neigh(i)) )
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*(atom_sigma_in + atom_sigma_scaling*rjs_in(k+1:k+n_neigh(i)))

! count number of atoms fully within soft region
nn = count( rjs_in(k+1:k+n_neigh(i)) - atom_widths < rcut_soft_in .and. mask(k+1:k+n_neigh(i)) )

if( nn > 0 )then

! These need to be allocated immediately
allocate( rjs(1:nn) )
allocate( rjs_idx(1:nn) )

! make assignment between original neighbors and reduced neighbors
k2 = 1
do j = 1, n_neigh(i)
!  if( rjs_in(k+j) < rcut_hard_in .and. mask(k+j) )then
  if( rjs_in(k+j) - atom_widths(j) < rcut_soft_in .and. mask(k+j) )then
    rjs(k2) = rjs_in(k+j)
    rjs_idx(k2) = k+j
    k2 = k2 + 1
  end if
end do

deallocate( atom_widths )




! allocate arrays that depend on the number of neighbors inside soap cutoff
allocate( atom_sigma_scaleds(1:nn) )
allocate( s2s(1:nn) )
allocate( atom_widths(1:nn) )
allocate( lim_soft_array(1:nn, 1:3) )
allocate( soft_weights(1:nn) )
!allocate( buffer_weights(1:nn) )
! we might be able to get some extra speed by swapping the order of some of these dimensions
allocate( I0_array(1:nn, 1:alpha_max, 1:3) )
allocate( g_aux_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( g_aux_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( M_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( M_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( I_left_array(1:nn, 1:alpha_max) )
allocate( I_right_array(1:nn, 1:alpha_max) )
allocate( exp_coeff_soft_array(1:nn, 1:alpha_max) )
!allocate( exp_coeff_buffer_array(1:nn, 1:alpha_max) )
allocate( amplitudes(1:nn) )
allocate( amplitudes_der(1:nn) )
!allocate( B(1:nn, 1:7) )
allocate( B(1:7, 1:nn) )

if( do_derivatives )then
allocate( mask_der_array(1:nn) )
allocate( lower_limit_der_array(1:nn, 1:1) )
allocate( upper_limit_der_array(1:nn, 1:1) )
allocate( g_aux_left_der_array(1:nn, 1:alpha_max, 1:2) )
allocate( g_aux_right_der_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( M_left_der_array(1:nn, 1:alpha_max, 1:2) )
allocate( M_right_der_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( I_left_der_array(1:nn, 1:alpha_max) )
allocate( I_right_der_array(1:nn, 1:alpha_max) )
allocate( poly_left(1:nn, 1:2) )
allocate( poly_right(1:nn, 1:2) )
allocate( I0_block(1:nn, 1:alpha_max) )
allocate( exp_coeff_soft_der_array(1:nn, 1:alpha_max) )
end if

! define atom parameters in vector form
!   **************** New basis ******************
rjs = rjs/rcut_hard_in
!   *********************************************
atom_sigma_scaleds = atom_sigma + atom_sigma_scaling*rjs
s2s = atom_sigma_scaleds**2
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaleds
atom_width_scaling = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaling
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitudes = 1.d0 / atom_sigma_scaleds
              amplitudes_der = - atom_sigma_scaling / s2s 
!            else if( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 <= 1.d-10 )then ! FIX THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              amplitudes = 0.d0
!              amplitude_ders = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )
                amplitudes_der = 6.d0 / atom_sigma_scaleds * (rjs**2 - rjs) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              else
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**amplitude_scaling
                amplitudes_der = 6.d0*amplitude_scaling / atom_sigma_scaleds * (rjs**2 - rjs) &
                                * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( size(amplitudes) > 0 )then
            if( rjs_idx(1) == k+1 )then
              amplitudes(1) = central_weight * amplitudes(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              amplitudes_der(1) = central_weight * amplitudes_der(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
          end if

!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitudes_der = amplitudes * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                             amplitudes_der * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
            amplitudes = amplitudes * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
          else if( radial_enhancement == 2 )then
            amplitudes_der = amplitudes*( 2.d0*rjs + 2.d0*atom_sigma_scaleds*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaleds + dsqrt(8.d0/pi)*rjs*atom_sigma_scaling ) + &
                            amplitudes_der*( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
            amplitudes = amplitudes * ( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
          end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! NOTE: this below is design to be used with arrays that contain all the atoms. In some cases, the left_weights or
! right_weights might contain lots of zeros, in which case a speed up of close to x2 could be achieved by not doing
! all atoms but only those for which the weights are 1. One could refactor this code with the nn construction to
! have, e.g., nn_soft and nn_buffer


! define arrays that give contribution of atoms within left and right of rcut_soft
! note that these are integer arrays (0|1 = .false.|.true.)
!soft_weights = rjs - atom_widths < rcut_soft
!nn_soft = count(soft_weights)
!buffer_weights = rcut_soft < rcut_hard .and. rjs + atom_widths > rcut_soft
!nn_buffer = count(buffer_weights)




! handle the central atom contribution; NOTE: I'M NOT SURE ABOUT THE NEED FOR THIS DO_CENTRAL TAG HERE, IT SHOULD BE HANDLED VIA THE CENTRAL_WEIGHT
if( .not. do_central )then
  if( size(amplitudes) > 0 )then
    if( rjs_idx(1) == k+1 )then
      amplitudes(1) = 0.d0
    end if
  end if
end if



! integration limits inside r_soft
lim_soft_array(:, 1) = max( 0.d0, rjs - atom_widths ) ! lower limit left
!lim_soft_array(:, 2) = min( max(0.d0, rjs), rcut_soft )   ! upper limit left / lower limit right ! rj should always be > 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
lim_soft_array(:, 2) = min( rjs, rcut_soft )               ! upper limit left / lower limit right
lim_soft_array(:, 3) = min( rcut_soft, rjs + atom_widths ) ! upper limit right



! dimensions of M_radial_poly_array are (1:nn, 1:alpha_max, 1:3)
I0_array = M_radial_poly_array(lim_soft_array, alpha_max + 4, rcut_hard)
g_aux_left_array(1:nn, 1:4, 1:2) = g_aux_array(lim_soft_array(:, 1:2), rjs(:), atom_widths(:), "left")
g_aux_right_array(1:nn, 1:4, 2:3) = g_aux_array(lim_soft_array(:, 2:3), rjs(:), atom_widths(:), "right")


vect = [-1.d0, -1.d0, -2.d0, -6.d0]
do j = 1, 4
  M_left_array(1:nn, j, 1) = I0_array(1:nn, j, 1) * g_aux_left_array(1:nn, j, 1) * vect(j)
end do
do j = 1, 4
  M_left_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_left_array(1:nn, j, 2) * vect(j)
end do
do j = 1, 4
  M_right_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_right_array(1:nn, j, 2) * vect(j)
end do
do j = 1, 4
  M_right_array(1:nn, j, 3) = I0_array(1:nn, j, 3) * g_aux_right_array(1:nn, j, 3) * vect(j)
end do

I_left_array(1:nn, 1:alpha_max) = matmul( M_left_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                  matmul( M_left_array(1:nn, 1:4, 1), transpose(A(1:alpha_max, 1:4)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 1)
I_right_array(1:nn, 1:alpha_max) = matmul( M_right_array(1:nn, 1:4, 3), transpose(A(1:alpha_max, 1:4)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                   matmul( M_right_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 2)

exp_coeff_soft_array = I_left_array + I_right_array

do j = 1, nn
  k2 = rjs_idx(j)
  exp_coeff(1:alpha_max, k2) = amplitudes(j) * exp_coeff_soft_array(j, 1:alpha_max)
end do


if( do_derivatives )then
! contribution to derivatives --> (3) in notes
  g_aux_left_der_array(1:nn, 1:4, 1:2) = g_aux_der_array(lim_soft_array(:, 1:2), rjs(:), atom_widths(:), &
                                                         atom_width_scaling, "left")
  g_aux_right_der_array(1:nn, 1:4, 2:3) = g_aux_der_array(lim_soft_array(:, 2:3), rjs(:), atom_widths(:), &
                                                          atom_width_scaling, "right")

  do j = 1, 4
    M_left_der_array(1:nn, j, 1) = I0_array(1:nn, j, 1) * g_aux_left_der_array(1:nn, j, 1) * vect(j)
  end do
  do j = 1, 4
    M_left_der_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_left_der_array(1:nn, j, 2) * vect(j)
  end do
  do j = 1, 4
    M_right_der_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_right_der_array(1:nn, j, 2) * vect(j)
  end do
  do j = 1, 4
    M_right_der_array(1:nn, j, 3) = I0_array(1:nn, j, 3) * g_aux_right_der_array(1:nn, j, 3) * vect(j)
  end do

  I_left_der_array(1:nn, 1:alpha_max) = matmul( M_left_der_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                        I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                        matmul( M_left_der_array(1:nn, 1:4, 1), transpose(A(1:alpha_max, 1:4)) ) * &
                                      I0_array(1:nn, 5:alpha_max + 4, 1)
  I_right_der_array(1:nn, 1:alpha_max) = matmul( M_right_der_array(1:nn, 1:4, 3), transpose(A(1:alpha_max, 1:4)) ) * &
                                         I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                         matmul( M_right_der_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                         I0_array(1:nn, 5:alpha_max + 4, 2)
  exp_coeff_soft_der_array = I_left_der_array + I_right_der_array
! add derivatives of the integral limits if necessary --> (1), (2) in notes
! check dimensions!!!!!!!!!!!!!!!!
  mask_der_array = (lim_soft_array(:, 1) == rjs - atom_widths)
  nn_lower_lim_der = count(mask_der_array)
  if( nn_lower_lim_der > 0 )then
    lower_limit_der_array(1:nn, 1) = merge(1, 0, mask_der_array)
    poly_left = g_aux_left_der_array(1:nn, 1, 1:1)
    I0_block = I0_array(1:nn, 4:alpha_max + 4, 1)
    exp_coeff_soft_der_array = exp_coeff_soft_der_array - (1.d0 - atom_width_scaling) * &
                               matmul( lower_limit_der_array * poly_left, I_ones ) * I0_block
  end if 
  mask_der_array = (lim_soft_array(:, 3) == rjs + atom_widths)
  nn_upper_lim_der = count(mask_der_array)
  if( nn_upper_lim_der > 0 )then
    upper_limit_der_array(1:nn, 1) = merge(1, 0, mask_der_array)
    poly_right = g_aux_right_der_array(1:nn, 1, 3:3)
    I0_block = I0_array(1:nn, 4:alpha_max + 4, 3)
    exp_coeff_soft_der_array = exp_coeff_soft_der_array + (1.d0 + atom_width_scaling) * &
                               matmul( upper_limit_der_array * poly_right, I_ones ) * I0_block
  end if  
  do j = 1, nn
    k2 = rjs_idx(j)
    exp_coeff_der(1:alpha_max, k2) = amplitudes(j) * exp_coeff_soft_der_array(j, 1:alpha_max) + &
                                     amplitudes_der(j) * exp_coeff_soft_array(j, 1:alpha_max)
  end do

  deallocate( g_aux_left_der_array, g_aux_right_der_array, M_left_der_array, M_right_der_array, & 
              lower_limit_der_array, upper_limit_der_array, poly_left, poly_right, I0_block, &
              I_left_der_array, I_right_der_array, exp_coeff_soft_der_array, mask_der_array )

end if


deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_soft_array, soft_weights, &
            I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
            I_left_array, I_right_array, exp_coeff_soft_array, amplitudes, amplitudes_der, B)

else

deallocate( atom_widths )

end if









! BUFFER REGION
! Temporarily allocate this as oversized array:
allocate( atom_widths(1:n_neigh(i)) )
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*(atom_sigma_in + atom_sigma_scaling*rjs_in(k+1:k+n_neigh(i)))

! count number of atoms within buffer region
nn = count( rcut_soft_in < rcut_hard_in .and. rjs_in(k+1:k+n_neigh(i)) + atom_widths > rcut_soft_in .and. mask(k+1:k+n_neigh(i)) )

if( nn > 0 )then

! These need to be allocated immediately
allocate( rjs(1:nn) )
allocate( rjs_idx(1:nn) )

! make assignment between original neighbors and reduced neighbors
k2 = 1
do j = 1, n_neigh(i)
!  if( rjs_in(k+j) < rcut_hard_in .and. mask(k+j) )then
  if( rcut_soft_in < rcut_hard_in .and. rjs_in(k+j) + atom_widths(j) > rcut_soft_in .and. mask(k+j) )then
    rjs(k2) = rjs_in(k+j)
    rjs_idx(k2) = k+j
    k2 = k2 + 1
  end if
end do

deallocate( atom_widths )




! allocate arrays that depend on the number of neighbors inside soap cutoff
allocate( atom_sigma_scaleds(1:nn) )
allocate( s2s(1:nn) )
allocate( atom_widths(1:nn) )
allocate( lim_buffer_array(1:nn, 1:3) )
!allocate( soft_weights(1:nn) )
allocate( buffer_weights(1:nn) )
! we might be able to get some extra speed by swapping the order of some of these dimensions
allocate( I0_array(1:nn, 1:alpha_max, 1:3) )
allocate( g_aux_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( g_aux_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( M_left_array(1:nn, 1:alpha_max, 1:2) )
allocate( M_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( I_left_array(1:nn, 1:alpha_max) )
allocate( I_right_array(1:nn, 1:alpha_max) )
!allocate( exp_coeff_soft_array(1:nn, 1:alpha_max) )
allocate( exp_coeff_buffer_array(1:nn, 1:alpha_max) )
allocate( amplitudes(1:nn) )
allocate( amplitudes_der(1:nn) )
!allocate( B(1:nn, 1:7) )
allocate( B(1:7, 1:nn) )

if( do_derivatives )then
allocate( B_der(1:7, 1:nn) )
allocate( lower_limit_der_array(1:nn, 1:1) )
allocate( upper_limit_der_array(1:nn, 1:1) )
allocate( g_aux_left_der_array(1:nn, 1:alpha_max, 1:2) )
allocate( g_aux_right_der_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( M_left_der_array(1:nn, 1:alpha_max, 1:2) )
allocate( M_right_der_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
allocate( I_left_der_array(1:nn, 1:alpha_max) )
allocate( I_right_der_array(1:nn, 1:alpha_max) )
allocate( poly_left(1:nn, 1:2) )
allocate( poly_right(1:nn, 1:2) )
allocate( I0_block(1:nn, 1:alpha_max) )
allocate( exp_coeff_buffer_der_array(1:nn, 1:alpha_max) )
end if


! define atom parameters in vector form
!   **************** New basis ******************
rjs = rjs/rcut_hard_in
!   *********************************************
atom_sigma_scaleds = atom_sigma + atom_sigma_scaling*rjs
s2s = atom_sigma_scaleds**2
atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaleds
atom_width_scaling = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaling
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitudes = 1.d0 / atom_sigma_scaleds
              amplitudes_der = - atom_sigma_scaling / s2s 
!            else if( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 <= 1.d-10 )then ! FIX THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              amplitudes = 0.d0
!              amplitude_ders = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )
                amplitudes_der = 6.d0 / atom_sigma_scaleds * (rjs**2 - rjs) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              else
                amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**amplitude_scaling
                amplitudes_der = 6.d0*amplitude_scaling / atom_sigma_scaleds * (rjs**2 - rjs) &
                                * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( size(amplitudes) > 0 )then
            if( rjs_idx(1) == k+1 )then
              amplitudes(1) = central_weight * amplitudes(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              amplitudes_der(1) = central_weight * amplitudes_der(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            end if
          end if

!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitudes_der = amplitudes * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                             amplitudes_der * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
            amplitudes = amplitudes * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
          else if( radial_enhancement == 2 )then
            amplitudes_der = amplitudes*( 2.d0*rjs + 2.d0*atom_sigma_scaleds*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaleds + dsqrt(8.d0/pi)*rjs*atom_sigma_scaling ) + &
                            amplitudes_der*( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
            amplitudes = amplitudes * ( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
          end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! NOTE: this below is design to be used with arrays that contain all the atoms. In some cases, the left_weights or
! right_weights might contain lots of zeros, in which case a speed up of close to x2 could be achieved by not doing
! all atoms but only those for which the weights are 1. One could refactor this code with the nn construction to
! have, e.g., nn_soft and nn_buffer


! define arrays that give contribution of atoms within left and right of rcut_soft
! note that these are integer arrays (0|1 = .false.|.true.)
!soft_weights = rjs - atom_widths < rcut_soft
!nn_soft = count(soft_weights)
!buffer_weights = rcut_soft < rcut_hard .and. rjs + atom_widths > rcut_soft
!nn_buffer = count(buffer_weights)




! handle the central atom contribution; NOTE: I'M NOT SURE ABOUT THE NEED FOR THIS DO_CENTRAL TAG HERE, IT SHOULD BE HANDLED VIA THE CENTRAL_WEIGHT
if( .not. do_central )then
  if( size(amplitudes) > 0 )then
    if( rjs_idx(1) == k+1 )then
      amplitudes(1) = 0.d0
    end if
  end if
end if


! Atoms within the buffer region

lim_buffer_array(:, 1) = max( rcut_soft, rjs - atom_widths ) ! lower limit left
lim_buffer_array(:, 2) = max( rjs, rcut_soft ) ! upper limit left / lower limit right
lim_buffer_array(:, 3) = min( rcut_hard, rjs + atom_widths ) ! upper limit right
I0_array = M_radial_poly_array(lim_buffer_array, max(7, alpha_max + 4), rcut_hard)

call get_constant_poly_filter_coeff_array(rjs, atom_widths, rcut_soft, filter_width, 'left', B)

! We should try to figure out a more "vectorized way" of doing this
do k2 = 1, nn
M_left_array(k2, 1:7, 1) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 1), 6) ) * I0_array(k2, 1:7, 1)
M_left_array(k2, 1:7, 2) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * I0_array(k2, 1:7, 2)
end do

I_left_array(1:nn, 1:alpha_max) = matmul( M_left_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                  matmul( M_left_array(1:nn, 1:7, 1), transpose(A(1:alpha_max, 1:7)) ) * &
                                  I0_array(1:nn, 5:alpha_max + 4, 1)

call get_constant_poly_filter_coeff_array(rjs, atom_widths, rcut_soft, filter_width, 'right', B)

! We should try to figure out a more "vectorized way" of doing this
do k2 = 1, nn
M_right_array(k2, 1:7, 2) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * I0_array(k2, 1:7, 2)
M_right_array(k2, 1:7, 3) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 3), 6) ) * I0_array(k2, 1:7, 3)
end do

I_right_array(1:nn, 1:alpha_max) = matmul( M_right_array(1:nn, 1:7, 3), transpose(A(1:alpha_max, 1:7)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                   matmul( M_right_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                   I0_array(1:nn, 5:alpha_max + 4, 2)

exp_coeff_buffer_array = I_left_array + I_right_array


do j = 1, nn
  k2 = rjs_idx(j)
  exp_coeff(1:alpha_max, k2) = exp_coeff(1:alpha_max, k2) + amplitudes(j) * exp_coeff_buffer_array(j, 1:alpha_max)
end do

if( do_derivatives )then
! contribution to derivatives --> (3) in notes
  call get_constant_poly_filter_coeff_der_array(rjs, atom_widths, atom_width_scaling, rcut_soft, filter_width, 'left', B_der)

! We should try to figure out a more "vectorized way" of doing this
  do k2 = 1, nn
    M_left_der_array(k2, 1:7, 1) = matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 1), 6) ) * I0_array(k2, 1:7, 1)
    M_left_der_array(k2, 1:7, 2) = matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * I0_array(k2, 1:7, 2)
  end do
  I_left_der_array(1:nn, 1:alpha_max) = matmul( M_left_der_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                        I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                        matmul( M_left_der_array(1:nn, 1:7, 1), transpose(A(1:alpha_max, 1:7)) ) * &
                                        I0_array(1:nn, 5:alpha_max + 4, 1)

  call get_constant_poly_filter_coeff_der_array(rjs, atom_widths, atom_width_scaling, rcut_soft, filter_width, 'right', B_der)

! We should try to figure out a more "vectorized way" of doing this
  do k2 = 1, nn
    M_right_der_array(k2, 1:7, 2) = matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * I0_array(k2, 1:7, 2)
    M_right_der_array(k2, 1:7, 3) = matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 3), 6) ) * I0_array(k2, 1:7, 3)
  end do
  I_right_der_array(1:nn, 1:alpha_max) = matmul( M_right_der_array(1:nn, 1:7, 3), transpose(A(1:alpha_max, 1:7)) ) * &
                                         I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                         matmul( M_right_der_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                         I0_array(1:nn, 5:alpha_max + 4, 2)

  exp_coeff_buffer_der_array = I_left_der_array + I_right_der_array

! add derivatives of the integral limits if necessary --> (1), (2) in notes
! check dimensions!!!!!!!!!!!!!!!!
  mask_der_array = (lim_buffer_array(:, 1) == rjs - atom_widths)
  nn_lower_lim_der = count(mask_der_array)
  if( nn_lower_lim_der > 0 )then
    lower_limit_der_array(1:nn, 1) = merge(1, 0, mask_der_array)
    g_aux_left_der_array(1:nn, 1:4, 1:2) = g_aux_der_array(lim_buffer_array(:, 1:2), rjs(:), atom_widths(:), &
                                                           atom_width_scaling, "left")
    poly_left = g_aux_left_der_array(1:nn, 1, 1:1)
    I0_block = I0_array(1:nn, 4:alpha_max + 4, 1)
    exp_coeff_buffer_der_array = exp_coeff_buffer_der_array - (1.d0 - atom_width_scaling) * &
                                 matmul( lower_limit_der_array * poly_left, I_ones ) * I0_block
  end if

  mask_der_array = (lim_buffer_array(:, 3) == rjs + atom_widths)
  nn_upper_lim_der = count(mask_der_array)
  if( nn_upper_lim_der > 0 )then
    upper_limit_der_array(1:nn, 1) = merge(1, 0, mask_der_array)
    g_aux_right_der_array(1:nn, 1:4, 2:3) = g_aux_der_array(lim_buffer_array(:, 2:3), rjs(:), atom_widths(:), &
                                                            atom_width_scaling, "right")
    poly_right = g_aux_right_der_array(1:nn, 1, 3:3)
    I0_block = I0_array(1:nn, 4:alpha_max + 4, 3)
    exp_coeff_buffer_der_array = exp_coeff_buffer_der_array + (1.d0 + atom_width_scaling) * &
                                 matmul( upper_limit_der_array * poly_right, I_ones ) * I0_block
  end if


  do j = 1, nn
    k2 = rjs_idx(j)
    exp_coeff_der(1:alpha_max, k2) = exp_coeff_der(1:alpha_max, k2) + &
                                     amplitudes(j) * exp_coeff_buffer_der_array(j, 1:alpha_max) + &
                                     amplitudes_der(j) * exp_coeff_buffer_array(j, 1:alpha_max)
  end do



  deallocate( g_aux_left_der_array, g_aux_right_der_array, M_left_der_array, M_right_der_array, &
              lower_limit_der_array, upper_limit_der_array, poly_left, poly_right, I0_block, & 
              I_left_der_array, I_right_der_array, exp_coeff_buffer_der_array, B_der, mask_der_array )

end if


deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_buffer_array, buffer_weights, &
            I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
            I_left_array, I_right_array, amplitudes, amplitudes_der, B, exp_coeff_buffer_array )

else

deallocate( atom_widths )

end if





! Transform from g_alpha to g_n (the orthonormal basis)
!exp_coeff_soft_array(1:nn, 1:alpha_max) = matmul( exp_coeff_soft_array(1:nn, 1:alpha_max), W )

!exp_coeff_buffer_array(1:nn, 1:alpha_max) = matmul( exp_coeff_buffer_array(1:nn, 1:alpha_max), W )



!do j = 1, nn
!  k2 = rjs_idx(j)
!  exp_coeff(1:alpha_max, k2) = amplitudes(j) * soft_weights(j) * exp_coeff_soft_array(j, 1:alpha_max) + &
!                               amplitudes(j) * buffer_weights(j) * exp_coeff_buffer_array(j, 1:alpha_max)
!end do


! Transform from g_alpha to g_n (the orthonormal basis)
exp_coeff(1:alpha_max, k+1:k+n_neigh(i)) = matmul( W, exp_coeff(1:alpha_max, k+1:k+n_neigh(i)) )
if( do_derivatives )then
  exp_coeff_der(1:alpha_max, k+1:k+n_neigh(i)) = matmul( W, exp_coeff_der(1:alpha_max, k+1:k+n_neigh(i)) )
end if

k = k + n_neigh(i)

!deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_soft_array, soft_weights, &
!            buffer_weights, I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
!            I_left_array, I_right_array, exp_coeff_soft_array, amplitudes, amplitudes_der, B, &
!            exp_coeff_buffer_array )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do

!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
    if( do_derivatives )then
      exp_coeff_der = exp_coeff_der / dsqrt(rcut_hard_in)
    end if
!   *********************************************

!   This is for debugging
    if( .true. )then
      open(10, file="radial_expansion_coefficients_operator-fast-der.dat", status="unknown")
      do i=1, size(exp_coeff, 2)
        write(10,*) exp_coeff(1:alpha_max, i)
      end do
      close(10)
      if( do_derivatives )then
        open(10, file="radial_expansion_coefficients_der_operator-fast-der.dat", status="unknown")
        do i=1, size(exp_coeff, 2)
          write(10,*) exp_coeff_der(1:alpha_max, i)
        end do
        close(10)
      end if
    end if

    deallocate( A )

  return
  end subroutine
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and a polynomial piecewise representation for 
! the atomic density
!
  subroutine get_radial_expansion_coefficients_poly3operator(n_sites, n_neigh, rjs_in, alpha_max, &
                                                             rcut_soft_in, rcut_hard_in, atom_sigma_in, &
                                                             atom_sigma_scaling, amplitude_scaling, W, &
                                                             scaling_mode, mask, radial_enhancement, &
                                                             do_derivatives, do_central, central_weight, &
                                                             exp_coeff, exp_coeff_der)

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling, central_weight
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
!
    integer :: i, j, k
    real*8 :: pi, rj, s2, atom_width, atom_width_scaling, filter_width, x
    real*8 :: lim_soft(1:3), lim_buffer(1:3), B(1:7)
    real*8 :: M_left(1:7, 1:2), M_right(1:7, 1:2)
    real*8, allocatable :: A(:,:), I0(:,:), I_left(:), I_right(:)
    real*8 :: W(:,:)
!   derivatives
    real*8 :: amplitude_der
    real*8 :: poly_left(1:4), poly_right(1:4), poly_filter(1:4), B_der(1:7)
    real*8 :: M_left_der(1:7, 1:2), M_right_der(1:7, 1:2)
    real*8, allocatable :: I_left_der(:), I_right_der(:), upper_limit_der(:), lower_limit_der(:)
!   Results will be stored in exp_coeff, which is an array of dimension (alpha_max, n_atom_pairs)
    real*8 :: exp_coeff(:,:), exp_coeff_der(:,:)
    real*8, allocatable :: exp_coeff_soft(:), exp_coeff_buffer(:)
    real*8, allocatable :: exp_coeff_soft_der(:), exp_coeff_buffer_der(:) 
    logical, save :: print_basis = .false.


    allocate( exp_coeff_soft(1:alpha_max) )
    allocate( exp_coeff_buffer(1:alpha_max) )
    allocate( A(1:alpha_max, 1:7) )
    allocate( I0(1:alpha_max + 4, 1:3) )
    allocate( I_left(1:alpha_max) )
    allocate( I_right(1:alpha_max) )

    allocate( exp_coeff_soft_der(1:alpha_max) )
    allocate( exp_coeff_buffer_der(1:alpha_max) )
    allocate( I_left_der(1:alpha_max) )
    allocate( I_right_der(1:alpha_max) )
    allocate( lower_limit_der(1:alpha_max) )
    allocate( upper_limit_der(1:alpha_max) )


!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      do j=1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,*) "0."
      end do
    end if

!   *************** New basis *******************
!   Normalized with rcut_hard to avoid instabilities
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
!   *********************************************
    pi = dacos(-1.d0)
    filter_width = 2.d0*sqrt(2.d0*log(2.d0))*(rcut_hard - rcut_soft)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if

    call get_constant_poly_coeff(alpha_max, rcut_hard, A)

    k = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
!       Check if we need to do the central atom
        if( j == 1 .and. .not. do_central )then
          cycle
        end if
        if( rjs_in(k) < rcut_hard_in .and. mask(k) )then
!   **************** New basis ******************
          rj = rjs_in(k)/rcut_hard_in
!   *********************************************
          atom_sigma_scaled = atom_sigma + atom_sigma_scaling*rj
          s2 = atom_sigma_scaled**2
          atom_width = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaled
          atom_width_scaling = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaling
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
          if( scaling_mode == "polynomial" )then
!           WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!           be taken into account in the corresponding part of the code.
!       
!           WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!           from the expressions
            if( amplitude_scaling == 0.d0 )then
              amplitude = 1.d0 / atom_sigma_scaled
              amplitude_der = - atom_sigma_scaling / s2 
            else if( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 <= 1.d-10 )then
              amplitude = 0.d0
              amplitude_der = 0.d0
            else
              if( amplitude_scaling == 1.d0 )then
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )
                amplitude_der = 6.d0 / atom_sigma_scaled * (rj**2 - rj) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              else
                amplitude = 1.d0 / atom_sigma_scaled * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**amplitude_scaling
                amplitude_der = 6.d0*amplitude_scaling / atom_sigma_scaled * (rj**2 - rj) &
                                * ( 1.d0 + 2.d0*rj**3 - 3.d0*rj**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaled * amplitude
              end if
            end if
          end if
!         The central atom needs to be scaled by central_weight
          if( j == 1 )then
            amplitude = central_weight * amplitude
            amplitude_der = central_weight * amplitude_der
          end if
!         The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!         of atom j.
          if( radial_enhancement == 1 )then
            amplitude_der = amplitude * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                            amplitude_der * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
            amplitude = amplitude * ( rj + dsqrt(2.d0/pi)*atom_sigma_scaled )
          else if( radial_enhancement == 2 )then
            amplitude_der = amplitude*( 2.d0*rj + 2.d0*atom_sigma_scaled*atom_sigma_scaling + &
                                        dsqrt(8.d0/pi)*atom_sigma_scaled + dsqrt(8.d0/pi)*rj*atom_sigma_scaling ) + &
                            amplitude_der*( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
            amplitude = amplitude * ( rj**2 + s2 + dsqrt(8.d0/pi)*atom_sigma_scaled*rj )
          end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          exp_coeff_soft = 0.d0
          exp_coeff_buffer = 0.d0
          exp_coeff_soft_der = 0.d0
          exp_coeff_buffer_der = 0.d0
          lower_limit_der = 0.d0
          upper_limit_der = 0.d0
!         contribution inside rcut_soft from left and right part of the piecewise polynomial
          lim_soft = 0.d0
          if ( rj - atom_width < rcut_soft )then
!           integration limits inside r_soft
            lim_soft(1) = max( 0.d0, rj - atom_width )      ! lower limit left
            lim_soft(2) = min( max(0.d0, rj), rcut_soft )   ! upper limit left / lower limit right
            lim_soft(3) = min( rcut_soft, rj + atom_width ) ! upper limit right
!            write(*,'(A)') '- Inner zone:'
!            print *, lim_soft
!           contribution to the expansion coeff.
            I0 = transpose( M_radial_poly(lim_soft, alpha_max + 4, rcut_hard) ) ! +4 because we include poly alpha_max
            M_left(1:4, 1) = I0(1:4, 1) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(1), rj, atom_width, "left" )
            M_left(1:4, 2) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(2), rj, atom_width, "left" )
            I_left = matmul( A(1:alpha_max, 1:4), M_left(1:4, 2) ) * I0(5:alpha_max + 4, 2) - &
                     matmul( A(1:alpha_max, 1:4), M_left(1:4, 1) ) * I0(5:alpha_max + 4, 1)
            M_right(1:4, 1) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(2), rj, atom_width, "right" )
            M_right(1:4, 2) = I0(1:4, 3) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(3), rj, atom_width, "right" )
            I_right = matmul( A(1:alpha_max, 1:4), M_right(1:4, 2) ) * I0(5:alpha_max + 4, 3) - &
                      matmul( A(1:alpha_max, 1:4), M_right(1:4, 1) ) * I0(5:alpha_max + 4, 2)
            exp_coeff_soft = I_left + I_right
            if( do_derivatives )then
!             contribution to derivatives --> (3) in notes
              M_left_der(1:4, 1) = I0(1:4, 1) * [-1.d0, -1.d0, -2.d0, -6.d0] * &
                                   g_aux_der( lim_soft(1), rj, atom_width, atom_width_scaling, "left" )
              M_left_der(1:4, 2) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * &
                                   g_aux_der( lim_soft(2), rj, atom_width, atom_width_scaling, "left" )
              I_left_der = matmul( A(1:alpha_max, 1:4), M_left_der(1:4, 2) ) * I0(5:alpha_max + 4, 2) - &
                           matmul( A(1:alpha_max, 1:4), M_left_der(1:4, 1) ) * I0(5:alpha_max + 4, 1)
              M_right_der(1:4, 1) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * &
                                    g_aux_der( lim_soft(2), rj, atom_width, atom_width_scaling, "right" )
              M_right_der(1:4, 2) = I0(1:4, 3) * [-1.d0, -1.d0, -2.d0, -6.d0] * &
                                    g_aux_der( lim_soft(3), rj, atom_width, atom_width_scaling, "right" )
              I_right_der = matmul( A(1:alpha_max, 1:4), M_right_der(1:4, 2) ) * I0(5:alpha_max + 4, 3) - &
                            matmul( A(1:alpha_max, 1:4), M_right_der(1:4, 1) ) * I0(5:alpha_max + 4, 2)
              exp_coeff_soft_der = I_left_der + I_right_der
!             add derivatives of the integral limits if necessary --> (1), (2) in notes
              if ( lim_soft(1) == rj - atom_width )then
                x = (lim_soft(1) - rj)/atom_width
                poly_left = 1.d0 - 3.d0*x**2 - 2.d0*x**3
                lower_limit_der = poly_left * (1.d0 - atom_width_scaling) * I0(4:alpha_max + 4, 1)
              end if
              if ( lim_soft(3) == rj + atom_width )then
                x = (lim_soft(3) - rj)/atom_width
                poly_right = 1.d0 - 3.d0*x**2 + 2.d0*x**3
                upper_limit_der = poly_right * (1.d0 + atom_width_scaling) * I0(4:alpha_max + 4, 3)
              end if 
            end if
          end if
!         contribution in the buffer zone from left and right part of the piecewise polynomial
          lim_buffer = 0.d0
          if ( rcut_soft < rcut_hard .and. rj + atom_width > rcut_soft )then
!           integration limits inside r_buffer
            lim_buffer(1) = max( rcut_soft, rj - atom_width ) ! lower limit left
            lim_buffer(2) = max( rj, rcut_soft )                     ! upper limit left / lower limit right
            lim_buffer(3) = min( rcut_hard, rj + atom_width ) ! upper limit right
!            write(*,'(A)') '* Buffer zone:' 
!           contribution to the expansion coeff.
            I0 = transpose( M_radial_poly(lim_buffer, max(7, alpha_max + 4), rcut_hard) )
            call get_constant_poly_filter_coeff(rj, atom_width, rcut_soft, filter_width, 'left', B)
            M_left(1:7, 1) = matmul( -B(1:7), M_radial_monomial(lim_buffer(1), 6) ) * I0(1:7, 1)
            M_left(1:7, 2) = matmul( -B(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
            I_left = matmul( A(1:alpha_max, 1:7), M_left(1:7, 2) ) * I0(5:alpha_max + 4, 2) - &
                     matmul( A(1:alpha_max, 1:7), M_left(1:7, 1) ) * I0(5:alpha_max + 4, 1)
            call get_constant_poly_filter_coeff(rj, atom_width, rcut_soft, filter_width, 'right', B)
            M_right(1:7, 1) = matmul( -B(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
            M_right(1:7, 2) = matmul( -B(1:7), M_radial_monomial(lim_buffer(3), 6) ) * I0(1:7, 3)
            I_right = matmul( A(1:alpha_max, 1:7), M_right(1:7, 2) ) * I0(5:alpha_max + 4, 3) - &
                      matmul( A(1:alpha_max, 1:7), M_right(1:7, 1) ) * I0(5:alpha_max + 4, 2)
            exp_coeff_buffer = I_left + I_right
!           contribution to derivatives --> (3) in notes
            if( do_derivatives )then
!             contribution to derivatives --> (3) in notes
              call get_constant_poly_filter_coeff_der(rj, atom_width, atom_width_scaling, rcut_soft, filter_width, &
                                                      'left', B_der)
              M_left_der(1:7, 1) = matmul( -B_der(1:7), M_radial_monomial(lim_buffer(1), 6) ) * I0(1:7, 1)
              M_left_der(1:7, 2) = matmul( -B_der(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
              I_left_der = matmul( A(1:alpha_max, 1:7), M_left_der(1:7, 2) ) * I0(5:alpha_max + 4, 2) - &
                           matmul( A(1:alpha_max, 1:7), M_left_der(1:7, 1) ) * I0(5:alpha_max + 4, 1)
              call get_constant_poly_filter_coeff_der(rj, atom_width, atom_width_scaling, rcut_soft, filter_width, &
                                                      'right', B_der)
              M_right_der(1:7, 1) = matmul( -B_der(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
              M_right_der(1:7, 2) = matmul( -B_der(1:7), M_radial_monomial(lim_buffer(3), 6) ) * I0(1:7, 3)
              I_right_der = matmul( A(1:alpha_max, 1:7), M_right_der(1:7, 2) ) * I0(5:alpha_max + 4, 3) - &
                            matmul( A(1:alpha_max, 1:7), M_right_der(1:7, 1) ) * I0(5:alpha_max + 4, 2)
              exp_coeff_buffer_der = I_left_der + I_right_der
!             add derivatives of the integral limits if necessary --> (1), (2) in notes
              if ( lim_buffer(1) == rj - atom_width )then
                x = (lim_buffer(1) - rj)/atom_width
                poly_left = 1.d0 - 3.d0*x**2 - 2.d0*x**3
                x = (lim_buffer(1) - rcut_soft)/filter_width
                poly_filter = 1.d0 - 3.d0*x**2 + 2.d0*x**3
                lower_limit_der = poly_left * (1.d0 - atom_width_scaling) * poly_filter * I0(4:alpha_max + 4, 1)
              end if
              if ( lim_buffer(3) == rj + atom_width )then
                x = (lim_buffer(3) - rj)/atom_width
                poly_right = 1.d0 - 3.d0*x**2 + 2.d0*x**3
                x = (lim_buffer(3) - rcut_soft)/filter_width
                poly_filter = 1.d0 - 3.d0*x**2 + 2.d0*x**3
                upper_limit_der = poly_right * (1.d0 + atom_width_scaling) * poly_filter * I0(4:alpha_max + 4, 3)
              end if 
            end if
          end if
!         Transform from g_alpha to g_n (the orthonormal basis)
          exp_coeff(1:alpha_max, k) = amplitude * matmul( W, exp_coeff_soft(1:alpha_max) + exp_coeff_buffer(1:alpha_max) )
          if( do_derivatives )then
            exp_coeff_der(1:alpha_max, k) = amplitude * ( exp_coeff_soft_der(1:alpha_max) + upper_limit_der + &
                                            exp_coeff_buffer_der(1:alpha_max) - lower_limit_der ) + &
                                            amplitude_der * ( exp_coeff_soft(1:alpha_max) + exp_coeff_buffer(1:alpha_max) )
            exp_coeff_der(1:alpha_max, k) = matmul( W,  exp_coeff_der(1:alpha_max, k))
          end if
        end if
      end do
    end do

!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
    if( do_derivatives )then
      exp_coeff_der = exp_coeff_der / dsqrt(rcut_hard_in)
    end if
!   *********************************************

!   This is for debugging
    if( .true. )then
      open(10, file="radial_expansion_coefficients_operator.dat", status="unknown")
      do i=1, size(exp_coeff, 2)
        write(10,*) exp_coeff(1:alpha_max, i)
      end do
      close(10)
      if( do_derivatives )then
        open(10, file="radial_expansion_coefficients_der_operator.dat", status="unknown")
        do i=1, size(exp_coeff, 2)
          write(10,*) exp_coeff_der(1:alpha_max, i)
        end do
        close(10)
      end if
    end if

    deallocate( exp_coeff_soft, exp_coeff_buffer, A, I0, I_left, I_right )
    deallocate( exp_coeff_soft_der, exp_coeff_buffer_der )

  return
  end subroutine
!**************************************************************************


!**************************************************************************
!
! Auxiliary function that computes a piecewise polynomial function and
! its derivatives with respect to r
!
  function g_aux(r, r0, width, piece) result(poly)
    implicit none
 
    real*8, intent(in) :: r, r0, width
    character(*), intent(in) :: piece
    real*8 :: x
    real*8 :: poly(1:4)

    x = (r - r0)/width

    if ( piece == "left" )then
      poly = [ 1.d0 - 3.d0*x**2 - 2.d0*x**3, -6.d0*(x**2 + x)/width, -3.d0*(2.d0*x + 1)/width**2, -2.d0/width**3 ]
    else if ( piece == "right" )then
      poly = [ 1.d0 - 3.d0*x**2 + 2.d0*x**3, 6.d0*(x**2 - x)/width, 3.d0*(2.d0*x - 1)/width**2, 2.d0/width**3 ]
    end if

  return
  end function
!**************************************************************************

!**************************************************************************
!
! Auxiliary function that computes a piecewise polynomial function and
! its derivatives with respect to r
!
  function g_aux_array(r, r0, width, piece) result(poly)
    implicit none
 
    real*8, intent(in) :: r(:,:), r0(:), width(:)
    character(*), intent(in) :: piece
    real*8, dimension( 1:size(r,1), 1:size(r,2) ) :: x, x2, x3
    real*8, dimension( 1:size(width) ) :: w, w2, w3
    real*8, dimension( 1:size(r,1), 1:4, 1:size(r,2) ) :: poly
    integer :: i, j

    w(:) = 1.d0/width(:)

    do i = 1, size(r,2)
!      x(:, i)  = (r(:, i) - r0(:))/width(:)
      x(:, i)  = (r(:, i) - r0(:)) * w
    end do
    x2 = x*x
    x3 = x2*x

    w2 = w*w
    w3 = w2*w

    if ( piece == "left" )then
      do i = 1, size(r,2)
!        poly(:, 1, i) = 1.d0 - 3.d0*x(:, i)**2 - 2.d0*x(:, i)**3
!        poly(:, 2, i) = -6.d0*(x(:, i)**2 + x(:, i))/width(:)
!        poly(:, 3, i) = -3.d0*(2.d0*x(:, i) + 1)/width(:)**2
!        poly(:, 4, i) = -2.d0/width(:)**3
        poly(:, 1, i) = 1.d0 - 3.d0*x2(:, i) - 2.d0*x3(:, i)
        poly(:, 2, i) = -6.d0*(x2(:, i) + x(:, i)) * w
        poly(:, 3, i) = -3.d0*(2.d0*x(:, i) + 1) * w2
        poly(:, 4, i) = -2.d0 * w3
      end do
    else if ( piece == "right" )then
      do i = 1, size(r,2)
        poly(:, 1, i) = 1.d0 - 3.d0*x2(:, i) + 2.d0*x3(:, i)
        poly(:, 2, i) = 6.d0*(x2(:, i) - x(:, i)) * w
        poly(:, 3, i) = 3.d0*(2.d0*x(:, i) - 1) * w2
        poly(:, 4, i) = 2.d0 * w3
      end do
    end if

  return
  end function
!**************************************************************************

!**************************************************************************
!
! Auxiliary function derivatives of g_aux() with respect to r0, where
!                width = width0 + width_scaling * r0
! and
!   width_scaling = 2.d0*sqrt(2.d0*log(2.d0)) * atom_sigma_scaling
!
  function g_aux_der(r, r0, width, width_scaling, piece) result(poly)
    implicit none
 
    real*8, intent(in) :: r, r0, width, width_scaling
    character(*), intent(in) :: piece
    real*8 :: x
    real*8 :: poly(1:4)

    x = (r - r0)/width

    if ( piece == "left" )then
      poly = [ 6.d0*(x + (width_scaling + 1.d0)*x**2 + width_scaling*x**3)/width,              &
               6.d0*(1.d0 + 2.d0*(width_scaling + 1.d0)*x + 3.d0*width_scaling*x**2)/width**2, &
               6.d0*(width_scaling + 1.d0 + 3.d0*width_scaling*x)/width**3,                       &
               6.d0*width_scaling/width**4 ]
    else if ( piece == "right" )then
      poly = [ 6.d0*(x + (width_scaling - 1.d0)*x**2 - width_scaling*x**3)/width,              &
               6.d0*(1.d0 + 2.d0*(width_scaling - 1.d0)*x - 3.d0*width_scaling*x**2)/width**2, &
               6.d0*(width_scaling - 1.d0 - 3.d0*width_scaling*x)/width**3,                       &
              -6.d0*width_scaling/width**4 ]
    end if

  return
  end function
!**************************************************************************

!**************************************************************************
!
! Auxiliary function derivatives of g_aux() with respect to r0, where
!                width = width0 + width_scaling * r0
! and
!   width_scaling = 2.d0*sqrt(2.d0*log(2.d0)) * atom_sigma_scaling
!
  function g_aux_der_array(r, r0, width, width_scaling, piece) result(poly)
    implicit none
 
    real*8, intent(in) :: r(:,:), r0(:), width(:), width_scaling
    character(*), intent(in) :: piece
    real*8, dimension( 1:size(r,1), 1:size(r,2) ) :: x, x2, x3
    real*8, dimension( 1:size(width) ) :: w, w2, w3, w4
    real*8, dimension( 1:size(r,1), 1:4, 1:size(r,2) ) :: poly
    integer :: i, j

    w(:) = 1.d0/width(:)

    do i = 1, size(r,2)
      x(:, i)  = (r(:, i) - r0(:)) * w
    end do
    x2 = x*x
    x3 = x2*x

    w2 = w*w
    w3 = w2*w
    w4 = w3*w

    if ( piece == "left" )then
      do i = 1, size(r,2)
        poly(:, 1, i) = 6.d0*(x(:, i) + (width_scaling + 1.d0)*x2(:, i) + width_scaling*x3(:, i)) * w
        poly(:, 2, i) = 6.d0*(1.d0 + 2.d0*(width_scaling + 1.d0)*x(:, i) + 3.d0*width_scaling*x2(:, i)) * w2
        poly(:, 3, i) = 6.d0*(width_scaling + 1.d0 + 3.d0*width_scaling*x(:, i)) * w3
        poly(:, 4, i) = 6.d0*width_scaling * w4
      end do
    else if ( piece == "right" )then
      do i = 1, size(r,2)
        poly(:, 1, i) = 6.d0*(x(:, i) + (width_scaling - 1.d0)*x2(:, i) - width_scaling*x3(:, i)) * w
        poly(:, 2, i) = 6.d0*(1.d0 + 2.d0*(width_scaling - 1.d0)*x(:, i) - 3.d0*width_scaling*x2(:, i)) * w2
        poly(:, 3, i) = 6.d0*(width_scaling - 1.d0 - 3.d0*width_scaling*x(:, i)) * w3
        poly(:, 4, i) = -6.d0*width_scaling * w4
      end do
    end if


  return
  end function
!**************************************************************************



!**************************************************************************
!
! This subroutine returns the matrix containing the atom-independent
! coefficients for the polynomial radial basis.
! (in the notes: A, A*)
!
  subroutine get_constant_poly_coeff(alpha_max, rcut, A)
    implicit none

    integer, intent(in) :: alpha_max
    real*8, intent(in) :: rcut
    integer :: i, j, l(1:2)
    integer, allocatable :: factor(:)
    real*8, intent(inout) :: A(:,:)

    allocate( factor(1:alpha_max) )
    
    do i = 1, alpha_max
      A(i, 1) = rcut / N_a( rcut, i ) / dfloat(i + 3)
      factor(i) = (i + 3)
    end do

    l = shape(A)
    do j = 2, l(2)
        A(:, j) = A(:, j - 1) * rcut / dfloat(factor + j - 1) 
    end do

  end subroutine
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the constant coefficients for the 
! contribution from the buffer zone, coming from the dot products
! of poly. radial basis and smoothing function (filter)
! (in the notes: B_l, B_r)
!
  subroutine get_constant_poly_filter_coeff(rj, width_j, r_filter, filter_width, piece, B)
    implicit none

    real*8, intent(in) :: rj, width_j, r_filter, filter_width
    character(*), intent(in) :: piece
    integer :: i
    real*8 :: C_filter(1:4), col_poly(1:4), C_poly(1:7,1:4)
    real*8, intent(inout) :: B(:)

!   coeff. from the filter
    C_filter = g_aux(0.d0, r_filter, filter_width, "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    C_poly = 0.d0
    col_poly(1:4) = g_aux(0.d0, rj, width_j, piece)

    do i = 1, 4
      C_poly(i:i+3, i) = col_poly
    end do

    B = matmul( C_poly, C_filter )

  end subroutine
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the constant coefficients for the 
! contribution from the buffer zone, coming from the dot products
! of poly. radial basis and smoothing function (filter)
! (in the notes: B_l, B_r)
!
  subroutine get_constant_poly_filter_coeff_array(rj, width_j, r_filter, filter_width, piece, B)
    implicit none

    real*8, intent(in) :: rj(:), width_j(:), r_filter, filter_width
    character(*), intent(in) :: piece
    integer :: i, j, k
    real*8 :: C_filter(1:4)
    real*8, allocatable :: col_poly(:, :, :), C_poly(:, :), rj_temp(:,:)
    real*8, intent(out) :: B(:, :)

    allocate( col_poly(1:size(rj), 1:4, 1:1) )
    allocate( C_poly(1:7*size(rj), 1:4) )
    allocate( rj_temp(1:size(rj), 1:1) )

!   coeff. from the filter
    C_filter = g_aux(0.d0, r_filter, filter_width, "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    C_poly = 0.d0

    rj_temp(:, 1) = -rj(:)

!    col_poly = g_aux_array(-reshape(rj, [size(rj), 1]), 0.d0*rj, width_j, piece)
    col_poly = g_aux_array(rj_temp, 0.d0*rj, width_j, piece)

    do i = 1, 4
      do j = 1, size(rj)
        k = (j-1)*7
        C_poly(k+i:k+i+3, i) = col_poly(j, 1:4, 1)
      end do
    end do

!    B = transpose( reshape(matmul( C_poly, C_filter ), [7, size(rj)]) )
    B = reshape(matmul( C_poly, C_filter ), [7, size(rj)])

    deallocate( col_poly, C_poly, rj_temp )

  end subroutine
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the constant coefficients for the 
! contribution from the buffer zone for the derivatives, coming from
! the dot products of poly. radial basis and smoothing function (filter)
! (in the notes: partial_der(B_l), partial_der(B_r))
!
  subroutine get_constant_poly_filter_coeff_der(rj, width_j, width_scaling, r_filter, filter_width, piece, B)
    implicit none

    real*8, intent(in) :: rj, width_j, r_filter, filter_width, width_scaling
    character(*), intent(in) :: piece
    integer :: i
    real*8 :: C_filter(1:4), col_poly(1:4), C_poly(1:7,1:4)
    real*8, intent(inout) :: B(:)

!   coeff. from the filter
    C_filter = g_aux(0.d0, r_filter, filter_width, "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    C_poly = 0.d0
    col_poly(1:4) = g_aux_der(0.d0, rj, width_j, width_scaling, piece)

    do i = 1, 4
      C_poly(i:i+3, i) = col_poly
    end do

    B = matmul( C_poly, C_filter )

  end subroutine
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the constant coefficients for the 
! contribution from the buffer zone for the derivatives, coming from
! the dot products of poly. radial basis and smoothing function (filter)
! (in the notes: partial_der(B_l), partial_der(B_r))
!
  subroutine get_constant_poly_filter_coeff_der_array(rj, width_j, width_scaling, r_filter, filter_width, piece, B)
    implicit none

    real*8, intent(in) :: rj(:), width_j(:), r_filter, filter_width, width_scaling
    character(*), intent(in) :: piece
    integer :: i, j, k
    real*8 :: C_filter(1:4)
    real*8, allocatable :: col_poly(:, :, :), C_poly(:, :), rj_temp(:,:)
    real*8, intent(out) :: B(:, :)

    allocate( col_poly(1:size(rj), 1:4, 1:1) )
    allocate( C_poly(1:7*size(rj), 1:4) )
    allocate( rj_temp(1:size(rj), 1:1) )

!   coeff. from the filter
    C_filter = g_aux(0.d0, r_filter, filter_width, "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    C_poly = 0.d0

    rj_temp(:, 1) = -rj(:)

    col_poly = g_aux_der_array(rj_temp, 0.d0*rj, width_j, width_scaling, piece)

    do i = 1, 4
      do j = 1, size(rj)
        k = (j-1)*7
        C_poly(k+i:k+i+3, i) = col_poly(j, 1:4, 1)
      end do
    end do

!    B = transpose( reshape(matmul( C_poly, C_filter ), [7, size(rj)]) )
    B = reshape(matmul( C_poly, C_filter ), [7, size(rj)])

    deallocate( col_poly, C_poly, rj_temp )

  end subroutine
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the radial terms coming from the
! polynomial radial basis 
! (in the notes: I0, R_hard)
!
  function M_radial_poly(r, alpha_max, rcut) result(radial_terms)
    implicit none

    integer, intent(in) :: alpha_max
    real*8, intent(in) :: rcut, r(:)
    integer :: i
    real*8, dimension(1:size(r), 1:alpha_max) :: radial_terms

    radial_terms = 1.d0
    do i = 2, alpha_max
      radial_terms(:, i) = radial_terms(:, i - 1) * (1.d0 - r/rcut)
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the radial terms coming from the
! polynomial radial basis 
! (in the notes: I0, R_hard)
!
  function M_radial_poly_array(r, alpha_max, rcut) result(radial_terms)
    implicit none

    integer, intent(in) :: alpha_max
    real*8, intent(in) :: rcut, r(:,:)
    integer :: i, j
!   The 1st dimension of r is neighbor index and the third dimension 1:3 for the 3 limits
    real*8, dimension(1:size(r,1), 1:alpha_max, 1:size(r,2)) :: radial_terms
    real*8, dimension(1:size(r,1), 1:size(r,2)) :: r_slice

    r_slice = 1.d0 - r/rcut

    do j = 1, size(r, 2)
      radial_terms(:, 1, j) = 1.d0
      do i = 2, alpha_max
        radial_terms(:, i, j) = radial_terms(:, i-1, j) * r_slice(:, j)
      end do
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns a vector of monomial terms (1, x, x**2, ...)
! for a given polynomial degree
! (in the notes: R, R_ext) 
!
  function radial_monomial(r, degree) result(radial_terms)
    implicit none

    integer, intent(in) :: degree
    real*8, intent(in) :: r
    integer :: p
    real*8, dimension(1:degree + 1) :: radial_terms 

    radial_terms = 1.d0
    do p = 2, degree + 1
      radial_terms(p) = r * radial_terms(p - 1) 
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns a vector of monomial terms (1, x, x**2, ...)
! for a given polynomial degree
! (in the notes: R, R_ext) 
!
  function radial_monomial_array(r, degree) result(radial_terms)
    implicit none

    integer, intent(in) :: degree
    real*8, intent(in) :: r(:)
    integer :: p
    real*8, dimension(1:size(r), 1:degree + 1) :: radial_terms 

    radial_terms = 1.d0
    do p = 2, degree + 1
      radial_terms(:, p) = r(:) * radial_terms(:, p - 1) 
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns the matrix of monomial terms and their
! derivatives, for a given polynomial degree:
!        ( 1  x  x**2    ...  x**degree              )
!        ( 0  1  2 * x   ...  degree * x**(degree-1) )
!        ( .   .                                     )
!        ( .       .                                 )
!        ( .             .    degree!                )
!
! (in the notes: R*) 
!
!
! NOTE: THIS SHOULD ONLY WORK IF DEGREE = 6, IN ALL OTHER CASES THERE IS A SIZE MISMATCH
! BETWEEN RADIAL_TERMS() AND M_RADIAL()!!!!!!!!!!!!
!
  function M_radial_monomial(r, degree) result(M_radial)
    implicit none

    integer, intent(in) :: degree
    real*8, intent(in) :: r
    real*8, dimension(1:degree + 1) :: radial_terms
    real*8, dimension(1:degree + 1, 1:degree + 1) :: M_radial

    radial_terms = radial_monomial(r, degree)

    M_radial = 0.d0
    M_radial(1:degree + 1, 1) = radial_terms(1:7)
    M_radial(2:degree + 1, 2) = radial_terms(1:6) * [1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0]
    M_radial(3:degree + 1, 3) = radial_terms(1:5) * [2.d0, 6.d0, 12.d0, 20.d0, 30.d0]
    M_radial(4:degree + 1, 4) = radial_terms(1:4) * [6.d0, 24.d0, 60.d0, 120.d0]
    M_radial(5:degree + 1, 5) = radial_terms(1:3) * [24.d0, 120.d0, 360.d0]
    M_radial(6:degree + 1, 6) = radial_terms(1:2) * [120.d0, 720.d0]
    M_radial(7:degree + 1, 7) = 720.d0

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns the matrix of monomial terms and their
! derivatives, for a given polynomial degree:
!        ( 1  x  x**2    ...  x**degree              )
!        ( 0  1  2 * x   ...  degree * x**(degree-1) )
!        ( .   .                                     )
!        ( .       .                                 )
!        ( .             .    degree!                )
!
! (in the notes: R*) 
!
! Hardcoded for degree = 6
!
  function M_radial_monomial_array(r) result(M_radial_t)
    implicit none

    integer :: i
    real*8, intent(in) :: r(:)
    real*8, dimension(1:size(r), 1:7) :: radial_terms
    real*8, dimension(1:size(r), 1:7, 1:7) :: M_radial
    real*8, dimension(1:7, 1:7, 1:size(r)) :: M_radial_t
    real*8, dimension(6) :: v1 = [1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0]
    real*8, dimension(5) :: v2 = [2.d0, 6.d0, 12.d0, 20.d0, 30.d0]
    real*8, dimension(4) :: v3 = [6.d0, 24.d0, 60.d0, 120.d0]
    real*8, dimension(3) :: v4 = [24.d0, 120.d0, 360.d0]
    real*8, dimension(2) :: v5 = [120.d0, 720.d0]

    radial_terms = radial_monomial_array(r, 6)

    M_radial = 0.d0
    M_radial(:, 1:7, 1) = radial_terms(:, 1:7)
    do i = 1, 6
      M_radial(:, i + 1, 2) = radial_terms(:, i) * v1(i)
    end do
    do i = 1, 5
      M_radial(:, i + 2, 3) = radial_terms(:, i) * v2(i)
    end do
    do i = 1, 4
      M_radial(:, i + 3, 4) = radial_terms(:, i) * v3(i)
    end do
    do i = 1, 3
      M_radial(:, i + 4, 5) = radial_terms(:, i) * v4(i)
    end do
    do i = 1, 2
      M_radial(:, i + 5, 6) = radial_terms(:, i) * v5(i)
    end do
    M_radial(:, 7, 7) = 720.d0

    do i = 1, size(r)
      M_radial_t(1:7, 1:7, i) = M_radial(i, 1:7, 1:7)
    end do

  end function
!**************************************************************************


end module soap_turbo_radial_op

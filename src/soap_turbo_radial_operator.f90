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
  subroutine get_radial_exp_coeff_operator_poly3(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                 rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                 amplitude_scaling, W, scaling_mode, mask, &
                                                 radial_enhancement, do_derivatives, do_central, &
                                                 central_weight, exp_coeff)

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real*8, intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, atom_sigma_scaling
    real*8, intent(in) :: amplitude_scaling, central_weight
    real*8 :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
!
    integer :: i, j, k, l
    real*8 :: pi, rj, s2
    real*8 :: lim_soft(1:3), lim_buffer(1:3), B(1:7)
    real*8, allocatable :: A(:,:), I0(:,:), M_left(:,:), M_right(:,:)
    real*8, allocatable :: I_left(:), I_right(:)
    real*8 :: W(:,:)
!   Results will be stored in exp_coeff, which is an array of dimension (alpha_max, n_atom_pairs)
    real*8 :: exp_coeff(:,:)
    real*8, allocatable :: exp_coeff_soft(:), exp_coeff_buffer(:)
    logical, save :: print_basis = .false.
    real*8 :: denom, amplitude_der, pref_f

!   NOTE: the derivatives ARE MISSING !!!!!!!!
    allocate( exp_coeff_soft(1:alpha_max) )
    allocate( exp_coeff_buffer(1:alpha_max) )
    allocate( A(1:alpha_max, 1:7) )
    allocate( I0(1:alpha_max + 4, 1:3) )
    allocate( I_left(1:alpha_max) )
    allocate( I_right(1:alpha_max) )
    allocate( M_left(1:7, 1:2) )
    allocate( M_right(1:7, 1:2) )

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
    exp_coeff = 0.d0


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
!         contribution inside rcut_soft from left and right part of the piecewise polynomial
          lim_soft = 0.d0
          if ( rj - atom_sigma_scaled < rcut_soft )then
!           integration limits inside r_soft
            lim_soft(1) = max( 0.d0, rj - atom_sigma_scaled )         ! lower limit left
            lim_soft(2) = min( max(0.d0, rj), rcut_soft )             ! upper limit left / lower limit right
            lim_soft(3) = min( rcut_soft, rj + atom_sigma_scaled ) ! upper limit right
!            write(*,'(A)') '- Inner zone:'
!            print *, lim_soft
!           contribution to the expansion coeff.
            I0 = transpose( M_radial_poly(lim_soft, alpha_max + 4, rcut_hard) ) ! +4 because we include poly alpha_max
            M_left(1:4, 1) = I0(1:4, 1) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(1), rj, atom_sigma_scaled, "left" )
            M_left(1:4, 2) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(2), rj, atom_sigma_scaled, "left" )
            I_left = matmul( A(1:alpha_max, 1:4), M_left(1:4, 2) ) * I0(5:alpha_max + 4, 2) - &
                     matmul( A(1:alpha_max, 1:4), M_left(1:4, 1) ) * I0(5:alpha_max + 4, 1)
            M_right(1:4, 1) = I0(1:4, 2) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(2), rj, atom_sigma_scaled, "right" )
            M_right(1:4, 2) = I0(1:4, 3) * [-1.d0, -1.d0, -2.d0, -6.d0] * g_aux( lim_soft(3), rj, atom_sigma_scaled, "right" )
            I_right = matmul( A(1:alpha_max, 1:4), M_right(1:4, 2) ) * I0(5:alpha_max + 4, 3) - &
                      matmul( A(1:alpha_max, 1:4), M_right(1:4, 1) ) * I0(5:alpha_max + 4, 2)
            exp_coeff_soft = I_left + I_right
          end if
!         contribution in the buffer zone from left and right part of the piecewise polynomial
          lim_buffer = 0.d0
          if ( rcut_soft < rcut_hard .and. rj + atom_sigma_scaled > rcut_soft )then
!           integration limits inside r_soft
            lim_buffer(1) = max( rcut_soft, rj - atom_sigma_scaled ) ! lower limit left
            lim_buffer(2) = max( rj, rcut_soft )                     ! upper limit left / lower limit right
            lim_buffer(3) = min( rcut_hard, rj + atom_sigma_scaled ) ! upper limit right
!            write(*,'(A)') '* Buffer zone:'
!           contribution to the expansion coeff.
            I0 = transpose( M_radial_poly(lim_buffer, max(7, alpha_max + 2), rcut_hard) )
            call get_constant_poly_filter_coeff(rj, atom_sigma_scaled, rcut_soft, rcut_hard, 'left', B)
            M_left(1:7, 1) = matmul(  B(1:7), M_radial_monomial(lim_buffer(1), 6) ) * I0(1:7, 1)
            M_left(1:7, 2) = matmul( B(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
            I_left = matmul( A(1:alpha_max, 1:7), M_left(1:7, 2) ) * I0(3:alpha_max + 2, 2) - &
                     matmul( A(1:alpha_max, 1:7), M_left(1:7, 1) ) * I0(3:alpha_max + 2, 1)
            call get_constant_poly_filter_coeff(rj, atom_sigma_scaled, rcut_soft, rcut_hard, 'right', B)
            M_right(1:7, 1) = matmul( B(1:7), M_radial_monomial(lim_buffer(2), 6) ) * I0(1:7, 2)
            M_right(1:7, 2) = matmul( B(1:7), M_radial_monomial(lim_buffer(3), 6) ) * I0(1:7, 3)
            I_right = matmul( A(1:alpha_max, 1:7), M_right(1:7, 2) ) * I0(3:alpha_max + 2, 3) - &
                      matmul( A(1:alpha_max, 1:7), M_right(1:7, 1) ) * I0(3:alpha_max + 2, 2)
            exp_coeff_buffer = I_left + I_right
          end if
!         Transform from g_alpha to g_n (the orthonormal basis)
          exp_coeff(1:alpha_max, k) = amplitude * matmul( W, exp_coeff_soft(1:alpha_max) + exp_coeff_buffer(1:alpha_max) )
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
!   *********************************************

!   This is for debugging
    if( .false. )then
      open(10, file="radial_expansion_coefficients.dat", status="unknown", access="append")
      write(10,*) exp_coeff(1:alpha_max, 1)
      close(10)
    end if

    deallocate( exp_coeff_soft, exp_coeff_buffer, A, I0, M_left, M_right, I_left, I_right)

  return
  end subroutine
!**************************************************************************


!**************************************************************************
!
! Auxiliary function that computes a piecewise polynomial function and
! its derivatives.
!
  function g_aux(r, r0, sigma, piece)
    implicit none
 
    real*8, intent(in) :: r, r0, sigma
    character(*), intent(in) :: piece
    real*8 :: x
    real*8 :: g_aux(1:4)

    x = (r - r0)/sigma

    if ( piece == "left" )then
      g_aux = [ 1.d0 - 3.d0*x**2 - 2.d0*x**3, -6.d0*(x**2 + x)/sigma, -3.d0*(2.d0*x + 1)/sigma**2, -2.d0/sigma**3 ]
    else if ( piece == "right" )then
      g_aux = [ 1.d0 - 3.d0*x**2 + 2.d0*x**3, 6.d0*(x**2 - x)/sigma, 3.d0*(2.d0*x - 1)/sigma**2, 2.d0/sigma**3 ]
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
    print *, l
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
  subroutine get_constant_poly_filter_coeff(rj, sigma_j, rcut_soft, rcut_hard, piece, B)
    implicit none

    real*8, intent(in) :: rj, sigma_j, rcut_soft, rcut_hard
    character(*), intent(in) :: piece
    integer :: i
    real*8, allocatable :: C_filter(:), col_poly(:), C_poly(:,:)
    real*8, intent(inout) :: B(:)

!   coeff. from the filter
    allocate( C_filter(1:4) )
    C_filter = g_aux(0.d0, rcut_soft, rcut_hard - rcut_soft, "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    allocate( C_poly(1:7, 1:4) )
    allocate( col_poly(1:4) )
    C_poly = 0.d0
    col_poly(1:4) = g_aux(0.d0, rj, sigma_j, piece)

    do i = 1, 4
      C_poly(i:i+4, i) = col_poly
    end do

    B = matmul( C_poly, C_filter )

    deallocate( C_filter, C_poly, col_poly )

  end subroutine
!**************************************************************************


!**************************************************************************
!
! This subroutine returns the radial terms coming from the
! polynomial radial basis 
! (in the notes: I0, R_hard)
!
  function M_radial_poly(r, alpha_max, rcut)
    implicit none

    integer, intent(in) :: alpha_max
    real*8, intent(in) :: rcut, r(:)
    integer :: i
    real*8, allocatable :: radial_terms(:,:), M_radial_poly(:,:)

    allocate( radial_terms(1:size(r), 1:alpha_max) )

    radial_terms = 1.d0
    do i = 2, alpha_max
      radial_terms(:, i) = radial_terms(:, i - 1) * (1.d0 - r/rcut)
    end do

    M_radial_poly = radial_terms

    deallocate( radial_terms )

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns a vector of monomial terms (1, x, x**2, ...)
! for a given polynomial degree
! (in the notes: R, R_ext) 
!
  function radial_monomial(r, degree)
    implicit none

    integer, intent(in) :: degree
    real*8, intent(in) :: r
    integer :: p
    real*8, allocatable :: radial_terms(:), radial_monomial(:) 
    
    allocate( radial_terms(1:degree + 1) )

    radial_terms = 1.d0
    do p = 2, degree + 1
      radial_terms(p) = r * radial_terms(p - 1) 
    end do

    radial_monomial = radial_terms

    deallocate( radial_terms )

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
  function M_radial_monomial(r, degree)
    implicit none

!    integer, intent(in) :: degree
    real*8, intent(in) :: r
    integer :: i, j, degree 
    real*8, allocatable :: coeff(:,:), radial_terms(:), M_radial_monomial(:,:)
    
    degree = 6
    allocate( coeff(1:degree, 1:degree + 1) )
    allocate( radial_terms(1:degree + 1) )
    allocate( M_radial_monomial(1:degree + 1, 1:degree + 1) )

!    coeff = [(dfloat(j), j=1, degree + 1)]
    coeff = reshape([0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0,    &
                     0.d0, 0.d0, 2.d0, 6.d0, 12.d0, 20.d0, 30.d0, &
                     0.d0, 0.d0, 0.d0, 6.d0, 24.d0, 60.d0, 120.d0, &
                     0.d0, 0.d0, 0.d0, 0.d0, 24.d0, 120.d0, 360.d0, &
                     0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 120.d0, 720.d0, &
                     0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 720.d0], shape(coeff))
    radial_terms = radial_monomial(r, degree)

    M_radial_monomial = 0.d0
    M_radial_monomial(1:degree + 1, 1) = radial_terms

    do i = 2, degree + 1
      M_radial_monomial(i:degree + 1, i) = radial_terms(1:degree + 2 - i) * coeff(i:degree + 1, i - 1)
    end do
    
    deallocate( coeff, radial_terms )

  end function
!**************************************************************************




!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set. It is also needed to transform between original
! basis and orthonormal basis. Unlike the original function, it does not
! require blas/lapack to work and instead relies on pretabulated values.
!
  subroutine get_orthonormalization_matrix_poly3_tabulated_copy(alpha_max, S, W)

    implicit none

    integer :: alpha_max
    real*8, intent(inout) :: W(:,:), S(:,:)

    select case(alpha_max)
      case(:0)
        write(*,*) "Bad value of alpha_max"
        stop
      case(1)
        W = reshape([1.0000000000000000d0], shape(W))
        S = reshape([1.0000000000000000d0], shape(S))
      case(2)
        W = reshape([5.9999999999999920d0, -5.2915026221291726d0, -5.2915026221291726d0,  &
                    5.9999999999999920d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0], shape(S))
      case(3)
        W = reshape([16.532871534059733d0, -29.078310649176263d0, 13.308493852760057d0, -29.078310649176263d0,  &
                    65.256761373659685d0, -36.000096455630107d0, 13.308493852760057d0, -36.000096455630107d0,  &
                    23.492063480194044d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0], shape(S))
      case(4)
        W = reshape([31.619473929163473d0, -82.899603631697133d0, 76.876069386827055d0, -24.858289195076921d0,  &
                    -82.899603631697133d0, 278.76362329353776d0, -309.81360154834925d0, 114.16667789081102d0,  &
                    76.876069386827055d0, -309.81360154834925d0, 401.66898510515369d0, -168.42692378140524d0,  &
                    -24.858289195076921d0, 114.16667789081102d0, -168.42692378140524d0, 79.877446522677999d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0,  &
                    0.98333216603563356d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0,  &
                    1.0000000000000000d0], shape(S))
      case(5)
        W = reshape([48.701840798683406d0, -167.01751270577688d0, 232.29114353203897d0, -152.17733678177117d0,  &
                    38.937949480633065d0, -167.01751270577688d0, 745.79170301284648d0, -1251.5498958156577d0,  &
                    936.75356599200450d0, -263.84749322361466d0, 232.29114353203897d0, -1251.5498958156577d0,  &
                    2447.7952379389576d0, -2072.0455670547381d0, 643.96375585768624d0, -152.17733678177117d0,  &
                    936.75356599200450d0, -2072.0455670547381d0, 1959.7497439691933d0, -672.11823364760107d0,  &
                    38.937949480633065d0, -263.84749322361466d0, 643.96375585768624d0, -672.11823364760107d0,  &
                    253.84463194408610d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.99215674164922152d0, 1.0000000000000000d0,  &
                    0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0,  &
                    0.99744571741206722d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0], shape(S))
      case(6)
        W = reshape([65.344980865540876d0, -270.96483658090318d0, 495.43665447034391d0, -487.28339480337212d0,  &
                    252.44181841585481d0, -54.246284684757562d0, -270.96483658090318d0, 1493.2157334161463d0,  &
                    -3338.9403306650820d0, 3783.7603589943619d0, -2167.7466416209732d0, 500.79512237094178d0,  &
                    495.43665447034391d0, -3338.9403306650820d0, 8760.0319024860873d0, -11249.630736045141d0,  &
                    7104.9320746203475d0, -1771.4461242822615d0, -487.28339480337212d0, 3783.7603589943619d0,  &
                    -11249.630736045141d0, 16097.101050349149d0, -11150.857982565711d0, 3007.1993994548748d0,  &
                    252.44181841585481d0, -2167.7466416209732d0, 7104.9320746203475d0, -11150.857982565711d0,  &
                    8419.3910413187914d0, -2457.9617037326980d0, -54.246284684757562d0, 500.79512237094178d0,  &
                    -1771.4461242822615d0, 3007.1993994548748d0, -2457.9617037326980d0, 776.42840231397884d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0,  &
                    0.95148591360407542d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.95393920141694566d0,  &
                    0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0,  &
                    0.99107124982123374d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.90905934288630952d0,  &
                    0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0,  &
                    1.0000000000000000d0], shape(S))
      case(7)
        W = reshape([80.112087744414481d0, -380.83753196337295d0, 846.68851603934706d0, -1097.4660724104394d0,  &
                    854.00203553098982d0, -371.15373926820087d0, 69.379320671488060d0, -380.83753196337295d0,  &
                    2460.4826695826819d0, -6806.5629280860976d0, 10292.638409627154d0, -8930.9770566136904d0,  &
                    4194.6906070515997d0, -829.33608862620451d0, 846.68851603934706d0, -6806.5629280860976d0,  &
                    22347.744248873638d0, -38562.685975984692d0, 37023.176288616021d0, -18793.587621101400d0,  &
                    3945.6376780219634d0, -1097.4660724104394d0, 10292.638409627154d0, -38562.685975984692d0,  &
                    74347.213817597105d0, -78269.402026367752d0, 42877.575814051277d0, -9587.7268556762156d0,  &
                    854.00203553098982d0, -8930.9770566136904d0, 37023.176288616021d0, -78269.402026367752d0,  &
                    89507.134104252254d0, -52778.316794353603d0, 12594.783311461077d0, -371.15373926820087d0,  &
                    4194.6906070515997d0, -18793.587621101400d0, 42877.575814051277d0, -52778.316794353603d0,  &
                    33381.624954509396d0, -8510.6924230665845d0, 69.379320671488060d0, -829.33608862620451d0,  &
                    3945.6376780219634d0, -9587.7268556762156d0, 12594.783311461077d0, -8510.6924230665845d0,  &
                    2318.7297663139648d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0,  &
                    0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.95393920141694566d0, 0.98333216603563356d0,  &
                    0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0,  &
                    0.98226460284385697d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0,  &
                    0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0], shape(S))
      case(8)
        W = reshape([92.539024757823469d0, -485.62295892374391d0, 1245.3098138485498d0, -1970.1315498890267d0,  &
                    2024.2063495743746d0, -1321.8047757365766d0, 499.34651114023347d0, -83.121761910171273d0,  &
                    -485.62295892374391d0, 3539.9434708348890d0, -11511.708266971060d0, 21569.895229194808d0,  &
                    -24986.945984915284d0, 17767.671463504408d0, -7134.3810386646874d0, 1241.2366889962939d0,  &
                    1245.3098138485498d0, -11511.708266971060d0, 45076.344590889283d0, -97380.089133774367d0,  &
                    125682.09286798064d0, -97005.646137706557d0, 41462.330953470708d0, -7568.2374319869477d0,  &
                    -1970.1315498890267d0, 21569.895229194808d0, -97380.089133774367d0, 236558.01704791249d0,  &
                    -335931.86898086715d0, 280178.31947277795d0, -127521.98045884802d0, 24497.997672243826d0,  &
                    2024.2063495743746d0, -24986.945984915284d0, 125682.09286798064d0, -335931.86898086715d0,  &
                    518562.62266775075d0, -464898.16962265247d0, 225191.00178569887d0, -45642.664557060503d0,  &
                    -1321.8047757365766d0, 17767.671463504408d0, -97005.646137706557d0, 280178.31947277795d0,  &
                    -464898.16962265247d0, 445491.35949089238d0, -229343.07661604194d0, 49131.675673804595d0,  &
                    499.34651114023347d0, -7134.3810386646874d0, 41462.330953470708d0, -127521.98045884802d0,  &
                    225191.00178569887d0, -229343.07661604194d0, 125237.47296520101d0, -28390.563739949452d0,  &
                    -83.121761910171273d0, 1241.2366889962939d0, -7568.2374319869477d0, 24497.997672243826d0,  &
                    -45642.664557060503d0, 49131.675673804595d0, -28390.563739949452d0, 6814.4475643740907d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0,  &
                    0.98333216603563356d0, 0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0,  &
                    0.91651513899116799d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.96378881965339736d0,  &
                    0.94991775959816649d0, 0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0,  &
                    1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0, 0.98226460284385697d0,  &
                    0.97192421422695907d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0,  &
                    0.99107124982123374d0, 0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0,  &
                    0.99444440145743085d0, 0.88712019959006128d0, 0.93404977361585861d0, 0.96378881965339736d0,  &
                    0.98226460284385697d0, 0.99305547153730200d0, 0.99845559753396829d0, 1.0000000000000000d0,  &
                    0.99874921777190884d0, 0.86602540378443860d0, 0.91651513899116799d0, 0.94991775959816649d0,  &
                    0.97192421422695907d0, 0.98601329718326935d0, 0.99444440145743085d0, 0.99874921777190884d0,  &
                    1.0000000000000000d0], shape(S))
      case(9)
        W = reshape([102.73996023131654d0, -579.28008477226797d0, 1649.6329933089637d0, -3023.3784007460226d0,  &
                    3801.6152070879762d0, -3281.1404539213454d0, 1861.8484048267730d0, -625.49270608226868d0,  &
                    94.172577622820413d0, -579.28008477226797d0, 4626.3192102608891d0, -17024.980936140477d0,  &
                    37570.528484501476d0, -53968.531249935149d0, 51220.175635520885d0, -31105.502641416926d0,  &
                    10973.451092982890d0, -1712.0993861885026d0, 1649.6329933089637d0, -17024.980936140477d0,  &
                    76648.552042974829d0, -197491.45685094723d0, 319059.39758025052d0, -330888.18983548041d0,  &
                    214867.87633591256d0, -79759.988337705232d0, 12939.555536288681d0, -3023.3784007460226d0,  &
                    37570.528484501476d0, -197491.45685094723d0, 577420.51966726186d0, -1032806.0301229986d0,  &
                    1161968.8118061803d0, -805283.01678246027d0, 314926.85573368822d0, -53282.713938116933d0,  &
                    3801.6152070879762d0, -53968.531249935149d0, 319059.39758025052d0, -1032806.0301229986d0,  &
                    2015313.2813086673d0, -2441081.2560675377d0, 1800654.0587436187d0, -742277.93893866800d0,  &
                    131305.73650891884d0, -3281.1404539213454d0, 51220.175635520885d0, -330888.18983548041d0,  &
                    1161968.8118061803d0, -2441081.2560675377d0, 3159490.1849562805d0, -2472607.0610428583d0,  &
                    1074278.1958661408d0, -199099.57483824741d0, 1861.8484048267730d0, -31105.502641416926d0,  &
                    214867.87633591256d0, -805283.01678246027d0, 1800654.0587436187d0, -2472607.0610428583d0,  &
                    2045686.8262865648d0, -936138.89101480751d0, 182064.24577301860d0, -625.49270608226868d0,  &
                    10973.451092982890d0, -79759.988337705232d0, 314926.85573368822d0, -742277.93893866800d0,  &
                    1074278.1958661408d0, -936138.89101480751d0, 450714.89061671734d0, -92090.959139143888d0,  &
                    94.172577622820413d0, -1712.0993861885026d0, 12939.555536288681d0, -53282.713938116933d0,  &
                    131305.73650891884d0, -199099.57483824741d0, 182064.24577301860d0, -92090.959139143888d0,  &
                    19782.408284583180d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.99215674164922152d0, 1.0000000000000000d0,  &
                    0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0, 0.95148591360407542d0,  &
                    0.93404977361585861d0, 0.91651513899116799d0, 0.89921841062113494d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.94991775959816649d0, 0.93564551297569798d0,  &
                    0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0,  &
                    0.99744571741206722d0, 0.99107124982123374d0, 0.98226460284385697d0, 0.97192421422695907d0,  &
                    0.96064535921058791d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.90905934288630952d0, 0.95148591360407542d0,  &
                    0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0, 1.0000000000000000d0,  &
                    0.99845559753396829d0, 0.99444440145743085d0, 0.98868599666425949d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0, 0.99874921777190884d0, 0.99545452192223205d0,  &
                    0.86602540378443860d0, 0.91651513899116799d0, 0.94991775959816649d0, 0.97192421422695907d0,  &
                    0.98601329718326935d0, 0.99444440145743085d0, 0.99874921777190884d0, 1.0000000000000000d0,  &
                    0.99896640799254144d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0], shape(S))
      case(10)
        W = reshape([110.85899098970923d0, -656.32724208730247d0, 2002.0561806262745d0, -4030.0430894678625d0,  &
                    5747.1082633813776d0, -5876.5673849153736d0, 4220.3047498562373d0, -2019.2279545723384d0,  &
                    576.83178827809900d0, -74.279463179713289d0, -656.32724208730247d0, 5580.5224445491212d0,  &
                    -22331.612043043937d0, 54995.143701985420d0, -91095.538134775474d0, 104230.69173430053d0,  &
                    -81616.918643421988d0, 41841.558030521599d0, -12662.655255751506d0, 1715.2096007504856d0,  &
                    2002.0561806262745d0, -22331.612043043937d0, 110954.23971315096d0, -323698.56673909881d0,  &
                    611814.99361075275d0, -776207.90529295243d0, 659769.38292835664d0, -361528.21759789030d0,  &
                    115656.71780639346d0, -16430.693998954230d0, -4030.0430894678625d0, 54995.143701985420d0,  &
                    -323698.56673909881d0, 1086072.8812834970d0, -2301554.0722570983d0, 3206429.9814929799d0,  &
                    -2943668.1793391723d0, 1719782.1674985890d0, -580745.53209077963d0, 86416.331611407732d0,  &
                    5747.1082633813776d0, -91095.538134775474d0, 611814.99361075275d0, -2301554.0722570983d0,  &
                    5380578.4815178504d0, -8152091.3190899873d0, 8040206.2774941670d0, -4994941.7070384044d0,  &
                    1778403.4588345746d0, -277067.37753186666d0, -5876.5673849153736d0, 104230.69173430053d0,  &
                    -776207.90529295243d0, 3206429.9814929799d0, -8152091.3190899873d0, 13309104.674707124d0,  &
                    -14025079.427185630d0, 9238846.3749270160d0, -3464500.8524018270d0, 565144.53718240885d0,  &
                    4220.3047498562373d0, -81616.918643421988d0, 659769.38292835664d0, -2943668.1793391723d0,  &
                    8040206.2774941670d0, -14025079.427185630d0, 15706401.821684649d0, -10938310.914757373d0,  &
                    4315299.1100801602d0, -737221.24287861609d0, -2019.2279545723384d0, 41841.558030521599d0,  &
                    -361528.21759789030d0, 1719782.1674985890d0, -4994941.7070384044d0, 9238846.3749270160d0,  &
                    -10938310.914757373d0, 8029013.8383061187d0, -3328353.6668580547d0, 595670.14588087215d0,  &
                    576.83178827809900d0, -12662.655255751506d0, 115656.71780639346d0, -580745.53209077963d0,  &
                    1778403.4588345746d0, -3464500.8524018270d0, 4315299.1100801602d0, -3328353.6668580547d0,  &
                    1447833.8533988285d0, -271507.14094415621d0, -74.279463179713289d0, 1715.2096007504856d0,  &
                    -16430.693998954230d0, 86416.331611407732d0, -277067.37753186666d0, 565144.53718240885d0,  &
                    -737221.24287861609d0, 595670.14588087215d0, -271507.14094415621d0, 53355.279479472047d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.82679728470768454d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0,  &
                    0.95148591360407542d0, 0.93404977361585861d0, 0.91651513899116799d0, 0.89921841062113494d0,  &
                    0.88235294117647056d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.96378881965339736d0,  &
                    0.94991775959816649d0, 0.93564551297569798d0, 0.92128466398761111d0, 0.95393920141694566d0,  &
                    0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0,  &
                    0.99107124982123374d0, 0.98226460284385697d0, 0.97192421422695907d0, 0.96064535921058791d0,  &
                    0.94882928301683922d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.96824583655185426d0, 0.90905934288630952d0,  &
                    0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0,  &
                    1.0000000000000000d0, 0.99845559753396829d0, 0.99444440145743085d0, 0.98868599666425949d0,  &
                    0.98169181562325258d0, 0.88712019959006128d0, 0.93404977361585861d0, 0.96378881965339736d0,  &
                    0.98226460284385697d0, 0.99305547153730200d0, 0.99845559753396829d0, 1.0000000000000000d0,  &
                    0.99874921777190884d0, 0.99545452192223205d0, 0.99065885080469862d0, 0.86602540378443860d0,  &
                    0.91651513899116799d0, 0.94991775959816649d0, 0.97192421422695907d0, 0.98601329718326935d0,  &
                    0.99444440145743085d0, 0.99874921777190884d0, 1.0000000000000000d0, 0.99896640799254144d0,  &
                    0.99621210759909562d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0, 0.99913156735681652d0, 0.82679728470768454d0,  &
                    0.88235294117647056d0, 0.92128466398761111d0, 0.94882928301683922d0, 0.96824583655185426d0,  &
                    0.98169181562325258d0, 0.99065885080469862d0, 0.99621210759909562d0, 0.99913156735681652d0,  &
                    1.0000000000000000d0], shape(S))
      case(11)
        W = reshape([117.21235214772904d0, -714.35806297792158d0, 2245.1643321280740d0, -4610.3933792100061d0,  &
                    6511.4009202362895d0, -6150.8302224825920d0, 3391.5937860140148d0, -436.01742880721537d0,  &
                    -738.92275647372799d0, 484.31715587976799d0, -98.454105559226775d0, -714.35806297792158d0,  &
                    6269.8952732304706d0, -25875.069939432633d0, 65110.259187545853d0, -107732.04536419440d0,  &
                    117390.86907425059d0, -78379.654709546943d0, 23457.309598262604d0, 5607.1440832265725d0,  &
                    -6684.8214165052414d0, 1550.5418761193073d0, 2245.1643321280740d0, -25875.069939432633d0,  &
                    132599.60954916809d0, -396361.31134008802d0, 756826.91023577540d0, -944113.58367922041d0,  &
                    746475.35184183146d0, -329819.32141719793d0, 40631.572575033351d0, 26062.477344374791d0,  &
                    -8671.4089555561513d0, -4610.3933792100061d0, 65110.259187545853d0, -396361.31134008802d0,  &
                    1372463.1491815776d0, -2988502.0270795608d0, 4244631.4027743703d0, -3916425.4633714431d0,  &
                    2232511.7325943527d0, -677804.31608534034d0, 55391.291470506068d0, 13595.784351342649d0,  &
                    6511.4009202362895d0, -107732.04536419440d0, 756826.91023577540d0, -2988502.0270795608d0,  &
                    7374201.5953415921d0, -11897330.128776217d0, 12685063.180279443d0, -8765094.4377386160d0,  &
                    3697769.2054832382d0, -829684.91928406083d0, 67971.552751104493d0, -6150.8302224825920d0,  &
                    117390.86907425059d0, -944113.58367922041d0, 4244631.4027743703d0, -11897330.128776217d0,  &
                    21851270.467823889d0, -26754279.226484302d0, 21633158.682814363d0, -11077651.903835505d0,  &
                    3246049.5859478437d0, -412975.14465735806d0, 3391.5937860140148d0, -78379.654709546943d0,  &
                    746475.35184183146d0, -3916425.4633714431d0, 12685063.180279443d0, -26754279.226484302d0,  &
                    37494848.639549591d0, -34689165.786953673d0, 20377863.358143974d0, -6891923.6598007092d0,  &
                    1022531.8799494690d0, -436.01742880721537d0, 23457.309598262604d0, -329819.32141719793d0,  &
                    2232511.7325943527d0, -8765094.4377386160d0, 21633158.682814363d0, -34689165.786953673d0,  &
                    36151306.746357322d0, -23647318.572119914d0, 8826384.4847903624d0, -1434984.5992588799d0,  &
                    -738.92275647372799d0, 5607.1440832265725d0, 40631.572575033351d0, -677804.31608534034d0,  &
                    3697769.2054832382d0, -11077651.903835505d0, 20377863.358143974d0, -23647318.572119914d0,  &
                    16913806.476479985d0, -6818990.6947852084d0, 1186826.9821223991d0, 484.31715587976799d0,  &
                    -6684.8214165052414d0, 26062.477344374791d0, 55391.291470506068d0, -829684.91928406083d0,  &
                    3246049.5859478437d0, -6891923.6598007092d0, 8826384.4847903624d0, -6818990.6947852084d0,  &
                    2933523.2297482090d0, -540611.16396153928d0, -98.454105559226775d0, 1550.5418761193073d0,  &
                    -8671.4089555561513d0, 13595.784351342649d0, 67971.552751104493d0, -412975.14465735806d0,  &
                    1022531.8799494690d0, -1434984.5992588799d0, 1186826.9821223991d0, -540611.16396153928d0,  &
                    104864.79597292397d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.82679728470768454d0, 0.80868982852161886d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0,  &
                    0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0, 0.91651513899116799d0,  &
                    0.89921841062113494d0, 0.88235294117647056d0, 0.86602540378443871d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.94991775959816649d0, 0.93564551297569798d0,  &
                    0.92128466398761111d0, 0.90703620734810975d0, 0.95393920141694566d0, 0.98333216603563356d0,  &
                    0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0,  &
                    0.98226460284385697d0, 0.97192421422695907d0, 0.96064535921058791d0, 0.94882928301683922d0,  &
                    0.93674969975975964d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.96824583655185426d0, 0.95831484749990992d0,  &
                    0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0,  &
                    0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0, 0.99444440145743085d0,  &
                    0.98868599666425949d0, 0.98169181562325258d0, 0.97383114934675230d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0, 0.99874921777190884d0, 0.99545452192223205d0,  &
                    0.99065885080469862d0, 0.98476101329618471d0, 0.86602540378443860d0, 0.91651513899116799d0,  &
                    0.94991775959816649d0, 0.97192421422695907d0, 0.98601329718326935d0, 0.99444440145743085d0,  &
                    0.99874921777190884d0, 1.0000000000000000d0, 0.99896640799254144d0, 0.99621210759909562d0,  &
                    0.99215674164922152d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0, 0.99913156735681652d0, 0.99679486355016889d0,  &
                    0.82679728470768454d0, 0.88235294117647056d0, 0.92128466398761111d0, 0.94882928301683922d0,  &
                    0.96824583655185426d0, 0.98169181562325258d0, 0.99065885080469862d0, 0.99621210759909562d0,  &
                    0.99913156735681652d0, 1.0000000000000000d0, 0.99926008128973698d0, 0.80868982852161886d0,  &
                    0.86602540378443871d0, 0.90703620734810975d0, 0.93674969975975964d0, 0.95831484749990992d0,  &
                    0.97383114934675230d0, 0.98476101329618471d0, 0.99215674164922152d0, 0.99679486355016889d0,  &
                    0.99926008128973698d0, 1.0000000000000000d0], shape(S))
      case(12:)
        write(*,*) "Bad value of alpha_max"
        stop
    end select

  return
  end subroutine
!**************************************************************************


end module soap_turbo_radial_op

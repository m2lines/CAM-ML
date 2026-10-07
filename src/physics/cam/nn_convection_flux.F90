module nn_convection_flux_mod
    !! Code to perform the Convection Flux parameterisation
    !! and interface to the nn_cf_net Neural Net
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
  ! Libraries to use
  ! XL XL XL
use shr_kind_mod, only: r8=>shr_kind_r8  
use nn_cf_net_mod, only: nn_cf_net_init, net_forward, nn_cf_net_finalize
use SAM_consts_mod, only: fac_cond, fac_fus, tprmin, a_pr, input_ver_dim, &
                          nrf, nrfq, ggr
use cam_logfile,    only: iulog

implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize, &
        esati, rsati, rsatw, dtrsatw, dtrsati


!---------------------------------------------------------------------
! local/private data

! Neural Net parameters used throughout module

integer :: n_inputs, n_outputs
    !! Length of input/output vector to the NN

logical :: do_init=.true.
    !! model initialisation is yet to be performed


!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_init(nn_filename)
        !! Initialise the NN module

        character(len=136), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights

        ! Initialise the Neural Net from file and get info
        call nn_cf_net_init(nn_filename, n_inputs, n_outputs, nrf)

        ! Set init flag as complete
        do_init = .false.

    end subroutine nn_convection_flux_init


    subroutine nn_convection_flux(tabs_i, q_i, ps_i, landfrac, &
                                  t_rad_rest_tend, t_flux_adv, &
                                  q_flux_adv, q_tend_auto, q_sed_flux)
      
      !! Interface to the neural net that applies physical constraints and reshaping
      !! of variables.
      !! Operates on subcycle of timestep dtn to update t, q, precsfc, and prec_xy

      ! -----------------------------------
      ! Input Variables
      ! -----------------------------------
      
      ! ---------------------
      ! Fields from beginning of time step used as NN inputs
      ! ---------------------
      != unit s :: tabs_i
      real(8), intent(in) :: tabs_i(:, :)
      !! Temperature
      
      != unit 1 :: q_i
      real(8), intent(in) :: q_i(:, :)
      !! Non-precipitating water mixing ratio

      != unit Pa :: ps_i
      real(8), intent(in) :: ps_i(:)
      !! Surface pressure

      != unit 1 :: landfrac
      real(8), intent(in) :: landfrac(:)
      !! landfraction
            
      ! -----------------------------------
      ! Local Variables
      ! -----------------------------------
      integer  i, k, dim_counter, out_dim_counter
      integer ncol
      real(8),   dimension(size(tabs_i,1),size(tabs_i, 2)) :: nn_mask      
      ! -----------------------------------
      ! variables for NN
      ! -----------------------------------
      real(4), dimension(n_inputs) :: features
      !! Vector of input features for the NN
      real(4), dimension(n_outputs) :: outputs
      !! vector of output features from the NN
      ! NN outputs
      real(8),  intent(out), dimension(size(tabs_i, 1), size(tabs_i, 2)) :: t_flux_adv, q_flux_adv, q_tend_auto, &
           q_sed_flux, t_rad_rest_tend

      real(8), allocatable, dimension(:) :: landfrac_bin

      ncol = size(tabs_i, 1)
      allocate(landfrac_bin(ncol))      
      
      ! Check that we have initialised all of the variables.
      if (do_init) call error_mesg('NN has not yet been initialised using nn_convection_flux_init.')
      
      ! The NN operates on atmospheric columns which have been flattened into 2D
      do i=1,ncol
         
         ! Initialize variables
         features = 0.
         dim_counter = 0
         outputs = 0.

         t_rad_rest_tend(i,:) = 0._r8
         t_flux_adv(i,:) = 0._r8
         q_flux_adv(i,:) = 0._r8
         q_tend_auto(i,:) = 0._r8
         q_sed_flux(i,:) = 0._r8
         landfrac_bin(i) = 0._r8
         
         ! transform landfrac from continuous to binary
         if (landfrac(i) .gt. 0.5_r8) then
            landfrac_bin(i)=1._r8
         else
            landfrac_bin(i)=0._r8
         end if
         
         !-----------------------------------------------------
         ! Combine all features into one vector
         
         ! Add temperature as input feature
         features(dim_counter+1:dim_counter + input_ver_dim) = real(tabs_i(i,1:input_ver_dim),4)
         dim_counter = dim_counter + input_ver_dim

         ! Add non-precipitating water mixing ratio
         features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,1:input_ver_dim),4)
         dim_counter =  dim_counter + input_ver_dim

         ! Add land fraction
         features(dim_counter+1:dim_counter+1) = real(landfrac_bin(i),4)
         dim_counter =  dim_counter + 1

         ! Add surface pressure
         features(dim_counter+1:dim_counter+1) = real(ps_i(i)/100._r8,4)
         dim_counter =  dim_counter + 1         
           
         !-----------------------------------------------------
         ! Call the forward method of the NN on the input features
         ! Scaling and normalisation done as layers in NN

         call net_forward(features, outputs)
         !-----------------------------------------------------
         ! Separate physical outputs from NN output vector
           
         out_dim_counter = 0
           
         ! Moist Static Energy radiative + rest(microphysical changes) tendency
         t_rad_rest_tend(i,1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
         out_dim_counter = out_dim_counter + nrf

         ! Moist Static Energy advective subgrid flux
         ! BC: advection surface flux is zero
         t_flux_adv(i,1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
         t_flux_adv(i,1) = 0._r8         
         out_dim_counter = out_dim_counter + nrf                     

         ! Total non-precip. water mix. ratio advective flux
         ! BC: advection surface flux is zero
         q_flux_adv(i,1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
         q_flux_adv(i,1) = 0._r8         
         out_dim_counter = out_dim_counter + nrf         

         ! Total non-precip. water autoconversion tendency
         q_tend_auto(i,1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
         out_dim_counter = out_dim_counter + nrf

         ! total non-precip. water mix. ratio ice-sedimenting flux
         q_sed_flux(i,1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
                  
      end do

      deallocate(landfrac_bin)
      
    end subroutine nn_convection_flux


    subroutine nn_convection_flux_finalize()
        !! Finalize the NN module

        ! Finalize the Neural Net deallocating arrays
        call nn_cf_net_finalize()

    end subroutine nn_convection_flux_finalize


    !-----------------------------------------------------------------

    subroutine error_mesg (message)
      character(len=*), intent(in) :: message
          !! message to be written to output   (character string)

      ! Since masterproc is a SAM variable adjust to print from all proc for now
      ! if(masterproc) print*, 'Neural network  module: ', message
      print*, 'Neural network  module: ', message
      stop

    end subroutine error_mesg


    !-----------------------------------------------------------------
    ! Need rsatw functions to:
    !     - run with rf_uses_rh option (currently unused)
    !     - convert variables in CAM interface
    ! Ripped from SAM model:
    ! https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/sat.f90
    ! Note r is used instead of q as q in CAM signifies a moist mixing ratio whilst SAM
    ! uses dry mixing ratios.
    !
    ! Saturation vapor pressure and mixing ratio.
    ! Based on Flatau et.al, (JAM, 1992:1507)

    != unit mb :: esatw
    real(8) function esatw(t)
      != unit K :: t
      real(8), intent(in) :: t  ! temperature (K)

      != unit :: a0
      != unit :: mb / k :: a1, a2, a3, a4, a5, a6, a7, a8
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
              6.105851, 0.4440316, 0.1430341e-1, &
              0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
              0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
      !       6.11239921, 0.443987641, 0.142986287e-1, &
      !       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
      !       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/

      != unit K :: dt
      real(8) :: dt
      dt = max(-80._r8,t-273.16_r8)
      esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function esatw


    != unit 1 :: rsatw
    real(8) function rsatw(t,p)
      != unit K :: t
      real(8), intent(in) :: t  ! temperature

      != unit mb :: p, esat
      real(8) :: p  ! pressure
      real(8) :: esat

      esat = esatw(t)
      rsatw = 0.622_r8 * esat/max(esat, p-esat)
    end function rsatw


    real(8) function dtesatw(t)
      real(8), intent(in) :: t  ! temperature (K)
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
                0.443956472, 0.285976452e-1, 0.794747212e-3, &
                0.121167162e-4, 0.103167413e-6, 0.385208005e-9, &
               -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
      real(8) :: dt
      dt = max(-80._r8,t-273.16_r8)
      dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function dtesatw


    real(8) function dtrsatw(t,p)
      real(8), intent(in) :: t  ! temperature (K)
      real(8), intent(in) :: p  ! pressure    (mb)
      dtrsatw=0.622_r8*dtesatw(t)/p
    end function dtrsatw


    real(8) function esati(t)
      != unit K :: t
      real(8), intent(in) :: t  ! temperature
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
              6.11147274, 0.503160820, 0.188439774e-1, &
              0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
              0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
      real(8) :: dt
      dt = max(-80._r8, t-273.16_r8)
      esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function esati


    != unit 1 :: rsati
    real(8) function rsati(t,p)
      != unit t :: K
      real(8), intent(in) :: t  ! temperature

      != unit mb :: p
      real(8), intent(in) :: p  ! pressure

      != unit mb :: esat
      real(8) :: esat
      esat = esati(t)
      rsati = 0.622_r8 * esat/max(esat,p-esat)
    end function rsati


    real(8) function dtesati(t)
      real(8), intent(in) :: t  ! temperature (K)
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 / &
              0.503223089, 0.377174432e-1,0.126710138e-2, &
              0.249065913e-4, 0.312668753e-6, 0.255653718e-8, &
              0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
      real(8) :: dt
      !dt = max(-800. ,t-273.16)
      dt = max(-80._r8 ,t-273.16_r8)      ! XL XL XL
      dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
      return
    end function dtesati


    real(8) function dtrsati(t,p)
      real(8), intent(in) :: t  ! temperature (K)
      real(8), intent(in) :: p  ! pressure    (mb)
      dtrsati = 0.622_r8 * dtesati(t) / p
    end function dtrsati


end module nn_convection_flux_mod

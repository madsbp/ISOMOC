program isopycnalMOC
    use netcdf
    use gsw_mod_kinds
    use gsw_mod_netcdf
    use gsw_mod_toolbox
    use gsw_mod_error_functions, only : gsw_error_code, gsw_error_limit    

    implicit none
      
    !Declarations for calculations
    integer, parameter :: nx = 3600, ny = 800, nz = 62, nsigma = 200, nt = 324 !Size of output file
    real (kind = 8) :: dsigma,sigmas,sigmae,c0
    real (kind = 8), allocatable, dimension(:) :: sigma, time, LAT, maxPD, surfavPD
    real (kind = 8), allocatable, dimension(:,:) :: MOCtmres, MOCtimemean, MOCzmtm, zmVVELmean,&
    zmPDmean, zsfu, area, MOCstanding, MOCtransient, vmPDmean, vsfu
    real (kind = 8), allocatable, dimension(:,:,:) :: fu, ft, ft_nopartial, w, PDint, MOCres, PDmean, VVELmean
    integer :: s,i,j,k,n

    !Declarations for grid input file
    character (len = *), parameter :: input_gridfile = "/home/nsz465/erda_ro/Output/ctrl.g.e11.G.T62_t12.002.pop.h.grid.nc"
    integer :: ncid_grid_in

    character (len = 2), parameter :: HT_name_in = "HT", HU_name_in = "HU"
    real (kind = 8), allocatable, dimension(:,:) :: HT,HU
    integer :: varid_HT, varid_HU

    character (len = 2), parameter :: dz_name_in = "dz"
    real (kind = 8), allocatable, dimension(:) :: dz 
    integer :: varid_dz

    character (len = 3), parameter :: DXU_name_in = "DXU"
    real (kind = 8), allocatable, dimension(:,:) :: DXU 
    integer :: varid_DXU

    character (len = 5), parameter :: TAREA_name_in = "TAREA"
    real (kind = 8), allocatable, dimension(:,:) :: TAREA 
    integer :: varid_TAREA

    character (len = 4), parameter :: ULAT_name_in = "ULAT"
    real (kind = 8), allocatable, dimension(:,:) :: ULAT 
    integer :: varid_ULAT

    character (len = 3), parameter :: z_t_name_in = "z_t"
    real (kind = 8), allocatable, dimension(:) :: z_t 
    integer :: varid_z_t

    character (len = 7), parameter :: z_w_bot_name_in = "z_w_bot"
    real (kind = 8), allocatable, dimension(:) :: z_w_bot 
    integer :: varid_z_w_bot

    character (len = 7), parameter :: z_w_top_name_in = "z_w_top"
    real (kind = 8), allocatable, dimension(:) :: z_w_top 
    integer :: varid_z_w_top

    !Declarations for model input file
    character (len = 130) :: input_file
    integer :: ncid_in, month, year

    real (kind = 8), allocatable, dimension(:,:,:,:) :: PD

    character (len = 4) :: SALT_name_in = "SALT", TEMP_name_in = "TEMP"
    real (kind = 8), allocatable, dimension(:,:,:,:) :: SALT, TEMP 
    integer :: varid_SALT, varid_TEMP

    character (len = 4) :: VVEL_name_in = "VVEL"
    real (kind = 8), allocatable, dimension(:,:,:,:) :: VVEL 
    integer :: varid_VVEL
    
    character (len = 4) :: t_name_in = "time"
    real (kind = 8), allocatable :: t
    integer :: varid_t

    !Declarations for output file
    character (len = *), parameter :: output_filename = "/cache/io.erda.dk+home-nsz465-erda+nsz465/&
    ctrl.g.e11.G.T62_t12.002.pop.h.isopycnalMOC.0016-0042.200siglevs.nc"
    integer :: ncid_out, dimid_lat, dimid_time, dimid_sigma !Dimension and file ID
    integer :: varid_lat, varid_sigma, varid_time !Coordinate variable ID
    integer :: varid_MOCres, varid_MOCtmres, varid_MOCzmtm, varid_MOCtransient,&
    varid_MOCstanding, varid_MOCtimemean, varid_maxPD, varid_surfavPD !Variable ID

    character (len = *), parameter :: lat_dim_name_out = "nlat"
    character (len = *), parameter :: time_dim_name_out = "time"
    character (len = *), parameter :: sigma_dim_name_out = "sigma"
    character (len = *), parameter :: units = "units"
    character (len = *), parameter :: longname = "long_name"

    character (len = *), parameter :: lat_name_out = "LAT"
    character (len = *), parameter :: units_lat = "Degrees north"
    character (len = *), parameter :: longname_lat = "array of t-grid latitudes"

    character (len = *), parameter :: sigma_name_out = "sigma"
    character (len = *), parameter :: units_sigma = "g/cm^3"
    character (len = *), parameter :: longname_sigma = "array of potential density levels "

    character (len = *), parameter :: time_name_out = "time"
    character (len = *), parameter :: units_time = "days since 0000-01-01 00:00:00"
    character (len = *), parameter :: longname_time = "time"

    character (len = *), parameter :: maxPD_name_out = "maxPD"
    character (len = *), parameter :: units_maxPD = "g/cm^3"
    character (len = *), parameter :: longname_maxPD = "Zonal maximum of time-mean PD field for the upper ocean column"

    character (len = *), parameter :: surfavPD_name_out = "surfavPD"
    character (len = *), parameter :: units_surfavPD = "g/cm^3"
    character (len = *), parameter :: longname_surfavPD = "Zonal mean of time-mean surface PD field"

    character (len = *), parameter :: MOCres_name_out = "MOCres"
    character (len = *), parameter :: units_MOCres = "Sverdrups"
    character (len = *), parameter :: longname_MOCres = "Isopycnal MOC"

    character (len = *), parameter :: MOCtmres_name_out = "MOCtmres"
    character (len = *), parameter :: units_MOCtmres = "Sverdrups"
    character (len = *), parameter :: longname_MOCtmres = "Time-mean isopycnal MOC"

    character (len = *), parameter :: MOCzmtm_name_out = "MOCzmtm"
    character (len = *), parameter :: units_MOCzmtm = "Sverdrups"
    character (len = *), parameter :: longname_MOCzmtm = "Isopycnal MOC derived from the zonal- and time-mean fields"

    character (len = *), parameter :: MOCtransient_name_out = "MOCtransient"
    character (len = *), parameter :: units_MOCtransient = "Sverdrups"
    character (len = *), parameter :: longname_MOCtransient = "Transient eddy component"

    character (len = *), parameter :: MOCstanding_name_out = "MOCstanding"
    character (len = *), parameter :: units_MOCstanding = "Sverdrups"
    character (len = *), parameter :: longname_MOCstanding = "Standing eddy component"

    character (len = *), parameter :: MOCtimemean_name_out = "MOCtimemean"
    character (len = *), parameter :: units_MOCtimemean = "Sverdrups"
    character (len = *), parameter :: longname_MOCtimemean = "Isopycnal MOC derived from the time-mean fields"

    !Declaration ended

    !Read from input grid file
    call check( nf90_open(input_gridfile, NF90_NOWRITE, ncid_grid_in) )

    allocate( HT(nx,ny + 1) )
    call check( nf90_inq_varid(ncid_grid_in, HT_name_in, varid_HT) )
    call check( nf90_get_var(ncid_grid_in, varid_HT, HT, (/1,1/), (/nx, ny + 1/)) )

    allocate( HU(nx,ny) )
    call check( nf90_inq_varid(ncid_grid_in, HU_name_in, varid_HU) )
    call check( nf90_get_var(ncid_grid_in, varid_HU, HU, (/1,1/), (/nx, ny/)) )

    allocate( TAREA(nx,ny + 1) )
    call check( nf90_inq_varid(ncid_grid_in, TAREA_name_in, varid_TAREA) )
    call check( nf90_get_var(ncid_grid_in, varid_TAREA, TAREA, (/1,1/), (/nx, ny + 1/)) )

    allocate( DXU(nx,ny) )
    call check( nf90_inq_varid(ncid_grid_in, DXU_name_in, varid_DXU) )
    call check( nf90_get_var(ncid_grid_in, varid_DXU, DXU, (/1,1/), (/nx, ny/)) )

    allocate( ULAT(nx,ny) )
    call check( nf90_inq_varid(ncid_grid_in, ULAT_name_in, varid_ULAT) )
    call check( nf90_get_var(ncid_grid_in, varid_ULAT, ULAT, (/1,1/), (/nx, ny/)) )

    allocate( z_t(nz) )
    call check( nf90_inq_varid(ncid_grid_in, z_t_name_in, varid_z_t) )
    call check( nf90_get_var(ncid_grid_in, varid_z_t, z_t) )

    allocate( dz(nz) )
    call check( nf90_inq_varid(ncid_grid_in, dz_name_in, varid_dz) )
    call check( nf90_get_var(ncid_grid_in, varid_dz, dz) )

    allocate( z_w_bot(nz) )
    call check( nf90_inq_varid(ncid_grid_in, z_w_bot_name_in, varid_z_w_bot) )
    call check( nf90_get_var(ncid_grid_in, varid_z_w_bot, z_w_bot) )

    allocate( z_w_top(nz) )
    call check( nf90_inq_varid(ncid_grid_in, z_w_top_name_in, varid_z_w_top) )
    call check( nf90_get_var(ncid_grid_in, varid_z_w_top, z_w_top) )

    call check( nf90_close(ncid_grid_in) )

    
    
    !!!Calculations!!!
    !Define sigma levels
    allocate( sigma(nsigma) );sigma=0
    sigmas = 1.0330_r8
    sigmae = 1.0375_r8 
    dsigma = (sigmae-sigmas)/(nsigma-1.)
    do s=1,nsigma
        sigma(s) = sigmas + dsigma * (s-1)
    enddo
    print "(a,f12.6)", "Sigma start = ", sigmas
    print "(a,f12.6)", "Sigma end = ", sigmae
    print *, "Sigma levels = ", nsigma
    
    !Calculate latitude axis
    allocate( LAT(ny) )
    do j=1,ny
        do i=1,nx
            if (ULAT(i,j) /= -1.0_r8) then
                LAT(j) = ULAT(i,j)
            end if
        enddo
    enddo
    print "(a)", "Latitude axis has been computed"

    !Compute mask on t and u grid
    allocate( fu(nx,ny,nz) )
    do k=1,nz
        do j=1,ny
            do i=1,nx
                fu(i,j,k) = (HU(i,j) - z_w_top(k))/(z_w_bot(k) - z_w_top(k))
                if (fu(i,j,k) > 1.0_r8) then
                    fu(i,j,k) = 1.0_r8
                else if (fu(i,j,k) < 0.01_r8) then
                    fu(i,j,k) = 0.0_r8
                end if
            enddo
        enddo
    enddo
    print "(a)", "The mask on the u-grid has been computed"

    allocate( ft(nx,ny + 1,nz) )
    allocate( ft_nopartial(nx,ny + 1,nz) )
    do k=1,nz
        do j=1,ny + 1
            do i=1,nx
                ft(i,j,k) = (HT(i,j) - z_w_top(k))/(z_w_bot(k) - z_w_top(k))
                ft_nopartial(i,j,k) = (HT(i,j) - z_w_top(k))/(z_w_bot(k) - z_w_top(k))
                if (ft(i,j,k) > 1.0_r8) then
                    ft(i,j,k) = 1.0_r8
                else if (ft(i,j,k) < 0.01_r8) then
                    ft(i,j,k) = 0.0_r8
                end if
                if (ft_nopartial(i,j,k) >= 0.01_r8) then
                    ft_nopartial(i,j,k) = 1.0_r8
                else if (ft_nopartial(i,j,k) < 0.01_r8) then
                    ft_nopartial(i,j,k) = 0.0_r8
                end if
            enddo
        enddo
    enddo
    print "(a)", "The masks on the t-grid have been computed"

    !Weight matrix for interpolation
    allocate( w(nx,ny,nz) )
    do k=1,nz
        do j=1,ny
            do i=1,nx
                if (i==nx) then
                    w(i,j,k) = ft(i,j,k) * TAREA(i,j) + &
                    ft(1,j,k) * TAREA(1,j) + &
                    ft(i,j + 1,k) * TAREA(i,j + 1) + &
                    ft(1,j + 1,k) * TAREA(1,j + 1)
                else
                    w(i,j,k) = ft(i,j,k) * TAREA(i,j) + &
                    ft(i + 1,j,k) * TAREA(i + 1,j) + &
                    ft(i,j + 1,k) * TAREA(i,j + 1) + &
                    ft(i + 1,j + 1,k) * TAREA(i + 1, j + 1)
                end if
            enddo
        enddo
    enddo
    print "(a)", "The interpolation weights have been computed"
    
    !Interpolation and calculation of isopycnal MOC
    allocate( PDint(nx,ny,nz) )
    allocate( MOCres(ny,nsigma,nt) ); MOCres=0.0_r8
    allocate( PD(nx,ny + 1,nz,1) ); PD=0.0_r8
    allocate( VVEL(nx,ny,nz,1) )
    allocate( SALT(nx,ny + 1,nz,1) )
    allocate( TEMP(nx,ny + 1,nz,1) )
    allocate( time(nt) )
    allocate( t )
    allocate( PDmean(nx,ny + 1,nz) ); PDmean=0.0_r8
    allocate( VVELmean(nx,ny,nz) ); VVELmean=0.0_r8
    month = 1
    year = 16 
    c0 = 1e3_r8 !Scale factor
    do n=1,nt 
        write(input_file, "(a101,i0.4,a1,i0.2,a15)") "/home/nsz465/erda_ro/Ocean/GCS_2015_2017/&
        PostProc_Output/monthly_mean/ctrl.g.e11.G.T62_t12.002.pop.h.", year, "-", month, ".monthlymean.nc"
        print "(a,a)","The current file is: ", input_file
        if (year==25 .and. month==1) then !Loop to change scale factor. At 0025-01, the files are means of 3-day mean fields and the scale factor disappears.
            c0 = 1e0_r8
        end if
        month = month + 1
        if (month==13) then
            month = 1
            year = year + 1
        end if
        call check( nf90_open(input_file, NF90_NOWRITE, ncid_in) )
        call check( nf90_inq_varid(ncid_in, VVEL_name_in, varid_VVEL) )
        call check( nf90_get_var(ncid_in, varid_VVEL, VVEL, (/1,1,1,1/),(/nx,ny,nz,1/)) )
        call check( nf90_inq_varid(ncid_in, TEMP_name_in, varid_TEMP) )
        call check( nf90_get_var(ncid_in, varid_TEMP, TEMP, (/1,1,1,1/),(/nx,ny + 1,nz,1/)) )
        call check( nf90_inq_varid(ncid_in, SALT_name_in, varid_SALT) )
        call check( nf90_get_var(ncid_in, varid_SALT, SALT, (/1,1,1,1/),(/nx,ny + 1,nz,1/)) )
        call check( nf90_inq_varid(ncid_in, t_name_in, varid_t) )
        call check( nf90_get_var(ncid_in, varid_t, t) )
        time(n) = t
        !Compute sigma2000
        do k=1,nz
            do j=1,ny + 1
                do i=1,nx
                    PD(i,j,k,1) = gsw_pd2_from_pt0(ft_nopartial(i,j,k) * SALT(i,j,k,1), &
                    ft_nopartial(i,j,k) * TEMP(i,j,k,1), c0)
                    PDmean(i,j,k) = PDmean(i,j,k) + PD(i,j,k,1)
                enddo
            enddo
        enddo
        !Interpolate and calculate isopycnal MOC
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    VVELmean(i,j,k) = VVELmean(i,j,k) + VVEL(i,j,k,1)
                    if (i==nx) then
                        PDint(i,j,k) = ft(i,j,k) * TAREA(i,j) * PD(i,j,k,1) + &
                        ft(1,j,k) * TAREA(1,j) * PD(1,j,k,1) + &
                        ft(i,j + 1,k) * TAREA(i,j + 1) * PD(i,j + 1,k,1) + &
                        ft(1,j + 1,k) * TAREA(1,j + 1) * PD(1,j + 1,k,1)
                    else
                        PDint(i,j,k) = ft(i,j,k) * TAREA(i,j) * PD(i,j,k,1) + &
                        ft(i + 1,j,k) * TAREA(i + 1,j) * PD(i + 1,j,k,1) + &
                        ft(i,j + 1,k) * TAREA(i,j + 1) * PD(i,j + 1,k,1) + &
                        ft(i + 1,j + 1,k) * TAREA(i + 1,j + 1) * PD(i + 1, j + 1,k,1)
                    end if
                    if (w(i,j,k)==0.0_r8) then
                        PDint(i,j,k) = 0.0_r8
                    else
                        PDint(i,j,k) = PDint(i,j,k) / w(i,j,k)
                    end if
                    do s=1,nsigma
                        if (PDint(i,j,k) > sigma(s)) then
                        MOCres(j,s,n) = MOCres(j,s,n) - fu(i,j,k) * DXU(i,j) * dz(k) * VVEL(i,j,k,1)
                        end if
                    enddo
                enddo
            enddo
        enddo
        call check( nf90_close(ncid_in) )
    enddo 
    print "(a)", "The isopycnal MOC has been calculated"
    
    !Calculate mean isopycnal MOC 
    allocate( MOCtmres(ny,nsigma) )
    MOCtmres = sum(MOCres, 3) / nt
    print "(a)", "The time-mean isopycnal MOC has been calculated"

    !Calculate the isopycnal MOC from the timemean fields
    allocate( MOCtimemean(ny,nsigma) ); MOCtimemean=0.0_r8
    PDmean = PDmean / nt
    VVELmean = VVELmean / nt
    do k=1,nz
        do j=1,ny
            do i=1,nx
                if (i==nx) then
                    PDint(i,j,k) = ft(i,j,k) * TAREA(i,j) * PDmean(i,j,k) + &
                    ft(1,j,k) * TAREA(1,j) * PDmean(1,j,k) + &
                    ft(i,j + 1,k) * TAREA(i,j + 1) * PDmean(i,j + 1,k) + &
                    ft(1,j + 1,k) * TAREA(1,j + 1) * PDmean(1,j + 1,k)
                else
                    PDint(i,j,k) = ft(i,j,k) * TAREA(i,j) * PDmean(i,j,k) + &
                    ft(i + 1,j,k) * TAREA(i + 1,j) * PDmean(i + 1,j,k) + &
                    ft(i,j + 1,k) * TAREA(i,j + 1) * PDmean(i,j + 1,k) + &
                    ft(i + 1,j + 1,k) * TAREA(i + 1,j + 1) * PDmean(i + 1, j + 1,k)
                end if
                if (w(i,j,k)==0.0_r8) then
                    PDint(i,j,k) = 0.0_r8
                else
                    PDint(i,j,k) = PDint(i,j,k) / w(i,j,k)
                end if
                do s=1,nsigma
                    if (PDint(i,j,k) > sigma(s)) then
                    MOCtimemean(j,s) = MOCtimemean(j,s) - fu(i,j,k) * DXU(i,j) * dz(k) * VVELmean(i,j,k)
                    end if
                enddo
            enddo
        enddo
    enddo
    print "(a)", "The isopycnal MOC derived from the time-mean fields has been calculated" 
    
    !Calculate zonal- and time-mean VVEL and PDint 
    allocate( zmVVELmean(ny,nz) ); zmVVELmean=0.0_r8
    allocate( zmPDmean(ny,nz) ); zmPDmean=0.0_r8
    allocate( zsfu(ny,nz) ); zsfu=0.0_r8
    allocate( area(ny,nz) ); area=0.0_r8
    do k=1,nz
        do j=1,ny
            do i=1,nx
                zmVVELmean(j,k) = zmVVELmean(j,k) + fu(i,j,k) * VVELmean(i,j,k)
                zmPDmean(j,k) = zmPDmean(j,k) + fu(i,j,k) * PDint(i,j,k)
                zsfu(j,k) = zsfu(j,k) + fu(i,j,k)
                area(j,k) = area(j,k) + fu(i,j,k) * DXU(i,j) * dz(k)
            enddo
            if (zsfu(j,k) == 0.0_r8) then
                zmVVELmean(j,k) = 0.0_r8
                zmPDmean(j,k) = 0.0_r8
            else 
                zmVVELmean(j,k) = zmVVELmean(j,k) / zsfu(j,k)
                zmPDmean(j,k) = zmPDmean(j,k) / zsfu(j,k)
            end if
        enddo
    enddo

    !Calculate time- and zonal mean isopycnal MOC
    allocate( MOCzmtm(ny,nsigma) ); MOCzmtm=0.0_r8
    do k=1,nz
        do j=1,ny
            do s=1,nsigma
                if (zmPDmean(j,k) > sigma(s)) then
                MOCzmtm(j,s) = MOCzmtm(j,s) - zmVVELmean(j,k) * area(j,k)
                end if
            enddo
        enddo
    enddo
    print "(a)", "The time- and zonal mean isopycnal MOC has been calculated"

    !Calculate standing MOC
    allocate( MOCstanding(ny,nsigma) ); MOCstanding=0.0_r8
    MOCstanding = MOCtimemean - MOCzmtm
    print "(a)", "The standing meander MOC has been calculated"

    !Calculate transient MOC
    allocate( MOCtransient(ny,nsigma) ); MOCtransient=0.0_r8
    MOCtransient = MOCtmres - MOCzmtm - MOCstanding
    print "(a)", "The transient eddy MOC has been calculated"

    !Calculate vert. mean of upper 50m of time mean PD field
    allocate( vmPDmean(nx,ny) ); vmPDmean=0.0_r8
    allocate( vsfu(nx,ny) ); vsfu=0.0_r8
    allocate( maxPD(ny) )
    allocate( surfavPD(ny) )
    do k=1,5
        do j=1,ny
            do i=1,nx
                vmPDmean(i,j) = vmPDmean(i,j) + ft(i,j,k) * PDmean(i,j,k)
                vsfu(i,j) = vsfu(i,j) + ft(i,j,k)
            enddo
        enddo 
    enddo
    do j=1,ny
        do i=1,nx
            if (vsfu(i,j) == 0.0_r8) then
                vmPDmean(i,j) = 0.0_r8
            else
                vmPDmean(i,j) = vmPDmean(i,j) / vsfu(i,j)
            end if
        enddo
    enddo
    maxPD = maxval(vmPDmean, 1)
    surfavPD(:) = zmPDmean(:,1)

    !!!Create output file!!!
    !Meta data
    call check( nf90_create(output_filename, NF90_64BIT_OFFSET, ncid_out) )
    call check( nf90_def_dim(ncid_out, lat_dim_name_out, ny, dimid_lat) )
    call check( nf90_def_dim(ncid_out, sigma_dim_name_out, nsigma, dimid_sigma) )
    call check( nf90_def_dim(ncid_out, time_dim_name_out, nt, dimid_time) )
    !Latitude
    call check( nf90_def_var(ncid_out, lat_name_out, NF90_REAL, (/dimid_lat/), varid_lat) )
    call check( nf90_put_att(ncid_out, varid_lat, units, units_lat) )
    call check( nf90_put_att(ncid_out, varid_lat, longname, longname_lat) )
    !Sigma
    call check( nf90_def_var(ncid_out, sigma_name_out, NF90_REAL, (/dimid_sigma/), varid_sigma) )
    call check( nf90_put_att(ncid_out, varid_sigma, units, units_sigma) )
    call check( nf90_put_att(ncid_out, varid_sigma, longname, longname_sigma) )
    !time
    call check( nf90_def_var(ncid_out, time_name_out, NF90_REAL, (/dimid_time/), varid_time) )
    call check( nf90_put_att(ncid_out, varid_time, units, units_time) )
    call check( nf90_put_att(ncid_out, varid_time, longname, longname_time) )
    !maxPD
    call check( nf90_def_var(ncid_out, maxPD_name_out, NF90_REAL, (/dimid_lat/), varid_maxPD) )
    call check( nf90_put_att(ncid_out, varid_maxPD, units, units_maxPD) )
    call check( nf90_put_att(ncid_out, varid_maxPD, longname, longname_maxPD) )
    !surfavPD
    call check( nf90_def_var(ncid_out, surfavPD_name_out, NF90_REAL, (/dimid_lat/), varid_surfavPD) )
    call check( nf90_put_att(ncid_out, varid_surfavPD, units, units_surfavPD) )
    call check( nf90_put_att(ncid_out, varid_surfavPD, longname, longname_surfavPD) )
    !Residual overturning circulation
    call check( nf90_def_var(ncid_out, MOCres_name_out, NF90_REAL, (/dimid_lat,dimid_sigma,dimid_time/), varid_MOCres) )
    call check( nf90_put_att(ncid_out, varid_MOCres, units, units_MOCres) )
    call check( nf90_put_att(ncid_out, varid_MOCres, longname, longname_MOCres) )
    !Time-mean residual overturning circulation
    call check( nf90_def_var(ncid_out, MOCtmres_name_out, NF90_REAL, (/dimid_lat,dimid_sigma/), varid_MOCtmres) )
    call check( nf90_put_att(ncid_out, varid_MOCtmres, units, units_MOCtmres) )
    call check( nf90_put_att(ncid_out, varid_MOCtmres, longname, longname_MOCtmres) )
    !Zonal- and time-mean overturning circulation
    call check( nf90_def_var(ncid_out, MOCzmtm_name_out, NF90_REAL, (/dimid_lat,dimid_sigma/), varid_MOCzmtm) )
    call check( nf90_put_att(ncid_out, varid_MOCzmtm, units, units_MOCzmtm) )
    call check( nf90_put_att(ncid_out, varid_MOCzmtm, longname, longname_MOCzmtm) )
    !Transient eddy component of overturning circulation
    call check( nf90_def_var(ncid_out, MOCtransient_name_out, NF90_REAL, (/dimid_lat,dimid_sigma/), varid_MOCtransient) )
    call check( nf90_put_att(ncid_out, varid_MOCtransient, units, units_MOCtransient) )
    call check( nf90_put_att(ncid_out, varid_MOCtransient, longname, longname_MOCtransient) )
    !Standing eddy component of overturning circulation
    call check( nf90_def_var(ncid_out, MOCstanding_name_out, NF90_REAL, (/dimid_lat,dimid_sigma/), varid_MOCstanding) )
    call check( nf90_put_att(ncid_out, varid_MOCstanding, units, units_MOCstanding) )
    call check( nf90_put_att(ncid_out, varid_MOCstanding, longname, longname_MOCstanding) )
    !Residual overturning circulation
    call check( nf90_def_var(ncid_out, MOCtimemean_name_out, NF90_REAL, (/dimid_lat,dimid_sigma/), varid_MOCtimemean) )
    call check( nf90_put_att(ncid_out, varid_MOCtimemean, units, units_MOCtimemean) )
    call check( nf90_put_att(ncid_out, varid_MOCtimemean, longname, longname_MOCtimemean) )


    call check( nf90_enddef(ncid_out) ) !This marks the end of the definition of the .nc file
    !Write output to file
    call check( nf90_put_var(ncid_out, varid_sigma, sigma) ) !Write sigma array to output file
    call check( nf90_put_var(ncid_out, varid_time, time) )
    call check( nf90_put_var(ncid_out, varid_lat, LAT) )
    call check( nf90_put_var(ncid_out, varid_maxPD, maxPD) )
    call check( nf90_put_var(ncid_out, varid_surfavPD, surfavPD) )
    call check( nf90_put_var(ncid_out, varid_MOCres, MOCres) ) 
    call check( nf90_put_var(ncid_out, varid_MOCtimemean, MOCtimemean) )
    call check( nf90_put_var(ncid_out, varid_MOCtmres, MOCtmres) )
    call check( nf90_put_var(ncid_out, varid_MOCzmtm, MOCzmtm) )
    call check( nf90_put_var(ncid_out, varid_MOCstanding, MOCstanding) )
    call check( nf90_put_var(ncid_out, varid_MOCtransient, MOCtransient) )
    call check( nf90_close(ncid_out) )
    
    contains
    subroutine check(status)
        integer, intent ( in) :: status
        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check
    function gsw_pd2_from_pt0(SALT,TEMP,c0)
    
    implicit none
    real (kind = 8) :: gsw_pd2_from_pt0
    real (kind = 8), intent(in) :: SALT, TEMP, c0
    real (kind = 8) :: CTEMP
    real (kind = 8), parameter :: pr = 2000.0

!    CTEMP = gsw_ct_from_pt(1e3_r8 * SALT,TEMP)
!    gsw_pd2_from_pt0 = 1e-3_r8 * gsw_rho(1e3_r8 * SALT,CTEMP,pr)
    CTEMP = gsw_ct_from_pt(c0 * SALT,TEMP)
    gsw_pd2_from_pt0 = 1e-3_r8 * gsw_rho(c0 * SALT,CTEMP,pr)
    end function gsw_pd2_from_pt0

end program isopycnalMOC

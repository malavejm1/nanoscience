subroutine parse_commands(file_line)
implicit none
character*255   :: file_line
character*255   :: the_key,value_string,temp_string
integer         :: i

    write(6,'(2a)') "INPUT::",trim(adjustl(file_line))
    ! Ignore anything on the line that appears after a #
    i=index(file_line,'#')
    if(i>0) then
      file_line=file_line(1:i-1)
    endif

    i=index(file_line,'=')
    if(i/=0) then
      read(file_line(1:i-1),'(a)')the_key
      read(file_line(i+1:),'(a)')value_string
      call convert_to_lowercase(the_key)
	  
	  select case(trim(adjustl(the_key)))
	  

            case("set_initial_fstate")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                set_initial_fstate=.TRUE.
              else
                set_initial_fstate=.FALSE.
              end if
              write(6,*) "SET:set_initial_fstate=",set_initial_fstate
              write(3,*) "SET:set_initial_fstate=",set_initial_fstate
            case("set_initial_f_e0")
               read(value_string,*) set_initial_f_e0
               write(6,*) "SET:set_initial_f_e0=",set_initial_f_e0
               write(3,*) "SET:set_initial_f_e0=",set_initial_f_e0
            case("set_initial_f_wid")
               read(value_string,*) set_initial_f_wid
               write(6,*) "SET:set_initial_f_wid=",set_initial_f_wid
               write(3,*) "SET:set_initial_f_wid=",set_initial_f_wid
            case("set_initial_f_shift")
               read(value_string,*) set_initial_f_shift
               write(6,*) "SET:set_initial_f_shift=",set_initial_f_shift
               write(3,*) "SET:set_initial_f_shift=",set_initial_f_shift
            case("guess_pib_states")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                guess_pib_states=.TRUE.
              else
                guess_pib_states=.FALSE.
              end if
              write(6,*) "SET:guess_pib_states=",guess_pib_states
              write(3,*) "SET:guess_pib_states=",guess_pib_states
            case("save_dump_files")
               read(value_string,*) save_dump_files
               write(6,*) "SET:save_dump_files=",save_dump_files
               write(3,*) "SET:save_dump_files=",save_dump_files
            case("switch_mix_cf_red_freq")
               read(value_string,*) switch_mix_cf_red_freq
               write(6,*) "SET:switch_mix_cf_red_freq=",switch_mix_cf_red_freq
               write(3,*) "SET:switch_mix_cf_red_freq=",switch_mix_cf_red_freq
            case("kick_psi")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                kick_psi=.TRUE.
              else
                kick_psi=.FALSE.
              end if
              write(6,*) "SET:kick_psi=",kick_psi
              write(3,*) "SET:kick_psi=",kick_psi
            case("kick_direction")
               read(value_string,*) kick_direction
               write(6,*) "SET:kick_direction=",kick_direction
               write(3,*) "SET:kick_direction=",kick_direction
            case("kick_plane")
               read(value_string,*) kick_plane
               write(6,*) "SET:kick_plane=",kick_plane
               write(3,*) "SET:kick_plane=",kick_plane
            case("kick_wid")
               read(value_string,*) kick_wid
               write(6,*) "SET:kick_wid=",kick_wid
               write(3,*) "SET:kick_wid=",kick_wid
            case("kick_k")
               read(value_string,*) kick_k
               write(6,*) "SET:kick_k=",kick_k
               write(3,*) "SET:kick_k=",kick_k
            case("kick_x0")
               read(value_string,*) kick_x0
               write(6,*) "SET:kick_x0=",kick_x0
               write(3,*) "SET:kick_x0=",kick_x0
            case("kick_y0")
               read(value_string,*) kick_y0
               write(6,*) "SET:kick_y0=",kick_y0
               write(3,*) "SET:kick_y0=",kick_y0
            case("kick_z0")
               read(value_string,*) kick_z0
               write(6,*) "SET:kick_z0=",kick_z0
               write(3,*) "SET:kick_z0=",kick_z0
            case("periodic")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                periodic=.TRUE.
              else
                periodic=.FALSE.
              end if
              write(6,*) "SET:periodic=",periodic
              write(3,*) "SET:periodic=",periodic
            case("nper_x")
               read(value_string,*) nper_x
               write(6,*) "SET:nper_x=",nper_x
               write(3,*) "SET:nper_x=",nper_x
            case("nper_y")
               read(value_string,*) nper_y
               write(6,*) "SET:nper_y=",nper_y
               write(3,*) "SET:nper_y=",nper_y
            case("nper_z")
               read(value_string,*) nper_z
               write(6,*) "SET:nper_z=",nper_z
               write(3,*) "SET:nper_z=",nper_z
            case("orbital_free")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                orbital_free=.TRUE.
              else
                orbital_free=.FALSE.
              end if
              write(6,*) "SET:orbital_free=",orbital_free
              write(3,*) "SET:orbital_free=",orbital_free
            case("tddft")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                tddft=.TRUE.
              else
                tddft=.FALSE.
              end if
              write(6,*) "SET:tddft=",tddft
              write(3,*) "SET:tddft=",tddft
            case("lambda_w")
               read(value_string,*) lambda_w
               write(6,*) "SET:lambda_w=",lambda_w
               write(3,*) "SET:lambda_w=",lambda_w
            case("lambda_tf")
               read(value_string,*) lambda_tf
               write(6,*) "SET:lambda_tf=",lambda_tf
               write(3,*) "SET:lambda_tf=",lambda_tf
            case("cdft")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                cdft=.TRUE.
              else
                cdft=.FALSE.
              end if
              write(6,*) "SET:cdft=",cdft
              write(3,*) "SET:cdft=",cdft
            case("addpot")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                addpot=.TRUE.
              else
                addpot=.FALSE.
              end if
              write(6,*) "SET:addpot=",addpot
              write(3,*) "SET:addpot=",addpot
            case("mgp")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                mgp=.TRUE.
              else
                mgp=.FALSE.
              end if
              write(6,*) "SET:mgp=",mgp
              write(3,*) "SET:mgp=",mgp
            case("rho0")
               read(value_string,*) rho0
               write(6,*) "SET:rho0=",rho0
               write(3,*) "SET:rho0=",rho0
            case("mgp_a")
               read(value_string,*) mgp_a
               write(6,*) "SET:mgp_a=",mgp_a
               write(3,*) "SET:mgp_a=",mgp_a
            case("mgp_b")
               read(value_string,*) mgp_b
               write(6,*) "SET:mgp_b=",mgp_b
               write(3,*) "SET:mgp_b=",mgp_b
            case("dipole_pert")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                dipole_pert=.TRUE.
              else
                dipole_pert=.FALSE.
              end if
              write(6,*) "SET:dipole_pert=",dipole_pert
              write(3,*) "SET:dipole_pert=",dipole_pert
            case("kick_strength")
               read(value_string,*) kick_strength
               write(6,*) "SET:kick_strength=",kick_strength
               write(3,*) "SET:kick_strength=",kick_strength
            case("kick_dir")
               read(value_string,*) kick_dir
               write(6,*) "SET:kick_dir=",kick_dir
               write(3,*) "SET:kick_dir=",kick_dir
            case("dkef")
               read(value_string,*) dkef
               write(6,*) "SET:dkef=",dkef
               write(3,*) "SET:dkef=",dkef
            case("n_electrons")
               read(value_string,*) n_electrons
               write(6,*) "SET:n_electrons=",n_electrons
               write(3,*) "SET:n_electrons=",n_electrons
            case("n_electrons_of")
               read(value_string,*) n_electrons_of
               write(6,*) "SET:n_electrons_of=",n_electrons_of
               write(3,*) "SET:n_electrons_of=",n_electrons_of
            case("n_x")
               read(value_string,*) n_x
               write(6,*) "SET:n_x=",n_x
               write(3,*) "SET:n_x=",n_x
            case("n_y")
               read(value_string,*) n_y
               write(6,*) "SET:n_y=",n_y
               write(3,*) "SET:n_y=",n_y
            case("n_z")
               read(value_string,*) n_z
               write(6,*) "SET:n_z=",n_z
               write(3,*) "SET:n_z=",n_z
            case("grid_step")
               read(value_string,*) grid_step
               write(6,*) "SET:grid_step=",grid_step
               write(3,*) "SET:grid_step=",grid_step
            case("n_cap_qm_x")
               read(value_string,*) n_cap_qm_x
               write(6,*) "SET:n_cap_qm_x=",n_cap_qm_x
               write(3,*) "SET:n_cap_qm_x=",n_cap_qm_x
            case("n_cap_qm_y")
               read(value_string,*) n_cap_qm_y
               write(6,*) "SET:n_cap_qm_y=",n_cap_qm_y
               write(3,*) "SET:n_cap_qm_y=",n_cap_qm_y
            case("n_cap_qm_z")
               read(value_string,*) n_cap_qm_z
               write(6,*) "SET:n_cap_qm_z=",n_cap_qm_z
               write(3,*) "SET:n_cap_qm_z=",n_cap_qm_z
            case("n_cap_rs_x")
               read(value_string,*) n_cap_rs_x
               write(6,*) "SET:n_cap_rs_x=",n_cap_rs_x
               write(3,*) "SET:n_cap_rs_x=",n_cap_rs_x
            case("n_cap_rs_y")
               read(value_string,*) n_cap_rs_y
               write(6,*) "SET:n_cap_rs_y=",n_cap_rs_y
               write(3,*) "SET:n_cap_rs_y=",n_cap_rs_y
            case("n_cap_rs_z")
               read(value_string,*) n_cap_rs_z
               write(6,*) "SET:n_cap_rs_z=",n_cap_rs_z
               write(3,*) "SET:n_cap_rs_z=",n_cap_rs_z
            case("load_gs")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                load_gs=.TRUE.
              else
                load_gs=.FALSE.
              end if
              write(6,*) "SET:load_gs=",load_gs
              write(3,*) "SET:load_gs=",load_gs
            case("num_gs_iter")
               read(value_string,*) num_gs_iter
               write(6,*) "SET:num_gs_iter=",num_gs_iter
               write(3,*) "SET:num_gs_iter=",num_gs_iter
            case("n_iteration_cg")
               read(value_string,*) n_iteration_cg
               write(6,*) "SET:n_iteration_cg=",n_iteration_cg
               write(3,*) "SET:n_iteration_cg=",n_iteration_cg
            case("n_iteration_switch")
               read(value_string,*) n_iteration_switch
               write(6,*) "SET:n_iteration_switch=",n_iteration_switch
               write(3,*) "SET:n_iteration_switch=",n_iteration_switch
            case("energy_convergence_tol")
               read(value_string,*) energy_convergence_tol
               write(6,*) "SET:energy_convergence_tol=",energy_convergence_tol
               write(3,*) "SET:energy_convergence_tol=",energy_convergence_tol
            case("switch_mix_coeff")
               read(value_string,*) switch_mix_coeff
               write(6,*) "SET:switch_mix_coeff=",switch_mix_coeff
               write(3,*) "SET:switch_mix_coeff=",switch_mix_coeff
            case("time_step")
               read(value_string,*) time_step
               write(6,*) "SET:time_step=",time_step
               write(3,*) "SET:time_step=",time_step
            case("n_time")
               read(value_string,*) n_time
               write(6,*) "SET:n_time=",n_time
               write(3,*) "SET:n_time=",n_time
            case("doutfreq")
               read(value_string,*) doutfreq
               write(6,*) "SET:doutfreq=",doutfreq
               write(3,*) "SET:doutfreq=",doutfreq
            case("e_0")
               read(value_string,*) e_0
               write(6,*) "SET:e_0=",e_0
               write(3,*) "SET:e_0=",e_0
            case("laser_wavelength")
               read(value_string,*) laser_wavelength
               write(6,*) "SET:laser_wavelength=",laser_wavelength
               write(3,*) "SET:laser_wavelength=",laser_wavelength
            case("laser_pulse_shift1")
               read(value_string,*) laser_pulse_shift1
               write(6,*) "SET:laser_pulse_shift1=",laser_pulse_shift1
               write(3,*) "SET:laser_pulse_shift1=",laser_pulse_shift1
            case("laser_pulse_width1")
               read(value_string,*) laser_pulse_width1
               write(6,*) "SET:laser_pulse_width1=",laser_pulse_width1
               write(3,*) "SET:laser_pulse_width1=",laser_pulse_width1
            case("ramp_time")
               read(value_string,*) ramp_time
               write(6,*) "SET:ramp_time=",ramp_time
               write(3,*) "SET:ramp_time=",ramp_time
            case("use_custom_laser")
               read(value_string,*) use_custom_laser
               write(6,*) "SET:use_custom_laser=",use_custom_laser
               write(3,*) "SET:use_custom_laser=",use_custom_laser
            case("custom_laser_type")
               read(value_string,*) custom_laser_type
               write(6,*) "SET:custom_laser_type=",custom_laser_type
               write(3,*) "SET:custom_laser_type=",custom_laser_type
            case("custom_laser_parm1")
               read(value_string,*) custom_laser_parm1
               write(6,*) "SET:custom_laser_parm1=",custom_laser_parm1
               write(3,*) "SET:custom_laser_parm1=",custom_laser_parm1
            case("custom_laser_parm2")
               read(value_string,*) custom_laser_parm2
               write(6,*) "SET:custom_laser_parm2=",custom_laser_parm2
               write(3,*) "SET:custom_laser_parm2=",custom_laser_parm2
            case("custom_laser_parm3")
               read(value_string,*) custom_laser_parm3
               write(6,*) "SET:custom_laser_parm3=",custom_laser_parm3
               write(3,*) "SET:custom_laser_parm3=",custom_laser_parm3
            case("custom_laser_parm4")
               read(value_string,*) custom_laser_parm4
               write(6,*) "SET:custom_laser_parm4=",custom_laser_parm4
               write(3,*) "SET:custom_laser_parm4=",custom_laser_parm4
            case("custom_laser_parm5")
               read(value_string,*) custom_laser_parm5
               write(6,*) "SET:custom_laser_parm5=",custom_laser_parm5
               write(3,*) "SET:custom_laser_parm5=",custom_laser_parm5
            case("custom_laser_parm6")
               read(value_string,*) custom_laser_parm6
               write(6,*) "SET:custom_laser_parm6=",custom_laser_parm6
               write(3,*) "SET:custom_laser_parm6=",custom_laser_parm6
            case("td_mode")
               read(value_string,*) td_mode
               write(6,*) "SET:td_mode=",td_mode
               write(3,*) "SET:td_mode=",td_mode
            case("rs_iter_per_ts")
               read(value_string,*) rs_iter_per_ts
               write(6,*) "SET:rs_iter_per_ts=",rs_iter_per_ts
               write(3,*) "SET:rs_iter_per_ts=",rs_iter_per_ts
            case("rs_helmhmode")
               read(value_string,*) rs_helmhmode
               write(6,*) "SET:rs_helmhmode=",rs_helmhmode
               write(3,*) "SET:rs_helmhmode=",rs_helmhmode
            case("scale_current")
               read(value_string,*) scale_current
               write(6,*) "SET:scale_current=",scale_current
               write(3,*) "SET:scale_current=",scale_current
            case("rs_freq_save_wall")
               read(value_string,*) rs_freq_save_wall
               write(6,*) "SET:rs_freq_save_wall=",rs_freq_save_wall
               write(3,*) "SET:rs_freq_save_wall=",rs_freq_save_wall
            case("time_order_mode")
               read(value_string,*) time_order_mode
               write(6,*) "SET:time_order_mode=",time_order_mode
               write(3,*) "SET:time_order_mode=",time_order_mode
            case("field_output_mode")
               read(value_string,*) field_output_mode
               write(6,*) "SET:field_output_mode=",field_output_mode
               write(3,*) "SET:field_output_mode=",field_output_mode
            case("moldyn")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                moldyn=.TRUE.
              else
                moldyn=.FALSE.
              end if
              write(6,*) "SET:moldyn=",moldyn
              write(3,*) "SET:moldyn=",moldyn
            case("atomic_movement_frequency")
               read(value_string,*) atomic_movement_frequency
               write(6,*) "SET:atomic_movement_frequency=",atomic_movement_frequency
               write(3,*) "SET:atomic_movement_frequency=",atomic_movement_frequency
            case("trajectory_output_frequency")
               read(value_string,*) trajectory_output_frequency
               write(6,*) "SET:trajectory_output_frequency=",trajectory_output_frequency
               write(3,*) "SET:trajectory_output_frequency=",trajectory_output_frequency
            case("calc_flux_z0")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                calc_flux_z0=.TRUE.
              else
                calc_flux_z0=.FALSE.
              end if
              write(6,*) "SET:calc_flux_z0=",calc_flux_z0
              write(3,*) "SET:calc_flux_z0=",calc_flux_z0
            case("use_core_conductivity")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                use_core_conductivity=.TRUE.
              else
                use_core_conductivity=.FALSE.
              end if
              write(6,*) "SET:use_core_conductivity=",use_core_conductivity
              write(3,*) "SET:use_core_conductivity=",use_core_conductivity
            case("read_in_switch_mix_coeff")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                read_in_switch_mix_coeff=.TRUE.
              else
                read_in_switch_mix_coeff=.FALSE.
              end if
              write(6,*) "SET:read_in_switch_mix_coeff=",read_in_switch_mix_coeff
              write(3,*) "SET:read_in_switch_mix_coeff=",read_in_switch_mix_coeff
            case("save_curr")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_curr=.TRUE.
              else
                save_curr=.FALSE.
              end if
              write(6,*) "SET:save_curr=",save_curr
              write(3,*) "SET:save_curr=",save_curr
            case("save_density")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_density=.TRUE.
              else
                save_density=.FALSE.
              end if
              write(6,*) "SET:save_density=",save_density
              write(3,*) "SET:save_density=",save_density
            case("save_den_slices")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_den_slices=.TRUE.
              else
                save_den_slices=.FALSE.
              end if
              write(6,*) "SET:save_den_slices=",save_den_slices
              write(3,*) "SET:save_den_slices=",save_den_slices
            case("save_e_induced")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_e_induced=.TRUE.
              else
                save_e_induced=.FALSE.
              end if
              write(6,*) "SET:save_e_induced=",save_e_induced
              write(3,*) "SET:save_e_induced=",save_e_induced
            case("save_e_full")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_e_full=.TRUE.
              else
                save_e_full=.FALSE.
              end if
              write(6,*) "SET:save_e_full=",save_e_full
              write(3,*) "SET:save_e_full=",save_e_full
            case("save_e_merged")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_e_merged=.TRUE.
              else
                save_e_merged=.FALSE.
              end if
              write(6,*) "SET:save_e_merged=",save_e_merged
              write(3,*) "SET:save_e_merged=",save_e_merged
            case("save_e_wv")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_e_wv=.TRUE.
              else
                save_e_wv=.FALSE.
              end if
              write(6,*) "SET:save_e_wv=",save_e_wv
              write(3,*) "SET:save_e_wv=",save_e_wv
            case("save_a_pot")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_a_pot=.TRUE.
              else
                save_a_pot=.FALSE.
              end if
              write(6,*) "SET:save_a_pot=",save_a_pot
              write(3,*) "SET:save_a_pot=",save_a_pot
            case("save_b_field")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                save_b_field=.TRUE.
              else
                save_b_field=.FALSE.
              end if
              write(6,*) "SET:save_b_field=",save_b_field
              write(3,*) "SET:save_b_field=",save_b_field
            case("output_stride")
               read(value_string,*) output_stride
               write(6,*) "SET:output_stride=",output_stride
               write(3,*) "SET:output_stride=",output_stride
            case("atom_crd_filename")
               write(atom_crd_filename,'(a)') trim(adjustl(value_string))
               write(6,*) "SET:atom_crd_filename=",trim(adjustl(atom_crd_filename))
               write(3,*) "SET:atom_crd_filename=",trim(adjustl(atom_crd_filename))
            case("load_atom_velocities")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                load_atom_velocities=.TRUE.
              else
                load_atom_velocities=.FALSE.
              end if
              write(6,*) "SET:load_atom_velocities=",load_atom_velocities
              write(3,*) "SET:load_atom_velocities=",load_atom_velocities
            case("conducting_shape")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                conducting_shape=.TRUE.
              else
                conducting_shape=.FALSE.
              end if
              write(6,*) "SET:conducting_shape=",conducting_shape
              write(3,*) "SET:conducting_shape=",conducting_shape
	  case("mix_conductor_schrodinger")
              call convert_to_lowercase(value_string)
              if((trim(adjustl(value_string))=="t").OR.(trim(adjustl(value_string))=="true").OR.(trim(adjustl(value_string))==".true.")) then
                mix_conductor_schrodinger=.TRUE.
              else
                mix_conductor_schrodinger=.FALSE.
              end if
              write(6,*) "SET:mix_conductor_schrodinger=",mix_conductor_schrodinger
              write(3,*) "SET:mix_conductor_schrodinger=",mix_conductor_schrodinger

		  
       case default
          write(6,'(2a)')"ERROR: Invalid variable name: ",trim(adjustl(the_key))
          stop		   
       end select	  
    end if
end subroutine parse_commands

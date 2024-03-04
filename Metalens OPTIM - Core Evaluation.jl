
################## MAIN DEFINITION #####################################################################################################################
#
function SM_main(disks_thickness, phase_shift, buffer_smooth, w0, focal_point, r_lens,gamma_prime,laser_detuning)
	#
	w0_x = w0_y = w0
	#
	#Lattice creation
	(r_atoms, n_atoms) = lattice_creation(r_lens, focal_point, disks_thickness,buffer_smooth, phase_shift)
	n_atoms>130000 ? error("Too many atoms: ", n_atoms) : nothing
	#
	#Core part (Green's function construction and inversion)
	E_field_in  = gaussian_beam_plus.(r_atoms[:,1],  r_atoms[:,2],  r_atoms[:,3],  w0_x, w0_y, k0)
	E_field_in  = [dot(dipoles_polarization,field_polarization) for index in 1:n_atoms].*E_field_in
	state_coeff = SM_inversion(r_atoms, n_atoms,dipoles_polarization, E_field_in,gamma_prime,laser_detuning)
	#
	#Initializes field arrays
	zR_in = k0*w0^2/2
	z0_in = 0.0
	magn  = focal_point/sqrt((z0_in-focal_point)^2+zR_in^2)
	w0_f  = w0*magn
	z0_f  = focal_point+(z0_in-focal_point)*magn^2
	E_field_focus     = gaussian_beam_plus.(r_atoms[:,1],  r_atoms[:,2],  r_atoms[:,3].-z0_f,  w0_f, w0_f, k0)
    atoms_mult        = [2*Int(r_atoms[i,1]!=0.0)+2*Int(r_atoms[i,2]!=0.0) for i in 1:n_atoms]
	#
	#Output field evaluation
	#(t_in, r_in, t_f)  = 
	t_f = SM_core_t(w0, conj.(E_field_focus), w0_f, z0_f, state_coeff,atoms_mult)
	#
	#return (t_in, r_in, t_f, n_atoms)
	return abs2(t_f)
end

################## SAVING FUNCTIONS DEFINITION #########################################################################################################
#
#Functions to write data files
function h5write_complex(file_name,data)
    file_h5=h5open(file_name*".h5", "w")
    file_h5["real"]=real.(data)
    file_h5["imag"]=imag.(data)
    close(file_h5)
end
function h5write_multiple(file_name,data_array)
    file_h5=h5open(file_name*".h5", "w")
    for index in 1:length(data_array)
		file_h5[(data_array[index])[1]]=(data_array[index])[2]
    end
    close(file_h5)
end
#

################## LATTICE FUNCTIONS DEFINITIONS #######################################################################################################
#
#Functions to construct the metalens
function nearest_choice(array, value)
    findmin(abs.(array.-value))[2]
end
#
function lattice_creation(r_lens_start, focal_point, disks_thickness, buffer_zone, phase_shift)
    phase_full_data=h5read("Data In/finalDataPhaseNOISE"*".h5", "/PhaseData")
	phase_full_data[4,:]=shift_phase.(phase_full_data[4,:])
	#
	phi_func = (r1,r2) -> mod(k0*(focal_point-sqrt(focal_point^2+((r1+r2)/2)^2))+phase_shift,2*pi)
	#
	if fill_until_r_lens_option
		r_lens = ceil(r_lens_start/disks_thickness)*disks_thickness
	else
		r_lens = r_lens_start
	end
	#
	r_range = collect(0.0:disks_thickness:r_lens)
	#
	n_phase_disks = length(r_range)-1
	buffer_range  = vcat(0.0,r_range[2:end-1].+(r_range[2]*buffer_zone))
	if phase_center_ring_option=="NO"
    	phase_range   = phi_func.(buffer_range, r_range[2:end])
	else
		phase_range   = phi_func.(r_range[1:end-1], r_range[2:end])
	end
    r_atoms       = Array{Float64}(undef,0,3)
	r_atoms_old   = Array{Float64}(undef,0,3)
	r_atoms_new   = Array{Float64}(undef,0,3)
    phase_array   = Array{Float64}(undef,n_phase_disks)
	lattice_array = Array{Float64}(undef,n_phase_disks,3)
	max_lattice_const = 1.0
	old_x_max = 0.0
	old_y_max = 0.0
    for i in 1:n_phase_disks
        current_phase_index = nearest_choice(phase_full_data[4,:], phase_range[i])
        phase_array[i]      = phase_full_data[4,   current_phase_index]
        lattice_array[i,:]  = phase_full_data[1:3, current_phase_index]
		if (old_x_max, old_y_max) == (0.0,0.0)
			old_x_max = lattice_array[i,1]/2
			old_y_max = lattice_array[i,2]/2
		end
		if lattice_array[i,1]<max_lattice_const && lattice_array[i,2]<max_lattice_const
			(r_atoms_new, old_x_max, old_y_max) = lattice_creation_core(r_lens, lattice_array[i,:],r_range[i+1],r_range[i], old_x_max, old_y_max,focal_point)
			if i>1 && length(r_atoms_old[:,1])>0
				a_x_variation = abs(lattice_array[i,1]-lattice_array[i-1,1])/lattice_array[i,1]
				a_y_variation = abs(lattice_array[i,2]-lattice_array[i-1,2])/lattice_array[i,2]
				if !(a_x_variation>0.01 && a_y_variation>0.01) && buffer_zone>0.0
					(r_atoms_new, old_x_max, old_y_max) = lattice_creation_core(r_lens, lattice_array[i,:],r_range[i+1],buffer_range[i], old_x_max, old_y_max,focal_point)
					if a_x_variation > a_y_variation
						indices_order=[1;2]
					else
						indices_order=[2;1]
					end
					(up_nodes_1_start, up_nodes_2_start)   = get_lower(r_atoms_new[:,indices_order[1]],r_atoms_new[:,indices_order[2]])
					(low_nodes_1, low_nodes_2) = get_upper(r_atoms_old[:,indices_order[1]],r_atoms_old[:,indices_order[2]])
					max_low = maximum(low_nodes_1)+lattice_array[i-1,indices_order[1]]*0.75
					select_up = up_nodes_1_start.<max_low
					up_nodes_1 = up_nodes_1_start[select_up]
					up_nodes_2 = up_nodes_2_start[select_up]
					select_up_axis = (up_nodes_1_start.>max_low).*(up_nodes_2_start.>lattice_array[i,indices_order[2]])
					up_nodes_1_to_axis = up_nodes_1_start[select_up_axis]
					up_nodes_2_to_axis = up_nodes_2_start[select_up_axis]
					if length(up_nodes_1)<=length(low_nodes_1)
						(up_nodes_1_buff, up_nodes_2_buff,low_nodes_1_buff, low_nodes_2_buff) = (up_nodes_1, up_nodes_2, low_nodes_1, low_nodes_2)
						index_up_down = 1
					else
						(up_nodes_1_buff, up_nodes_2_buff,low_nodes_1_buff, low_nodes_2_buff) = (low_nodes_1, low_nodes_2, up_nodes_1, up_nodes_2)
						index_up_down = 2
					end
					(r_atoms_buffer_x, r_atoms_buffer_y) = (lattice_creation_buffer(up_nodes_1_buff, up_nodes_2_buff,low_nodes_1_buff, low_nodes_2_buff,lattice_array[i,indices_order[2]] ,index_up_down))[indices_order]
					#
					(r_atoms_buffer_x_axis, r_atoms_buffer_y_axis) = (lattice_creation_buffer(up_nodes_1_to_axis, up_nodes_2_to_axis,up_nodes_1_to_axis, [-lattice_array[i,indices_order[2]]/2 for ii in 1:length(up_nodes_2_to_axis)],lattice_array[i,indices_order[2]],1))[indices_order]
					#
					r_atoms_buffer      = lattice_creation_buffer_z(r_atoms_buffer_x, r_atoms_buffer_y, lattice_array[i,3],focal_point)#!!!!!!!!!!!!!!!!!!!!!
					r_atoms = vcat(r_atoms,r_atoms_buffer)
					if z_fixed_option=="YES"
						r_atoms_buffer_axis = lattice_creation_buffer_z(r_atoms_buffer_x_axis, r_atoms_buffer_y_axis, lattice_array[i,3],focal_point)#!!!!!!!!!!!!!!!!!
					else
						r_atoms_buffer_axis = lattice_creation_buffer_z(r_atoms_buffer_x_axis, r_atoms_buffer_y_axis, 0.0,focal_point)
					end
					r_atoms = vcat(r_atoms,r_atoms_buffer_axis)
				else
					(r_atoms_new, old_x_max, old_y_max) = lattice_creation_core(r_lens, lattice_array[i,:],r_range[i+1],r_range[i], old_x_max, old_y_max,focal_point)
				end
			end
			if length(r_atoms_new)>0
				r_atoms = vcat(r_atoms,r_atoms_new)
				r_atoms_old=r_atoms_new[:,:]
			else
				r_atoms_old=Array{Float64}(undef,0,3)
				(old_x_max, old_y_max) == (0.0,0.0)
			end
		else
			r_atoms_old=Array{Float64}(undef,0,3)
			(old_x_max, old_y_max) == (0.0,0.0)
		end
    end
	#
	#
	if fill_until_r_lens_option
		selected_atoms = ((xx,yy)->sqrt(xx^2+yy^2)<=r_lens_start).(r_atoms[:,1],r_atoms[:,2])
		r_atoms = r_atoms[selected_atoms,:]
	end
	#
	#
    return (r_atoms[:,:],length(r_atoms[:,1]))
end
#
function shift_phase(phi)
	if phi<=0.0
		return phi+2*pi
	else
		return phi
	end
end
#
#Functions to construct a lattice
function boundaries_choice(nAtoms)
    if nAtoms%2==0
        return [-nAtoms/2; nAtoms/2-1]
    else
        return [-(nAtoms-1)/2; (nAtoms-1)/2]
    end
end
#
function lattice_creation_core(r_lens, lattice_constants, r_max,r_min, x_shift, y_shift,focal_point)
    a_x=lattice_constants[1] 
    a_y=lattice_constants[2] 
    a_z=lattice_constants[3] 
    naX = Int(floor(2*r_lens/a_x))+1
    naY = Int(floor(2*r_lens/a_y))+1
	naX%2!=0 ? naX+=1 : nothing
	naY%2!=0 ? naY+=1 : nothing
    boundX=boundaries_choice(naX)
    boundY=boundaries_choice(naY)
    xOption=(naX+1)%2
    yOption=(naY+1)%2
    lattice_array_x=[(i + xOption/2)*a_x for i in boundX[1]:boundX[2] for j in boundY[1]:boundY[2]]
    lattice_array_y=[(j + yOption/2)*a_y for i in boundX[1]:boundX[2] for j in boundY[1]:boundY[2]]
	#if shift_to_match_option=="YES"
	#	lattice_array_x = lattice_array_x.+(x_shift-xOption/2*a_x)
	#	lattice_array_y = lattice_array_y.+(y_shift-yOption/2*a_y)
	#end
    lattice_norms=((xx,yy)->sqrt(xx^2+yy^2)).(lattice_array_x,lattice_array_y)
    selected_points=(lattice_norms.<r_max).*(lattice_norms.>=r_min).*(lattice_array_x.>=0.0).*(lattice_array_y.>=0.0)
    lattice_array_x=lattice_array_x[selected_points]
    lattice_array_y=lattice_array_y[selected_points]
    n_points_selected=length(lattice_array_x)
	#
	if n_points_selected>0
		max_x = maximum(lattice_array_x)
		max_y = maximum(lattice_array_y)
	else
		max_x = 0.0
		max_y = 0.0
	end
	#
	if z_fixed_option=="YES"
		return (hcat(repeat(lattice_array_x,3),repeat(lattice_array_y,3),[z*a_z for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ]), max_x,max_y )
	else
		phi_func(rr) = mod(k0*(focal_point-sqrt(focal_point^2+rr^2))+phase_shift,2*pi)
		z_func(rr)   = (2*pi-mod(phi_func(rr),pi) )/(6*pi)
		return ( hcat(repeat(lattice_array_x,3),repeat(lattice_array_y,3),[z*(z_func(sqrt(lattice_array_x[i]^2+lattice_array_y[i]^2))) for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ]) , max_x,max_y )
	end
	#
	#
end
#
#
function get_upper(pos_array_1,pos_array_2)
	index_sorted=sort(collect(1:length(pos_array_1)), by=(xx->pos_array_1[xx]))
	pos_array_x=pos_array_1[index_sorted]
	pos_array_y=pos_array_2[index_sorted]
	current_max=pos_array_y[1]
	current_x=pos_array_x[1]
	up_x=[current_x]
	up_y=[current_max]
	current_index=1
	for i in 2:length(pos_array_x)
		x_here = pos_array_x[i]
		y_here = pos_array_y[i]
		if x_here>current_x
			up_x=vcat(up_x,x_here)
			up_y=vcat(up_y,y_here)
			current_x=x_here
			current_max = y_here
			current_index+=1
			#println("STEP: ", i, ",   x = ", x_here,",   y = ", y_here , ". I add a new column.")
		else
			if y_here>current_max
				up_x[current_index]=x_here
				up_y[current_index]=y_here
				current_max = y_here
			end
		end
	end
	return (up_x,up_y)
end

function get_lower(pos_array_1,pos_array_2)
	index_sorted=sort(collect(1:length(pos_array_1)), by=(xx->pos_array_1[xx]))
	pos_array_x=pos_array_1[index_sorted]
	pos_array_y=pos_array_2[index_sorted]
	current_min=pos_array_y[1]
	current_x=pos_array_x[1]
	down_x=[current_x]
	down_y=[current_min]
	current_index=1
	for i in 2:length(pos_array_x)
		x_here = pos_array_x[i]
		y_here = pos_array_y[i]
		if x_here>current_x
			down_x=vcat(down_x,x_here)
			down_y=vcat(down_y,y_here)
			current_x=x_here
			current_min = y_here
			current_index+=1
		else
			if y_here<current_min
				down_x[current_index]=x_here
				down_y[current_index]=y_here
				current_min = y_here
			end
		end
	end
	return (down_x,down_y)
end


function lattice_creation_buffer(up_nodes_1, up_nodes_2,low_nodes_1, low_nodes_2, a_2, index_up_down)
	index_low = 1
	index_low_max = length(low_nodes_2)
	r_atoms_buffer_x = Array{Float64}(undef,0)
	r_atoms_buffer_y = Array{Float64}(undef,0)
	for i in 1:length(up_nodes_1)
		if index_low<=index_low_max
			p1_x=up_nodes_1[i]
			p1_y=up_nodes_2[i]
			p2_x=low_nodes_1[index_low]
			p2_y=low_nodes_2[index_low]
			final_index_low = index_low
			for j in (index_low+1):index_low_max
				p2_x_test = low_nodes_1[j]
				p2_y_test = low_nodes_2[j]
				if abs(p1_x-p2_x_test)<abs(p1_x-p2_x)
					final_index_low = j
					p2_x = p2_x_test
					p2_y = p2_y_test
				end
			end
			if (p1_y>p2_y,p1_y<p2_y )[index_up_down]
				index_low = final_index_low
				(x_points_to_add, y_points_to_add) = line_connect(p1_x, p1_y,p2_x,p2_y, a_2)
				r_atoms_buffer_x = vcat(r_atoms_buffer_x,x_points_to_add)
				r_atoms_buffer_y = vcat(r_atoms_buffer_y,y_points_to_add)
				index_low+=1
			end
		end
	end
	return (r_atoms_buffer_x, r_atoms_buffer_y)
end

function line_connect(p1_x, p1_y,p2_x,p2_y, a_2)
	#println("connect [",p1_x,",",p1_y,"] with [", p2_x,",",p2_y,"]" )
	y_min = min(p1_y, p2_y)+a_2*0.9
	y_max = max(p1_y, p2_y)-a_2*0.9
	y_points = collect(y_min:a_2:y_max)
	x_points = @. (p1_x*y_points - p2_x*y_points + p2_x*p1_y - p1_x*p2_y)/(p1_y - p2_y)
	return (x_points,y_points)
end

function lattice_creation_buffer_z(r_atoms_buffer_x, r_atoms_buffer_y, a_z_sharp,focal_point)
	n_points_selected = length(r_atoms_buffer_x)
	if a_z_sharp>0.0
		return hcat(repeat(r_atoms_buffer_x,3),repeat(r_atoms_buffer_y,3),[z*a_z_sharp for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ])
	else
		phi_func(rr) = mod(k0*(focal_point-sqrt(focal_point^2+rr^2))+phase_shift,2*pi)
		z_func(rr)   = (2*pi-mod(phi_func(rr),pi) )/(6*pi)
		return hcat(repeat(r_atoms_buffer_x,3),repeat(r_atoms_buffer_y,3),[z*(z_func(sqrt(r_atoms_buffer_x[i]^2+r_atoms_buffer_y[i]^2))) for z in [-1 ; 0 ; 1 ].+0.0 for i in 1:n_points_selected ])
	end
end



function debug_r_atoms(r_lens , focal_point, disks_thickness,buffer_smooth, phase_shift)
    (r_atoms,n_atoms) = lattice_creation(r_lens , focal_point, disks_thickness,buffer_smooth, phase_shift) 
    if fill_until_r_lens_option
        h5write_multiple("r_atoms_YESfill_"*name_add,[("r_atoms", r_atoms)])
    else
        h5write_multiple("r_atoms_NOfill_"*name_add,[("r_atoms", r_atoms)])
    end
    error("CHECK COMPLETED. CORRECT EXIT.s")
end
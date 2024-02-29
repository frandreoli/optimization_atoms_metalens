
################## SM CORE EVALUATION FUNCTIONS DEFINITIONS ############################################################################################
#
#Function to invert the effective Hamiltonian
function SM_inversion(r_atoms, n_atoms,dipoles_polarization,E_field_in)
	SM_green_matrix  =  Array{Complex{TN},2}(undef, n_atoms, n_atoms)
	time_temp=time()
	SM_initialize!(SM_green_matrix, r_atoms, dipoles_polarization)
    #println("Green's matrix initialized in         ", time()-time_temp)
	time_temp=time()
	state_coeff  = SM_green_matrix\E_field_in
	#println("SM core evaluation finished in        ", time()-time_temp)
	return state_coeff
end
#
#Function to evaluate the SM transmission coefficient
#function SM_core_t(E_field_in, w0_in,E_field_f_conj, w0_f, z0_f, state_coeff,atoms_mult) 
function SM_core_t(w0_in,E_field_f_conj, w0_f, z0_f, state_coeff,atoms_mult)
	#E_field_in_conj=conj.(E_field_in)
	t_alpha_term=(im*(3/(4*pi^2))*((lambda0^2)/(w0_in^2)))
    #t_in    = 1 + (t_alpha_term*sum(E_field_in_conj.*state_coeff.*atoms_mult))
	#r_in    =     (t_alpha_term*sum(E_field_in.*     state_coeff.*atoms_mult))
	coeff_f        = exp(2.0im*pi*z0_f/lambda0)*(2*pi*w0_in*w0_f)/(pi*(w0_in^2+w0_f^2)+1.0im*z0_f*lambda0)
	t_alpha_term_f = t_alpha_term*(w0_in/w0_f)
	t_f            = coeff_f + (t_alpha_term_f*sum(E_field_f_conj.*state_coeff.*atoms_mult))
    #return (t_in, r_in, t_f)
	return t_f
end
#

################## INITIALIZATION FUNCTIONS DEFINITIONS ################################################################################################
#
#Function to fill the green tensor for the SM
function SM_initialize!(SM_green_matrix, r_vecs, p)
    na = length(r_vecs[:,1])
	r_vecs_x = r_vecs[:,1].*k0
	r_vecs_y = r_vecs[:,2].*k0
	r_vecs_z = r_vecs[:,3].*k0
	n_steps = na^2#Int((na^2 + na)/2)
	Threads.@threads for index in  1:n_steps
		#i = Int(floor( -(1/2) + sqrt(1/4 + 2*(index - 1 ))   )) + 1
		#j = Int(index - (i - 1)*i/2)
		i = Int(floor(  (index-1)/na  )) + 1
		j = Int( index-(i-1)*na)
		x_i=r_vecs_x[i]
		y_i=r_vecs_y[i]
		x_j=r_vecs_x[j]
		y_j=r_vecs_y[j]
		z = r_vecs_z[i]- r_vecs_z[j]
		to_add=0.0+0.0im
		#First step
		if i==j
			to_add+=-(0.5im*(1+gamma_prime)+laser_detuning) 
		else
			x = x_i-x_j
			y = y_i-y_j
			r =sqrt((x^2)+(y^2)+(z^2))
			cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
			to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
		end
		#Second step
		if x_j!=0.0
			x = x_i+x_j
			y = y_i-y_j
			r =sqrt((x^2)+(y^2)+(z^2))
			cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
			to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
		end
		#Third step
		if y_j!=0.0
			x = x_i-x_j
			y = y_i+y_j
			r =sqrt((x^2)+(y^2)+(z^2))
			cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
			to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
		end
		#Fourth step
		if x_j!=0.0 && y_j!=0.0
			x = x_i+x_j
			y = y_i+y_j
			r =sqrt((x^2)+(y^2)+(z^2))
			cos_theta_square = ((x*p[1])+(y*p[2])+(p[3]*z))^2
			to_add+=-(3*exp(im*r)/(4*r^3)*((r^2 + im*r - 1)*1+(3 - r^2 - 3im*r)*cos_theta_square/r^2))
		end
		#End
		SM_green_matrix[i, j] = Complex{TN}(to_add) #SM_green_matrix[i, j] = SM_green_matrix[j, i] = to_add
	end
end
#


################## PHYSICAL FUNCTIONS DEFINITIONS ######################################################################################################
#
# Functions to calculate beam field at position r_vec, where z_vec is the propagation direction z and the beam has waists w0_x and w0_y in the x and y directions
function zr(w0, lambda)
    pi*w0^2/lambda
end
function w(z, w0, lambda)
    w0*sqrt(1 + (z/zr(w0, lambda))^2) # beam width w(z)
end
function invr(z, w0, lambda)
    z/(z^2 + zr(w0, lambda)^2) # inverse of Radius of curvature of beam R(z)
end
function gaussian_beam_plus(x_vec,y_vec, z_vec, w0_x, w0_y, k)
    #This function wants r_vec NOT in units of lambda_0!
    lambda = 2*pi/k
    w_x = w(z_vec, w0_x, lambda)
    w_y = w(z_vec, w0_y, lambda)
    norm_factor = sqrt(w0_x*w0_y/(w_x*w_y))
    real_arg    = -(x_vec^2/w_x^2 + y_vec^2/w_y^2)
    imag_arg    = k*z_vec + k*(x_vec^2*invr(z_vec, w0_x, lambda) + y_vec^2*invr(z_vec, w0_y, lambda))/2 - atan(z_vec/zr(w0_x, lambda))/2 - atan(z_vec/zr(w0_y, lambda))/2
    return norm_factor*exp(real_arg + im*imag_arg)
end
#

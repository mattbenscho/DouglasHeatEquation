##### VARIABLE SETUP #####

roomtemp = 300;
watertemp = 290;
crystallength = 1.5e-2; 
endcaplength = 2e-3;
crystalwidth = 3e-3;
samples = 30;
global P;
load pumppower.dat P;
P
nddoping = 0.004; 
cyvo4 = 530;
rho = 4.22 * 1e-3 * 1e6;
kyvo4 = 2.34e-6;
kair = 21.6e-6;                   
kcopper = 115e-6;               
rhoA = 1.247e28;
lambdap = 808e-9;
lambdal = 1064e-9;
clight = 3e8;
abscrosss = 2.3e-23;
planck = 6.6e-34; # in J*s
pumpbeamwidth = 840e-6
# load timestep.dat timestep;
timestep = 1e-4
# load threshold.dat threshold;
threshold = 1
load eta.dat eta;
eta

##### DERIVED VALUES SETUP #####

domainwidth = 2 * crystalwidth;
domainheight = domainwidth;
domainlength = 2 * crystallength;
global h = domainwidth / samples;
h
DeltaV = h^3;
capi = int32( endcaplength / h );
crystalfront = ceil( (domainlength / 2. - crystallength / 2.)/h );
crystalback = floor( (domainlength / 2. + crystallength / 2.)/h );
crystalleft = ceil( (domainwidth / 2. - crystalwidth / 2.)/h );
crystalright = ceil( (domainwidth / 2. + crystalwidth / 2.)/h );
crystaltop = crystalright;
crystalbottom = crystalleft;
tcount = 10 * timestep;
global r = timestep / h^2;
r
global M = floor(domainwidth / h);
global L = floor(domainheight / h);
global J = floor(domainlength / h);
crystalmiddle = int32(L / 2.);
Un = roomtemp * ones(M,L,J);
Unplusonealt = roomtemp;
NNd = nddoping * rhoA * DeltaV;
ndexcited = zeros(M,L,J);
ndrelaxed = zeros(M,L,J);
nup = clight / lambdap;
nul = clight / lambdal;
global K = kair * ones(M,L,J)
global Kmatrixx
global Kmatrixy
global Kmatrixz
global tridiagonalmatrixx
global tridiagonalmatrixy
global tridiagonalmatrixz
global gamma = 2 / (pi * pumpbeamwidth^2);

##### FUNCTION DEFINITIONS #####

function f = pumplight(x, y, mux, muy, wp)
  global gamma
  global P
  f = P * gamma * exp( -2 * ( (x-mux)^2 + (y-muy)^2 ) ./ wp^2 );
endfunction

function setupmatrices
  global r
  global K
  global Kmatrixx
  global Kmatrixy
  global Kmatrixz
  global tridiagonalmatrixx
  global tridiagonalmatrixy
  global tridiagonalmatrixz
  global M L J

  printf("Setting up tridiagonalmatrixx... \n");
  for l = 1 : L
    for j = 1 : J
      Kmatrixx{l,j} = zeros(M-2,M-2);
      for i = 1 : M-3
	Kmatrixx{l,j}(i,i)   = -(K(i+1,l,j)+K(i,l,j));
	Kmatrixx{l,j}(i+1,i) = K(i+1,l,j);
	Kmatrixx{l,j}(i,i+1) = K(i+1,l,j);
      endfor
      Kmatrixx{l,j}(M-2,M-2) = -(K(M-1,l,j)+K(M-2,l,j));
      tridiagonalmatrixx{l,j} = eye(M-2) - r/2. * Kmatrixx{l,j};
    endfor
  endfor
  
  printf("Setting up tridiagonalmatrixy... \n");
  for m = 1 : M
    for j = 1 : J
      Kmatrixy{m,j} = zeros(L-2,L-2);
      for i = 1 : L-3
	Kmatrixy{m,j}(i,i)   = -(K(m,i+1,j)+K(m,i,j));
	Kmatrixy{m,j}(i+1,i) = K(m,i+1,j);
	Kmatrixy{m,j}(i,i+1) = K(m,i+1,j);
      endfor
      Kmatrixy{m,j}(L-2,L-2) = -(K(m,L-1,j)+K(m,L-2,j));
      tridiagonalmatrixy{m,j} = eye(L-2) - r/2. * Kmatrixy{m,j};
    endfor
  endfor
  
  printf("Setting up tridiagonalmatrixz... \n");
  for m = 1 : M
    for l = 1 : L
      Kmatrixz{m,l} = zeros(J-2,J-2);
      for i = 1 : J-3
	Kmatrixz{m,l}(i,i)   = -(K(m,l,i+1)+K(m,l,i));
	Kmatrixz{m,l}(i+1,i) = K(m,l,i+1);
	Kmatrixz{m,l}(i,i+1) = K(m,l,i+1);
      endfor
      Kmatrixz{m,l}(J-2,J-2) = -(K(m,l,J-1)+K(m,l,J-2));
      tridiagonalmatrixz{m,l} = eye(J-2) - r/2. * Kmatrixz{m,l};
    endfor
  endfor
endfunction

# Create K-Matrix
for j = crystalfront : crystalback
  K(:,:,j) = kcopper * ones(M,L);
endfor
for j = crystalfront : crystalback
  for m = crystalleft : crystalright
    for l = crystalbottom : crystaltop
      K(m,l,j) = kyvo4;
    endfor
  endfor
endfor

# Setup boundary conditions
for j = crystalfront : crystalback
  Un(1,:,j) = watertemp * ones(L,1);
  Un(M,:,j) = watertemp * ones(L,1);
  Un(:,1,j) = watertemp * ones(M,1);
  Un(:,L,j) = watertemp * ones(M,1);
endfor

Unplusone = Ustrichstrich = Ustrich = Un;
load (sprintf ("unplusone-3-eta=%01i-P=%02i.dat", eta, P), "Unplusone");

setupmatrices;

# Setup relaxed Nd atoms
for j = crystalfront + capi : crystalback - capi
  for m = crystalleft : crystalright
    for l = crystalbottom : crystaltop
      ndrelaxed(m,l,j) = NNd;
    endfor
  endfor
endfor

############# ITERATION START #############

timeelapsed = 0;

# datafile = sprintf("data/3D/%02i-%09i.dat", P, 10^6*timeelapsed);
# fid = fopen (datafile, "w");
# for m = 1 : M
#   for j = 1 : J
#     fprintf(fid,"%i %i %f\n",m,j,Unplusone(m,crystalmiddle,j));
#   endfor
#   fprintf(fid,"\n");
# endfor
# fclose (fid);

while (timeelapsed < 20)

  Un = Unplusone;
  Unplusonealt = ...
      Unplusone(crystalmiddle,crystalmiddle,crystalfront+capi);

  #### start heat deposition ####

  # First, fill the reservoir
  pumppower = 0;
  for m = crystalleft : crystalright
    for l = crystalbottom : crystaltop
      pumpphotons(m,l) = ...
          ( h^2 * timestep * ...
	    pumplight(h*m, h*l, h*M/2, h*L/2, pumpbeamwidth) ) ...
	  / ( planck * nup );
      pumppower += pumpphotons(m,l) * 1/timestep * planck * nup;
    endfor
  endfor

  # Distribute the energy from the reservoir
  for j = crystalfront + capi : crystalback - capi
    for m = crystalleft : crystalright
      for l = crystalbottom : crystaltop
  	if ( pumpphotons(m,l) > 0 )
  	  absphotons = h * pumpphotons(m,l) * abscrosss ...
  	      * (ndrelaxed(m,l,j)/NNd) * rhoA;
  	  if ( absphotons > pumpphotons(m,l) )
  	    ndexcited(m,l,j) = ndexcited(m,l,j) + pumpphotons(m,l);
  	    pumpphotons(m,l) = 0;
  	  else
  	    ndexcited(m,l,j) = ndexcited(m,l,j) + absphotons;
  	    pumpphotons(m,l) = pumpphotons(m,l) - absphotons;
  	  endif
  	endif
      endfor
    endfor
  endfor

  # relaxation to ground state generates energy:
  for j = crystalfront + capi : crystalback - capi
    for m = crystalleft : crystalright
      for l = crystalbottom : crystaltop
  	Un(m,l,j) = Un(m,l,j) ...
       	    + eta * ( ndexcited(m,l,j) * planck * (nup-nul) ) ...
       	      / ( cyvo4 * rho * DeltaV );
	ndexcited(m,l,j) = 0;
  	ndrelaxed(m,l,j) = NNd - ndexcited(m,l,j);
      endfor
    endfor
  endfor

  # calculate absorbed power
  leftpumppower = 0;
  for m = crystalleft : crystalright
    for l = crystalbottom : crystaltop
      leftpumppower += pumpphotons(m,l) * 1/timestep * planck * nup;
    endfor
  endfor
  abspumppower = pumppower - leftpumppower;

  #### end heat deposition ####  

  #### start heat transfer ####

  ### find the elements for the douglas equations ###
  dx2Un(2:M-1,2:L-1,2:J-1) = ...
      K(2:M-1,2:L-1,2:J-1) .* Un(3:M,2:L-1,2:J-1) - ... 
      ( K(2:M-1,2:L-1,2:J-1) + K(1:M-2,2:L-1,2:J-1) ) .* ...
        Un(2:M-1,2:L-1,2:J-1) + ...
      K(1:M-2,2:L-1,2:J-1) .* Un(1:M-2,2:L-1,2:J-1);
  dy2Un(2:M-1,2:L-1,2:J-1) = ...
      K(2:M-1,2:L-1,2:J-1) .* Un(2:M-1,3:L,2:J-1) - ...
      ( K(2:M-1,2:L-1,2:J-1) + K(2:M-1,1:L-2,2:J-1) ) .* ...
        Un(2:M-1,2:L-1,2:J-1) + K(2:M-1,1:L-2,2:J-1) .* ...
      Un(2:M-1,1:L-2,2:J-1);
  dz2Un(2:M-1,2:L-1,2:J-1) = ...
      K(2:M-1,2:L-1,2:J-1) .* Un(2:M-1,2:L-1,3:J) -...
      ( K(2:M-1,2:L-1,2:J-1) + K(2:M-1,2:L-1,1:J-2) ) .* ...
        Un(2:M-1,2:L-1,2:J-1) + K(2:M-1,2:L-1,1:J-2) .* ...
      Un(2:M-1,2:L-1,1:J-2);

  ### boundary elements ###
  dx2Un(2,  2:L-1,2:J-1) += K(1,  2:L-1,2:J-1) .* Un(1,2:L-1,2:J-1);
  dx2Un(M-1,2:L-1,2:J-1) += K(M-1,2:L-1,2:J-1) .* Un(M,2:L-1,2:J-1);
  dy2Un(2:M-1,2,  2:J-1) += K(2:M-1,1,  2:J-1) .* Un(2:M-1,1,2:J-1);
  dy2Un(2:M-1,L-1,2:J-1) += K(2:M-1,L-1,2:J-1) .* Un(2:M-1,L,2:J-1);
  dz2Un(2:M-1,2:L-1,2  ) += K(2:M-1,2:L-1,1  ) .* Un(2:M-1,2:L-1,1);
  dz2Un(2:M-1,2:L-1,J-1) += K(2:M-1,2:L-1,J-1) .* Un(2:M-1,2:L-1,J);

  for j = 2 : J-1 
    for l = 2 : L-1 # solve the first douglas equation
      Ustrich(2:M-1,l,j) = tridiagonalmatrixx{l,j} \ ...
	  ( Un(2:M-1,l,j) + r/2. * dx2Un(2:M-1,l,j) + ...
	   r * dy2Un(2:M-1,l,j) + r * dz2Un(2:M-1,l,j));
    endfor
    for m = 2 : M-1 # solve the second douglas equation
      Ustrichstrich(m,2:L-1,j) = tridiagonalmatrixy{m,j} \ ...
	  transpose( Ustrich(m,2:L-1,j) - r/2. * dy2Un(m,2:L-1,j));
    endfor
  endfor
 
  # solve the third douglas equation
  for m = 2 : M-1
    for l = 2 : L-1
      Unplusone(m,l,2:J-1) = tridiagonalmatrixz{m,l} \ ...
  	  (squeeze(Ustrichstrich(m,l,2:J-1)) ...
           - r/2. * squeeze(dz2Un(m,l,2:J-1)));
    endfor
  endfor

  #### end heat transfer ###

  tcount = tcount + timestep;
  timeelapsed = timeelapsed + timestep;

  #### start output ###

  printf(".");

  if (tcount >= 20 * timestep)
    tcount = 0;
    DeltaTrate = (- Unplusonealt ...
		  + Unplusone(crystalmiddle, crystalmiddle,
			      crystalfront + capi)) / timestep;
    printf(" DT=%5.3f K/s ", DeltaTrate);
    save (sprintf("unplusone-3-eta=%01i-P=%02i.dat",eta,P),"Unplusone");

    if ( sqrt(DeltaTrate^2) < threshold)
      printf("\n");
      printf("Steady state reached.\n");
      printf("Absorbed pump power: %6.4f W.\n", abspumppower);
      exit
    endif

  endif

  #### end output ####

endwhile

exit
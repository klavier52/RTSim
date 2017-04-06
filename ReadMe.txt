MIT Rocket Team- Simulation Environment
Project Start 4/5/17

Guidelines for editing
	- Project lives in git. First create a branch from main. Work on your branch and commit freely.
	- After any work, use the git add function, then commit and COMMENT on what you did
	- When have something complete and working, pull from main, merge changes, then push to main
	- Push to main if certain it won't break stuff. TEST
		We wont be unit testing much because of pace so we rely on good code
	- If there is no comment on a commit, you are fired
	- Add your name to a ToDo item you are working on. Then you OWN it!

ToDo Tier 1: Essencial items to reach next goal
	- Schedule weekly meetings
	- Schedule work plan- get goals from exec
	- Come up with plan for class structure (aka the integrator should be a class and so should dynamics)
	- Test integrator by comparing to analyticly solvable functions
	- Cart3D aerosim to get table of coefficients
		Step file of rocket body
	- RasAero aerosim to get table of coefficients
	- Compare Cart3D results to RasAero
	- Add viscous effects and skin friction drag? Ask Darmofal for best method given timeline
	- Come up with dynamics equations plans (should we 3dof then 6dof?)
		Density model
		Gravity model
		Allow input of aero coefficients
		Moment of Interia calc
	- Add comment to each function with specifications on input/output
		Add flags if unanticipated input/output
		Unit tests?
		How do we test the full 6dof? Compare to Open Rocket results?
	- Sensitiviy analysis on all variables
	- Better definition of main while loop that includes flags on apogee etc...
	- Look at when order 4 minus order 5 calc is 0- right now set to e-16
	- How to easily integrate monte carlo capabilities
		Make class that all vars go through
		Each var needs a mean and standard dev. or other probabolistic model
		How to deal with table skews (result can be wanky/discont. if done for each table var by step)

ToDo Tier 2: Cool ideas to have in the backburner
	- Automatic ballast adjuster to get mean of monte carlo at target altitude
	- Move integrator to FORTRAN? Call from python for speed
	- Wind payload- before launch send up payload to measure wind profile along altitude
		Upload file to monte-carlo. Adjust ballast to get mean at target altitude
	- Legendre Gaussian integrator?
	- Add events such as launch, coast, recovery, landing
	- Add anaylis tools that grab data like landing spot distance, vehicle vel after x seconds beyond apogee (for parachute force purposes), apogee, 'wobble' freq
	
How it works right now
 Run from main.py
	Setup located here: Defines properties and tables such as thrust curve and air density
	Functions to be included in dynamics selected here. Include each variable you want recorded
	Runs the integrator (Rungga Kutta Fehlsburg) for one time step
		Has variable step size, may run through loop multiple times before stepping forward
	Integrator calls specified dynamics functions as appropiate
	Outputs dynamics
	









